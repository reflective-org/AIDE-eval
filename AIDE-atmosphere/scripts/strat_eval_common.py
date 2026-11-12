"""
Shared constants, sources, readers and math of the stratospheric-dynamics evaluation
(the transport scorecard): 01_fields.py, 02_scorecard.py, 03_figures.py and
checks/port_check.py.

The scorecard is the jcm-strat relative scorecard (phase 29, phase 31's scored set)
ported into this repository: age of air (A), residual circulation (B), polar vortex (C)
and variability (D), each metric scored as its relative deviation from a reference.
docs/EVALUATION_PROTOCOL.md is the description; this module holds the parts the scripts
share, so that they cannot disagree.

It merges transport_common.py with the parts of aide_val_common.py the scorecard uses
(CESM constants, the h6 loader, interp_level); the rest of aide_val_common.py belonged to
the two-tier protocol, archived as git tag tier-protocol-final.

CESM h6 conventions, established against the files (evidence: scripts 01-01d at tag
tier-protocol-final):

  1. VTHzm is already the eddy flux v'theta'; do not subtract Vzm*THzm.
  2. Wzm is the log-pressure vertical velocity w = -H omega / p, H = 7 km.
  3. Below-surface points are the sentinel 1e35 with no _FillValue; masked with abs < 1e20.
  4. MSKtem is not a 0/1 mask (a fractional count of above-surface interfaces) and is
     not applied to the zonal-mean fields.
"""

from __future__ import annotations

import glob
import json
import os
import warnings
import zipfile

import netCDF4
import numpy as np
import xarray as xr

# --------------------------------------------------------------------- constants
# CESM values (CAM shr_const), used for model and references alike
A_EARTH = 6.37122e6      # m, CAM shr_const_rearth
G0 = 9.80616             # m/s2, CAM gravit
RGAS = 287.058           # J/kg/K, CAM rair
CP = 1004.64             # J/kg/K, CAM cpair
KAPPA = RGAS / CP        # 0.2858
H_SCALE = 7000.0         # m, standard log-pressure scale height (AHL 1987)
P_REF_HPA = 1000.0       # reference pressure for theta and log-pressure height
FILL_SENTINEL = 1e20     # anything larger is the 1e35 below-surface sentinel

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTDIR = os.path.join(_ROOT, "output")

DAYS_PER_YEAR = 365.25          # the jcm-strat clocks count days; the contract asks years

# --------------------------------------------------------------------- references
# Every reference the scorecard reads, and the years taken from each. Paths under
# /data stay there; the ERA5 files are copies made by scripts/copy_era5_refs.sh
# and checked against scripts/era5_refs.sha256.
DATA = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")
REFS = dict(
    # CLaMS v3.1 driven by ERA5: zonal-mean AGE [yr] on (month, press, lat), one
    # member per year - A2, A3, A4 (with the entry offset) and D7
    clams_zip="/data/CLaMS/CLaMS_v3/clams_v3.1_era5_zm_lat.zip",
    clams_member="clams_v3.1_era5_zm_lat_press_{year}.nc",
    # WACCM6 REF-D1 entry age [yr], base point 0.47N 103 hPa - A1
    waccm_aoa="/data/CESM2_REFD1_AOA/AoA_waccm6_refd1.04_AOA1mf_1970-2019_ba_0_100.0_ck_0_50.0.nc",
    # CESM2.1.5-WACCM6 histSST h6 tape (Vzm, VTHzm, THzm) through load_h6 - B1, B2
    waccm_tem_segment="1996-2014",
    # ERA5 monthly zonal means of u and T, 25 levels - C1, C2
    era5_zm_monthly=os.path.join(DATA, "era5", "era5_zm_monthly_{year}.nc"),
    # ERA5 CDS daily zonal means u, v, T, theta, v'theta', 17 levels 1-250 hPa - D1-D7
    era5_daily=os.path.join(DATA, "era5", "era5_{year}.npz"),
)
# Reference years follow the scored window (each result set's `window`): the CLaMS and
# WACCM6 age references (A), ERA5 (C, D) and CLaMS (D7) are taken over the window's
# own years. The WACCM6 TEM reference (B) is the whole 1996-2014 h6 segment, the
# only one, whatever the window.
WACCM_TEM_YEARS = (1996, 2014)   # all of the 1996-2014 h6 segment
WACCM_ILEV = slice(0, 61)        # h6 interfaces from the model top to 652 hPa: Psi* is
                                 # integrated from the top, so the top is kept


# --------------------------------------------------------------------- CESM h6 tape
# CESM2.1.5-WACCM6 fixed-SST histSST, daily zonal-mean TEM tape (cam.h6): Uzm Vzm Wzm
# THzm VTHzm UVzm UWzm on the 71 hybrid interfaces (ilev), noleap. The scorecard reads
# the 1996-2014 run (group B reference). Above ilev index 53 (~182 hPa) hybi = 0, so the
# interface pressure is hyai*P0 and needs no surface pressure.
H6_ROOT = "/data/cesm2.1.5_output/histSST"
H6_SEGMENTS = {
    "1996-2014": dict(
        case="f.e21.FWHIST.f09_f09_mg17.atmos-scale_fixedSST_1996-2014.001",
        subdir="archive/atm/proc/tseries/day_1",
        span="19960101-20141231"),
}


def h6_path(segment, var):
    s = H6_SEGMENTS[segment]
    return os.path.join(H6_ROOT, s["case"], s["subdir"], f"{s['case']}.cam.h6.{var}.{s['span']}.nc")


def load_h6(segment, var, ilev_slice):
    """One zonal-mean h6 variable as (time, ilev, lat), the 1e35 sentinel masked."""
    ds = xr.open_dataset(h6_path(segment, var), decode_times=True)
    da = ds[var].isel(zlon=0, ilev=ilev_slice)
    da = da.where(np.abs(da) < FILL_SENTINEL)
    da = da.load()
    ds.close()
    return da


def interp_level(field, p_hpa, p_target, axis=1):
    """Linear-in-log(p) interpolation of `field` onto p_target; p_hpa ascending."""
    lp = np.log(p_hpa)
    lt = np.log(p_target)
    k = int(np.searchsorted(lp, lt)) - 1
    k = max(0, min(k, len(lp) - 2))
    w = (lt - lp[k]) / (lp[k + 1] - lp[k])
    lo = np.take(field, k, axis=axis)
    hi = np.take(field, k + 1, axis=axis)
    return (1.0 - w) * lo + w * hi


# The field contract (docs/EVALUATION_PROTOCOL.md): required, then optional.
CONTRACT = ("Uzm", "Vzm", "THzm", "VTHzm", "AOA150")
OPTIONAL = ("T", "Wzm")
UNITS = {"Uzm": "m/s", "Vzm": "m/s", "THzm": "K", "VTHzm": "K m/s",
         "AOA150": "yr", "T": "K", "Wzm": "m/s"}

# --------------------------------------------------------------------- result sets
# One JSON file per result set in AIDE-atmosphere/sources/<set>.json: the runs scored
# together over one window, one scorecard each. A new model is one file with one run.
#
#   {"description": "...",                    what the set is
#    "reader": "jcm_npz" | "netcdf",          READERS below
#    "window": [Y0, Y1],                      scored years, inclusive
#    "runs": {"NAME": "path" | ["path", ...]},  a list is read as one continuous record
#    "notes": {"NAME": "..."}}                optional, one line per run
SOURCES_DIR = os.path.join(_ROOT, "sources")
SET_KEYS = ("description", "reader", "window", "runs")


def set_names():
    """Every registered result set, sorted."""
    return sorted(os.path.splitext(f)[0] for f in os.listdir(SOURCES_DIR) if f.endswith(".json"))


def load_set(name):
    """One result set from sources/<name>.json, checked against the schema above."""
    path = os.path.join(SOURCES_DIR, f"{name}.json")
    if not os.path.exists(path):
        raise SystemExit(f"no result set {name!r}: {path} missing; registered: {set_names()}")
    spec = json.load(open(path))
    missing = [k for k in SET_KEYS if k not in spec]
    if missing:
        raise SystemExit(f"{path}: missing {missing}")
    if spec["reader"] not in READERS:
        raise SystemExit(f"{path}: reader {spec['reader']!r} not in {list(READERS)}")
    y0, y1 = spec["window"]
    if not (isinstance(y0, int) and isinstance(y1, int) and y0 <= y1) or not spec["runs"]:
        raise SystemExit(f"{path}: window must be [Y0, Y1], runs non-empty")
    spec["runs"] = {n: ([p] if isinstance(p, str) else list(p)) for n, p in spec["runs"].items()}
    spec.setdefault("notes", {})
    spec["name"] = name
    return spec


def years_of(window):
    return list(range(window[0], window[1] + 1))


# jcm-strat npz names -> contract names. T is kept: C2 and D6 read it.
JCM_NAMES = {"u": "Uzm", "v": "Vzm", "th": "THzm", "vth": "VTHzm",
             "aoa150": "AOA150", "T": "T"}


# --------------------------------------------------------------------- readers
def orient(F):
    """Levels and latitudes ascending, as interp_level and band_mean assume."""
    fields = [k for k in CONTRACT + OPTIONAL if k in F]
    if F["plev"][0] > F["plev"][-1]:
        F["plev"] = F["plev"][::-1]
        F.update({k: F[k][:, ::-1, :] for k in fields})
    if F["lat"][0] > F["lat"][-1]:
        F["lat"] = F["lat"][::-1]
        F.update({k: F[k][:, :, ::-1] for k in fields})
    return F


def read_jcm_npz(paths, years):
    """jcm-strat daily_zm.py npz files: days, p_hpa, lat and (day, lev, lat) float32
    u, v, T, th, vth, aoa150 [, aoa_sfc]. The clock is in days; returned in years.
    Several files are one record: they must share the grid and join day to day, with
    no gap and no overlap."""
    parts = []
    for path in paths:
        z = np.load(path)
        days = z["days"].astype("datetime64[D]")
        keep = np.isin(days.astype("datetime64[Y]").astype(int) + 1970, list(years))
        F = {new: z[old][keep] for old, new in JCM_NAMES.items()}
        F.update(plev=z["p_hpa"].astype("f8"), lat=z["lat"].astype("f8"), days=days[keep])
        parts.append(F)
    F = parts[0]
    for path, G in zip(paths[1:], parts[1:]):
        if not (np.array_equal(F["plev"], G["plev"]) and np.array_equal(F["lat"], G["lat"])):
            raise SystemExit(f"{path}: grid differs from {paths[0]}")
        if G["days"].size and F["days"].size and G["days"][0] - F["days"][-1] != np.timedelta64(1, "D"):
            raise SystemExit(f"{path}: starts {G['days'][0]}, the record before ends "
                             f"{F['days'][-1]}: gap or overlap")
        F = {k: (np.concatenate([F[k], G[k]]) if k in JCM_NAMES.values() or k == "days" else F[k])
             for k in F}
    F["AOA150"] = F["AOA150"] / np.float32(DAYS_PER_YEAR)
    return orient(F)


def read_netcdf(paths, years):
    """Any model that writes the contract names: one or more netCDF files (each
    path may be a glob), dims (time, plev, lat), plev in hPa.

    Any calendar: each day is labelled by its own year, month and day, so a
    noleap record simply has no 29 February. Values beyond 1e20 are masked as
    fill (the 1e35 sentinel convention). T, if absent, is THzm (p/1000)^kappa.
    """
    files = sorted(f for p in paths for f in glob.glob(p))
    if not files:
        raise SystemExit(f"no files match {paths}")
    ds = xr.concat([xr.open_dataset(f) for f in files], dim="time")
    missing = [k for k in CONTRACT if k not in ds]
    if missing:
        raise SystemExit(f"{paths}: missing contract fields {missing}")
    yr, mo, dy = (ds.time.dt.year.values, ds.time.dt.month.values, ds.time.dt.day.values)
    keep = np.isin(yr, list(years))
    days = np.array([f"{y:04d}-{m:02d}-{d:02d}" for y, m, d in
                     zip(yr[keep], mo[keep], dy[keep])], dtype="datetime64[D]")
    F = {}
    for k in CONTRACT + OPTIONAL:
        if k in ds:
            x = ds[k].transpose("time", "plev", "lat").values[keep]
            F[k] = np.where(np.abs(x) < FILL_SENTINEL, x, np.nan)
    F.update(plev=ds.plev.values.astype("f8"), lat=ds.lat.values.astype("f8"), days=days)
    if "T" not in F:
        F["T"] = F["THzm"] * (F["plev"][None, :, None] / P_REF_HPA) ** KAPPA
    ds.close()
    return orient(F)


READERS = {"jcm_npz": read_jcm_npz, "netcdf": read_netcdf}


def load_run(spec, run):
    """The contract fields of one run of a result set over the set's window, levels and
    latitudes ascending. The record must cover the window from 1 January to 31 December."""
    w = spec["window"]
    F = READERS[spec["reader"]](spec["runs"][run], years_of(w))
    first, last = np.datetime64(f"{w[0]}-01-01"), np.datetime64(f"{w[1]}-12-31")
    if F["days"].size == 0 or F["days"][0] != first or F["days"][-1] != last:
        got = f"{F['days'][0]} .. {F['days'][-1]}" if F["days"].size else "no days"
        raise SystemExit(f"{spec['name']} / {run}: window {w[0]}-{w[1]} not covered ({got})")
    return F


# --------------------------------------------------------------------- scoring
# jcm-strat scripts/relative_scores.py: rel = model / ref - 1, err = |rel| / TOL.
# Metrics with no meaningful relative form keep their absolute bound and are
# reported as a difference in their own unit.
TOL = 0.25
A4_RANGE = (2.3, 3.2)            # yr; A4 is scored against the midpoint 2.75
ABS_BOUND = {"A2": 0.4, "A2-WACCM": 0.4, "C2-N": 2.0, "C2-S": 2.0, "D4-fwN": 10.0, "D4-fwS": 10.0,
             **{f"C1-jet{p}_{h}-lat": 4.0 for p in (50, 10) for h in "NS"}}
ABS_UNIT = {"A2": "yr", "A2-WACCM": "yr", "C2-N": "K", "C2-S": "K", "D4-fwN": "days", "D4-fwS": "days"}


# Two scorings of the same diagnostics, both scored as above (Susanne 2026-11-11,
# renamed 2026-11-12):
#   standard   the headline score and ranking: jcm-strat phase 31's scored set (the
#              former "reduced" scoring). A2-WACCM replaces A2, and the metrics of
#              make_phase31_report.py UNSCORED stay in the table but out of every group
#              score (winds at 10 hPa, vortex variability, final warmings, SSW count, and
#              the two relations D6, D7). Phase 31 scored them with absolute bounds
#              (success_metrics.py); here they get the relative scoring.
#   extensive  the full listing: the 34 metrics of the phase-29 scorecard (A2 against
#              shifted CLaMS); A2-WACCM is not part of it
# hidden: not shown; unscored: shown, not in the group scores.
PHASE31_UNSCORED = ("C1-jet10_N", "C1-jet10_N-lat", "C1-jet10_S", "C1-jet10_S-lat",
                    "D1-N", "D1-S", "D2-N", "D2-S", "D3", "D4-fwN", "D4-fwS",
                    "D6-N", "D6-S", "D7")
SCORINGS = {
    "standard": dict(hidden=(), unscored=("A2",) + PHASE31_UNSCORED),
    "extensive": dict(hidden=("A2-WACCM",), unscored=()),
}

# The standard window (Susanne 2026-11-12). A set scored on any other window is flagged
# in its reports and in both summaries.
STANDARD_WINDOW = (2006, 2019)


def window_flag(window):
    """'' on the standard window, else 'non-standard window (N of 14 standard years)'."""
    if tuple(window) == STANDARD_WINDOW:
        return ""
    n = len(set(years_of(window)) & set(years_of(STANDARD_WINDOW)))
    return f"non-standard window ({n} of {len(years_of(STANDARD_WINDOW))} standard years)"

def ref_value(ref):
    """A range reference (A4) is scored against its midpoint."""
    return float(np.mean(ref)) if isinstance(ref, (list, tuple)) else float(ref)


def score(mid, model, ref):
    """(value, err) of one metric. value is % off the reference for a relative
    metric, the difference in the metric's own unit for an absolute one (for a
    jet latitude |model| - |ref|). A metric that could not be formed fails (err 3)."""
    r = ref_value(ref)
    if not np.isfinite(model):
        return float("nan"), 3.0
    if mid in ABS_BOUND:
        d = abs(model) - abs(r) if mid.endswith("-lat") else model - r
        return float(d), float(abs(d) / ABS_BOUND[mid])
    rel = model / r - 1.0
    return float(100.0 * rel), float(abs(rel) / TOL)


def status(err):
    return "pass" if err <= 1.0 else ("near" if err <= 1.5 else "fail")


# --------------------------------------------------------------------- grid helpers
# Ported from jcm-strat sweep_score.py, which group A and the jcm-strat numbers use.
def nearest_level(p_hpa, target):
    """Index of the level nearest `target` in pressure (sweep_score.at_level)."""
    return int(np.argmin(np.abs(np.asarray(p_hpa) - target)))


def band_abs(lat, row, lo, hi):
    """cos(lat)-weighted mean over lo <= |lat| <= hi, both hemispheres, finite
    points only (sweep_score.band)."""
    m = (np.abs(lat) >= lo) & (np.abs(lat) <= hi) & np.isfinite(row)
    w = np.cos(np.deg2rad(lat[m]))
    return float(np.sum(row[m] * w) / np.sum(w))


def on_grid(field, p_src, lat_src, p_dst, lat_dst):
    """(p, lat) field onto (p_dst, lat_dst): linear in latitude, then in log p;
    NaN outside the source range (sweep_score.on_grid)."""
    o = np.argsort(p_src)
    lp, lpo = np.log(p_src[o]), np.log(np.asarray(p_dst, float))
    f = field[o]
    ol = np.argsort(lat_src)
    f, lat_s = f[:, ol], lat_src[ol]
    tmp = np.stack([np.interp(lat_dst, lat_s, f[k], left=np.nan, right=np.nan)
                    for k in range(f.shape[0])])
    return np.stack([np.interp(lpo, lp, tmp[:, j], left=np.nan, right=np.nan)
                     for j in range(tmp.shape[1])], axis=1)


def weighted_rmse(diff, lat, mask):
    """cos(lat)-weighted RMSE over the finite points of `mask` (sweep_score.weighted_rmse)."""
    w = np.cos(np.deg2rad(lat))[None, :] * np.ones_like(diff)
    ok = np.isfinite(diff) & mask
    return float(np.sqrt(np.sum((w * diff ** 2)[ok]) / np.sum(w[ok])))


def weighted_bias(diff, lat, mask):
    """cos(lat)-weighted mean over the finite points of `mask` (the bias that
    sweep_score.weighted_rmse returns beside the RMSE)."""
    w = np.cos(np.deg2rad(lat))[None, :] * np.ones_like(diff)
    ok = np.isfinite(diff) & mask
    return float(np.sum((w * diff)[ok]) / np.sum(w[ok]))


# --------------------------------------------------------------------- age references
def clams_member(year):
    """One CLaMS v3.1 year, read from the zip into memory (the extracted copies are
    not readable by every uid)."""
    with zipfile.ZipFile(REFS["clams_zip"]) as z:
        raw = z.read(REFS["clams_member"].format(year=year))
    return netCDF4.Dataset(f"clams_{year}", memory=raw)


def clams_age(years):
    """CLaMS mean AGE [yr] (press, lat): mean over months, then over years, -1e20
    fill excluded (aoa_vs_clams.clams_age)."""
    ages = []
    for y in years:
        with clams_member(y) as d:
            a = np.ma.filled(np.ma.asarray(d.variables["AGE"][:], dtype="f8"), np.nan)
            p = np.asarray(d.variables["press"][:], dtype="f8")
            lat = np.asarray(d.variables["lat"][:], dtype="f8")
        with warnings.catch_warnings():      # columns below the surface are all fill
            warnings.simplefilter("ignore", RuntimeWarning)
            ages.append(np.nanmean(np.where(a > -1e20, a, np.nan), axis=0))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return p, lat, np.nanmean(np.stack(ages), axis=0)


def waccm_age(years):
    """WACCM6 REF-D1 entry age [yr] (lev, lat), mean over the years' records,
    negative values (below the base point) masked (aoa_vs_clams.waccm_age).
    The hybrid lev is taken as pressure, exact above about 150 hPa."""
    with xr.open_dataset(REFS["waccm_aoa"], decode_times=False) as w:
        d = w.date.values
        sel = (d >= years[0] * 10000 + 101) & (d < (years[-1] + 1) * 10000 + 101)
        a = w.AOA.isel(time=np.where(sel)[0]).mean("time")
        a = a.where(a >= 0)
        return np.asarray(w.lev, "f8"), np.asarray(w.lat, "f8"), a.values.astype("f8")


# --------------------------------------------------------------------- residual circulation
# jcm-strat strat_circulation.py: the TEM residual streamfunction from time-mean
# [v], [v'theta'], [theta] in pressure coordinates, and w* from its latitude
# derivative - no w-bar needed. Constants are CESM's (above) where jcm-strat has
# a = 6.371e6 m, g = 9.80665 m/s2: w* scales with 1/a (-0.004 %), Psi* with a/g (+0.009 %).
def tem_streamfunction(vbar, vth, thbar, p_pa, lat_deg):
    """(p, lat) time means on pressure p_pa [Pa, any order] -> p ascending [Pa],
    v* [m/s] and Psi* [kg/s]:

        v*  = [v] - d/dp ( [v'theta'] / d[theta]/dp )
        Psi* = (2 pi a cos(phi) / g) * int_top^p v* dp'     (trapezoid from the top level)
    """
    order = np.argsort(p_pa)
    p = p_pa[order]
    vbar, vth, thbar = vbar[order], vth[order], thbar[order]
    dth_dp = np.gradient(thbar, p, axis=0)
    dth_dp = np.where(np.abs(dth_dp) < 1e-6, np.sign(dth_dp + 1e-30) * 1e-6, dth_dp)
    vstar = vbar - np.gradient(vth / dth_dp, p, axis=0)
    integrand = vstar * np.cos(np.deg2rad(lat_deg))[None, :]
    psi = np.zeros_like(vstar)
    psi[1:] = np.cumsum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(p)[:, None], axis=0)
    psi *= 2 * np.pi * A_EARTH / G0
    return p, vstar, psi


def residual_w(p, psi, lat_deg):
    """w* [mm/s] on (p, lat) from Psi*: omega* = -(g / (2 pi a^2 cos(phi))) dPsi*/dphi,
    w* = -H omega* / p, H = 7 km. Rows poleward of 88 deg are masked (1/cos)."""
    phi = np.deg2rad(lat_deg)
    dpsi = np.gradient(psi, phi, axis=1)
    omega = -(G0 / (2 * np.pi * A_EARTH ** 2 * np.cos(phi)[None, :])) * dpsi
    w = -H_SCALE * omega / p[:, None] * 1e3
    w[:, np.abs(lat_deg) > 88.0] = np.nan
    return w


def tropical_wstar(p_hpa, w, lat, target, nearest=False):
    """cos(lat)-weighted w* [mm/s] over |lat| <= 15 at `target` hPa: linear in
    log p between the bracketing levels (interp_level), or at the
    nearest level as jcm-strat does (nearest=True)."""
    row = (w[nearest_level(p_hpa, target)] if nearest
           else interp_level(w, p_hpa, target, axis=0))
    return band_abs(lat, row, 0, 15)


def upward_mass_flux(p_hpa, psi, lat, target, nearest=False):
    """Tropical upward mass flux [kg/s]: max - min of Psi* over |lat| <= 60 at
    `target` hPa, Psi* interpolated linearly in log p (or the nearest level)."""
    row = (psi[nearest_level(p_hpa, target)] if nearest
           else interp_level(psi, p_hpa, target, axis=0))
    row = row[np.abs(lat) <= 60.0]
    return float(row.max() - row.min())


def tropical_profile(w, lat, lat_max=15.0):
    """cos(lat)-weighted mean over |lat| <= lat_max at every level."""
    m = np.abs(lat) <= lat_max
    c = np.cos(np.deg2rad(lat[m]))
    return np.nansum(w[:, m] * c, axis=1) / np.sum(c)
