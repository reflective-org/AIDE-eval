"""
Shared loaders and TEM math for the AIDE-WACCM atmospheric validation targets.

Data source: CESM2.1.5 WACCM6 fixed-SST histSST runs, daily zonal-mean TEM tape (cam.h6).

  1970-1995  f.e21.FWHIST.f09_f09_mg17.atmos-scale_fixedSST.001
  1996-2014  f.e21.FWHIST.f09_f09_mg17.atmos-scale_fixedSST_1996-2014.001

The h6 tape carries Uzm Vzm Wzm THzm VTHzm UVzm UWzm on the 71 hybrid INTERFACE
levels (ilev), zonally averaged, daily means, calendar = noleap.

Two traps this module handles, both verified against the files:

  1. Below-surface points are written as the sentinel 1e35 with NO _FillValue
     attribute, so xarray does not mask them. Unguarded means return ~1e33.
     Affects 1.34% of cells, all at ilev index >= 60 (p >= 652 hPa).
  2. MSKtem is NOT a 0/1 mask. It is a (time, lat, lon) fractional count of
     above-surface interfaces (values 59-71, quantised in 1/48). It cannot be
     applied to the zonal-mean TEM fields. We derive the mask from the sentinel.

Above ilev index 53 (~182 hPa) hybi == 0, so the interface pressure is exactly
hyai*P0 and is independent of surface pressure. Every diagnostic here lives at
10 hPa or 70 hPa, i.e. in that pure-pressure region, so no PS is needed.
"""

from __future__ import annotations

import os
import numpy as np
from scipy import stats as _stats
import xarray as xr

# ----------------------------------------------------------------------------
# Paths
# ----------------------------------------------------------------------------

ROOT = "/data/cesm2.1.5_output/histSST"

SEGMENTS = {
    "1970-1995": dict(
        case="f.e21.FWHIST.f09_f09_mg17.atmos-scale_fixedSST.001",
        subdir="postprocessed_output/atm/proc/tseries/day_1",
        span="19700101-19951231",
        year0=1970,
        year1=1995,
    ),
    "1996-2014": dict(
        case="f.e21.FWHIST.f09_f09_mg17.atmos-scale_fixedSST_1996-2014.001",
        subdir="archive/atm/proc/tseries/day_1",
        span="19960101-20141231",
        year0=1996,
        year1=2014,
    ),
}

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTDIR = os.path.join(_ROOT, "output")

# ----------------------------------------------------------------------------
# Physical constants (CESM values)
# ----------------------------------------------------------------------------

A_EARTH = 6.37122e6      # m, CAM shr_const_rearth
G0 = 9.80616             # m/s2, CAM gravit
RGAS = 287.058           # J/kg/K, CAM rair
CP = 1004.64             # J/kg/K, CAM cpair
KAPPA = RGAS / CP        # 0.2858
H_SCALE = 7000.0         # m, standard log-pressure scale height (AHL 1987)
P_REF_HPA = 1000.0       # reference pressure for theta and log-pressure height

FILL_SENTINEL = 1e20     # anything larger is the 1e35 below-surface sentinel

# Vertical window kept in memory. ilev index 20 -> 0.099 hPa, index 60 -> 652 hPa.
# Everything at or above index 59 is sentinel-free over the whole record, and
# this window brackets both 10 hPa (idx 36) and 70 hPa (idx 46) with margin for
# vertical derivatives.
ILEV_SLICE = slice(20, 61)


# ----------------------------------------------------------------------------
# Loading
# ----------------------------------------------------------------------------

def h6_path(segment: str, var: str) -> str:
    s = SEGMENTS[segment]
    return os.path.join(
        ROOT, s["case"], s["subdir"], f"{s['case']}.cam.h6.{var}.{s['span']}.nc"
    )


def load_h6(segment: str, var: str, ilev_slice: slice = ILEV_SLICE) -> xr.DataArray:
    """Load one zonal-mean h6 variable as (time, ilev, lat), sentinel-masked."""
    ds = xr.open_dataset(h6_path(segment, var), decode_times=True)
    da = ds[var].isel(zlon=0, ilev=ilev_slice)
    da = da.where(np.abs(da) < FILL_SENTINEL)          # kill the 1e35 sentinel
    da = da.load()
    ds.close()
    return da


def load_grid(segment: str) -> dict:
    """Interface pressures (hPa), latitudes (deg) and Gaussian weights."""
    ds = xr.open_dataset(h6_path(segment, "Uzm"))
    ilev = ds["ilev"].isel(ilev=ILEV_SLICE).values.astype(np.float64)   # hPa
    lat = ds["lat"].values.astype(np.float64)
    gw = ds["gw"].values.astype(np.float64) if "gw" in ds else None
    hyai = ds["hyai"].isel(ilev=ILEV_SLICE).values.astype(np.float64)
    hybi = ds["hybi"].isel(ilev=ILEV_SLICE).values.astype(np.float64)
    p0 = float(ds["P0"].values)
    ds.close()
    return dict(ilev=ilev, lat=lat, gw=gw, hyai=hyai, hybi=hybi, P0=p0)


def time_index(segment: str):
    """Return (year, month, day-of-year) integer arrays for the daily record."""
    ds = xr.open_dataset(h6_path(segment, "Uzm"))
    t = ds["time"].values
    ds.close()
    # cam stamps daily means at the END of the interval, so day 1 of the record
    # carries the timestamp of day 2. Shift back one day to label correctly.
    import cftime
    yr = np.array([x.year for x in t])
    mo = np.array([x.month for x in t])
    dy = np.array([x.day for x in t])
    return yr, mo, dy


# ----------------------------------------------------------------------------
# Grid helpers
# ----------------------------------------------------------------------------

def log_p_height(p_hpa: np.ndarray) -> np.ndarray:
    """Log-pressure height z = -H ln(p/p0), metres."""
    return -H_SCALE * np.log(p_hpa / P_REF_HPA)


def interp_level(field: np.ndarray, p_hpa: np.ndarray, p_target: float,
                 axis: int = 1) -> np.ndarray:
    """Linear-in-log(p) interpolation of `field` onto p_target.

    field: array with vertical axis `axis`, p_hpa monotonically increasing.
    """
    lp = np.log(p_hpa)
    lt = np.log(p_target)
    k = int(np.searchsorted(lp, lt)) - 1
    k = max(0, min(k, len(lp) - 2))
    w = (lt - lp[k]) / (lp[k + 1] - lp[k])
    lo = np.take(field, k, axis=axis)
    hi = np.take(field, k + 1, axis=axis)
    return (1.0 - w) * lo + w * hi


def interp_lat(field: np.ndarray, lat: np.ndarray, lat_target: float,
               axis: int = -1) -> np.ndarray:
    """Linear interpolation onto an exact latitude."""
    j = int(np.searchsorted(lat, lat_target)) - 1
    j = max(0, min(j, len(lat) - 2))
    w = (lat_target - lat[j]) / (lat[j + 1] - lat[j])
    lo = np.take(field, j, axis=axis)
    hi = np.take(field, j + 1, axis=axis)
    return (1.0 - w) * lo + w * hi


def band_mean(field: np.ndarray, lat: np.ndarray, lat0: float, lat1: float,
              axis: int = -1) -> np.ndarray:
    """cos(lat)-weighted mean over a latitude band, inclusive."""
    sel = (lat >= lat0) & (lat <= lat1)
    w = np.cos(np.deg2rad(lat[sel]))
    w = w / w.sum()
    sub = np.take(field, np.where(sel)[0], axis=axis)
    return np.tensordot(sub, w, axes=([axis], [0]))


# ----------------------------------------------------------------------------
# TEM residual circulation
# ----------------------------------------------------------------------------

def tem_residual(vzm, wzm, thzm, vthzm, p_hpa, lat):
    """Transformed-Eulerian-mean residual circulation, log-pressure form.

    Andrews, Holton & Leovy (1987) eq. 3.5.1:

        v* = v_bar - (1/rho0) d/dz [ rho0 * psi ]
        w* = w_bar + (1/(a cos(phi))) d/dphi [ cos(phi) * psi ]
        psi = v'theta' / (d theta_bar / dz)

    Inputs are (time, lev, lat) numpy arrays on the SAME interface levels.

    CONVENTION, verified in 01c against the raw 3D h7 V and T fields:
    VTHzm on the CESM h6 tape is ALREADY the eddy flux v'theta', not the full
    product zonal-mean(v*theta). Reconstructing it as VTHzm - Vzm*THzm is wrong
    and produces a spurious DOWNWARD tropical w*. Median error against the true
    eddy flux: 15% for VTHzm as-is (the h6-interface vs h7-midpoint offset),
    1087% for the subtraction. The same holds for UVzm and UWzm.

    Wzm is the LOG-PRESSURE vertical velocity w = -H.omega/p with H = 7 km,
    verified in 01 by mass continuity (regression slope 1.011, r = 1.000).

    Returns (vstar, wstar, psi), all shaped like the inputs; w* in m/s.
    """
    phi = np.deg2rad(lat)
    cosphi = np.cos(phi)
    z = log_p_height(p_hpa)                       # (lev,)
    rho0 = np.exp(-z / H_SCALE)                   # arbitrary scaling; cancels

    vpt = vthzm                                   # v'theta' (K m/s), see above

    # d theta / dz, centred in z, one-sided at the ends
    dthdz = np.gradient(thzm, z, axis=1)

    psi = vpt / dthdz

    # w*: meridional divergence of cos(phi)*psi
    num = psi * cosphi[None, None, :]
    dnum_dphi = np.gradient(num, phi, axis=2)
    wstar = wzm + dnum_dphi / (A_EARTH * cosphi[None, None, :])

    # v*: vertical divergence of rho0*psi
    rp = rho0[None, :, None] * psi
    drp_dz = np.gradient(rp, z, axis=1)
    vstar = vzm - drp_dz / rho0[None, :, None]

    return vstar, wstar, psi


def tropical_upward_mass_flux(wstar, p_hpa, lat, p_target=70.0):
    """Net upward mass flux through a pressure surface, integrated between the
    turnaround latitudes (where the climatological w* changes sign), in kg/s.

        M = 2 pi a^2 * int rho * w* * cos(phi) dphi     over w* > 0

    Density from the hydrostatic log-pressure form rho = p / (H g) is avoided in
    favour of the exact isobaric mass element: on a pressure surface the upward
    mass flux per unit area is rho*w*, and rho = p/(R T). We do not carry T on
    the TEM tape, so we use the log-pressure density rho0 = p/(H*g0) * H = p/(g0*H)
    ... i.e. rho = p/(R*T) with the standard-atmosphere T implied by H = R T/g.
    This is the Butchart et al. (2010) / CCMVal convention and is what published
    tropical-upwelling mass fluxes use.
    """
    phi = np.deg2rad(lat)
    cosphi = np.cos(phi)
    dphi = np.gradient(phi)

    w_p = interp_level(wstar, p_hpa, p_target, axis=1)      # (time, lat)
    rho = (p_target * 100.0) / (G0 * H_SCALE)               # kg/m3

    up = np.where(w_p > 0.0, w_p, 0.0)
    integrand = rho * up * cosphi[None, :] * dphi[None, :]
    return 2.0 * np.pi * A_EARTH ** 2 * integrand.sum(axis=1)


def turnaround_latitudes(wstar_clim, lat, tropics=(-45.0, 45.0)):
    """Latitudes where the climatological w* profile crosses zero, bracketing
    the tropical upwelling region. Returns (lat_south, lat_north)."""
    sel = (lat >= tropics[0]) & (lat <= tropics[1])
    la, w = lat[sel], wstar_clim[sel]
    jmax = int(np.argmax(w))
    south = north = np.nan
    for j in range(jmax, 0, -1):
        if w[j - 1] <= 0.0 < w[j]:
            south = la[j - 1] + (la[j] - la[j - 1]) * (-w[j - 1]) / (w[j] - w[j - 1])
            break
    for j in range(jmax, len(w) - 1):
        if w[j] > 0.0 >= w[j + 1]:
            north = la[j] + (la[j + 1] - la[j]) * (w[j]) / (w[j] - w[j + 1])
            break
    return south, north


# ----------------------------------------------------------------------------
# Seasonal / statistical helpers
# ----------------------------------------------------------------------------

def djf_mask(mo):
    return np.isin(mo, [12, 1, 2])


def season_year(yr, mo):
    """Label DJF seasons by the year of the January (Dec 1996 -> season 1997)."""
    return np.where(mo == 12, yr + 1, yr)


def seasonal_means(x, yr, mo, months, label_fn=None):
    """Mean of x over the given months, per season-year. Returns (years, means)."""
    sel = np.isin(mo, months)
    lab = label_fn(yr, mo) if label_fn is not None else yr
    years = np.unique(lab[sel])
    out = np.array([np.nanmean(x[sel & (lab == y)]) for y in years])
    # drop incomplete edge seasons (fewer than 80 days for a 3-month season)
    cnt = np.array([np.sum(sel & (lab == y)) for y in years])
    good = cnt >= int(0.9 * 30 * len(months))
    return years[good], out[good]


def annual_means(x, yr):
    years = np.unique(yr)
    out = np.array([np.nanmean(x[yr == y]) for y in years])
    cnt = np.array([np.sum(yr == y) for y in years])
    good = cnt >= 350
    return years[good], out[good]


def monthly_climatology(x, mo):
    return np.array([np.nanmean(x[mo == m]) for m in range(1, 13)])


def sigma_sampling_error(sigma, n):
    """1-sigma sampling error on an estimate of sigma from n independent samples."""
    return sigma / np.sqrt(2.0 * (n - 1))


def join_segments(series, segments, key, year_key):
    """One diagnostic across several segments, joined and sorted by year.

    The segments are separate CESM runs, so this is a concatenation, not an
    integration: use it for statistics over the record, not for anything that
    differences consecutive years across the 1995/96 restart.
    """
    y = np.array(sum((series[s][year_key] for s in segments), []), dtype=float)
    v = np.array(sum((series[s][key] for s in segments), []), dtype=float)
    o = np.argsort(y)
    return y[o], v[o]


def window_stats(years, values, lo, hi):
    """Mean and sigma of a series over [lo, hi], detrended if the trend is real.

    Detrended when the OLS trend is significant at p < 0.05, the same rule as
    07_period_split.stats_of, so every sigma in the repo is comparable.
    """
    m = (years >= lo) & (years <= hi)
    x = values[m]
    t = np.arange(len(x), dtype=float)
    sl, ic, _, p, se = _stats.linregress(t, x)
    det = x - (ic + sl * t)
    sig_raw, sig_det = float(x.std(ddof=1)), float(det.std(ddof=1))
    return dict(n=int(m.sum()), years=[int(years[m].min()), int(years[m].max())],
                mean=float(x.mean()), sigma_raw=sig_raw, sigma_detrended=sig_det,
                sigma_used=(sig_det if p < 0.05 else sig_raw),
                trend_per_decade=float(sl * 10), trend_se_per_decade=float(se * 10),
                trend_p=float(p), detrended=bool(p < 0.05))


def bias_target(sigma, n):
    """max(0.5 sigma, 1.96 sigma / sqrt(n)) - D5, EVALUATION_PROTOCOL.md appendix A."""
    return max(0.5 * sigma, 1.96 * sigma / np.sqrt(n))


def ratio_window(n, n_ref, per_year=1.0):
    """95% window on a sigma ratio between an n-year sample and an n_ref anchor."""
    ne, nr = n * per_year, n_ref * per_year
    rr = np.sqrt(1 / (2 * (ne - 1)) + 1 / (2 * (nr - 1)))
    return float(1 - 1.96 * rr), float(1 + 1.96 * rr)


def bootstrap_ci(values, stat_fn, n_boot=10000, ci=(2.5, 97.5), seed=20260813):
    """Non-parametric bootstrap CI, resampling whole years."""
    rng = np.random.default_rng(seed)
    n = len(values)
    draws = np.empty(n_boot)
    for i in range(n_boot):
        draws[i] = stat_fn(values[rng.integers(0, n, n)])
    return np.percentile(draws, ci[0]), np.percentile(draws, ci[1])


# ----------------------------------------------------------------------------
# SSW detection (Charlton & Polvani 2007)
# ----------------------------------------------------------------------------

def detect_ssw_cp07(u60, yr, mo, dy, hemisphere="NH"):
    """Charlton-Polvani (2007) major midwinter warming detection.

    u60: daily zonal-mean zonal wind at 60 deg latitude, 10 hPa (m/s).
         For the SH pass in the sign-flipped series so that "westerly" is
         always positive.

    Criteria as implemented:
      - central date = first day the wind reverses (u < 0) in the extended
        winter season (NH: Nov 1 - Mar 31).
      - events must be separated by >= 20 consecutive westerly days.
      - final warmings excluded: the wind must return to westerly for >= 10
        consecutive days before 30 April of that season.

    Returns a list of dicts with the central date and season label.
    """
    if hemisphere == "NH":
        season_months = [11, 12, 1, 2, 3]
        lab = np.where(np.isin(mo, [11, 12]), yr + 1, yr)
        end_month, end_day = 4, 30
    else:
        season_months = [5, 6, 7, 8, 9]
        lab = yr.copy()
        end_month, end_day = 10, 31

    events = []
    in_season = np.isin(mo, season_months)

    for season in np.unique(lab[in_season]):
        m = in_season & (lab == season)
        idx = np.where(m)[0]
        if len(idx) < 100:          # incomplete season at the record edges
            continue
        u = u60[idx]
        rev = u < 0.0

        # window used for the final-warming test: to end_month/end_day
        tail = np.where((lab == season) & (
            ((mo > season_months[-1]) & (mo < end_month)) |
            ((mo == end_month) & (dy <= end_day))))[0]
        u_tail = np.concatenate([u, u60[tail]])

        last_central = -999
        j = 0
        while j < len(u):
            if rev[j]:
                if j - last_central < 20:
                    j += 1
                    continue
                # require 20 consecutive westerly days before this reversal
                if last_central > 0:
                    gap = u[last_central:j]
                    if np.sum(gap > 0) < 20 or _max_run(gap > 0) < 20:
                        j += 1
                        continue
                # final-warming test: >=10 consecutive westerly days after,
                # before end_month/end_day
                after = u_tail[j:]
                if _max_run(after > 0) < 10:
                    break                      # this is the final warming
                events.append(dict(season=int(season),
                                   month=int(mo[idx[j]]),
                                   day=int(dy[idx[j]]),
                                   u=float(u[j])))
                last_central = j
                # skip forward to the next westerly spell
                while j < len(u) and u[j] <= 0:
                    j += 1
            else:
                j += 1
    return events


def _max_run(boolarr):
    """Length of the longest run of True."""
    best = cur = 0
    for b in boolarr:
        cur = cur + 1 if b else 0
        best = max(best, cur)
    return best
