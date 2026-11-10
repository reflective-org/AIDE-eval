"""
Shared sources, readers and math for the transport scorecard (22, 23, 24).

The transport scorecard is the jcm-strat relative scorecard (phase 29) ported into
this repository: age of air (A), residual circulation (B), polar vortex (C) and
variability (D), each metric scored as its relative deviation from a reference.
docs/TRANSPORT_SCORECARD.md is the description; this module holds the parts
that 22, 23 and 24 share, so that they cannot disagree.

SOURCES AND READERS follow 20_series_from_zonal_mean.py: a source is a dict
naming its reader and where its data are, a reader is `reader(spec, years)`,
and READERS maps the name to the function. The transport readers return the
field contract of docs/TRANSPORT_SCORECARD.md section 1 rather than 20's
(U, T, V, OMEGA, VPTP) set, so they are registered here and not in 20.
"""

from __future__ import annotations

import glob
import os

import numpy as np
import xarray as xr

import aide_val_common as C

DAYS_PER_YEAR = 365.25          # the jcm-strat clocks count days; the contract asks years

# The contract (docs/TRANSPORT_SCORECARD.md section 1): required, then optional.
CONTRACT = ("Uzm", "Vzm", "THzm", "VTHzm", "AOA150")
OPTIONAL = ("T", "Wzm")
UNITS = {"Uzm": "m/s", "Vzm": "m/s", "THzm": "K", "VTHzm": "K m/s",
         "AOA150": "yr", "T": "K", "Wzm": "m/s"}

# --------------------------------------------------------------------- sources
# The eight jcm-strat phase 29 scorecard runs (QBO30 and QBO10dL76 excluded):
# T63 L95 stratosphere-only JCM, troposphere nudged to ERA5, daily zonal means
# 2006-2010 reduced from 6-hourly output by jcm-strat scripts/daily_zm.py
# (v'theta' formed within each 6-hourly snapshot). Read in place, never copied.
JCM_CACHE = "/data/JCM_stripped/jcm-strat-phase29/cache/variability"
JCM = dict(reader="jcm_npz", time_resolution="daily", calendar="gregorian",
           independent=True)
_RUNS = {
    "P20b": "baseline; ERA5 nudging top 400 hPa everywhere, no T_eq correction",
    "AL7": "P20b with the AL7 equilibrium-temperature correction (jcm-strat phase 22)",
    "AL7h": "AL7 correction at half amplitude (phase 25)",
    "AL7t": "AL7 correction tapered to zero at 35 deg, full from 50 deg (phase 24)",
    "L0": "P20b with the nudging top at 150 hPa in the tropics, 400 hPa in the "
          "extratropics (phase 23)",
    "LAL7": "L0 nudging top with the AL7 correction (phase 23)",
    "QBO10d": "AL7h with the QBO nudging timescale 10 d instead of 1 d (phase 27)",
    "LQBO10d": "QBO10d with L0's nudging top (phase 29)",
}
SOURCES = {n: dict(JCM, path=os.path.join(JCM_CACHE, f"{n}_0610.npz"),
                   long_name=f"jcm-strat {n} (T63 L95)", note=note)
           for n, note in _RUNS.items()}

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


def read_jcm_npz(spec, years):
    """One jcm-strat daily_zm.py npz: days, p_hpa, lat and (day, lev, lat) float32
    u, v, T, th, vth, aoa150 [, aoa_sfc]. The clock is in days; returned in years."""
    z = np.load(spec["path"])
    days = z["days"].astype("datetime64[D]")
    keep = np.isin(days.astype("datetime64[Y]").astype(int) + 1970, list(years))
    F = {new: z[old][keep] for old, new in JCM_NAMES.items()}
    F["AOA150"] = F["AOA150"] / np.float32(DAYS_PER_YEAR)
    F.update(plev=z["p_hpa"].astype("f8"), lat=z["lat"].astype("f8"), days=days[keep])
    return orient(F)


def read_netcdf(spec, years):
    """Any model that writes the contract names: one or more netCDF files
    (spec["path"] is a glob), dims (time, plev, lat), plev in hPa.

    Any calendar: each day is labelled by its own year, month and day, so a
    noleap record simply has no 29 February. Values beyond 1e20 are masked as
    fill (protocol section 5, rule 4). T, if absent, is THzm (p/1000)^kappa.
    """
    files = sorted(glob.glob(spec["path"]))
    if not files:
        raise SystemExit(f"no files match {spec['path']}")
    ds = xr.concat([xr.open_dataset(f) for f in files], dim="time")
    missing = [k for k in CONTRACT if k not in ds]
    if missing:
        raise SystemExit(f"{spec['path']}: missing contract fields {missing}")
    yr, mo, dy = (ds.time.dt.year.values, ds.time.dt.month.values, ds.time.dt.day.values)
    keep = np.isin(yr, list(years))
    days = np.array([f"{y:04d}-{m:02d}-{d:02d}" for y, m, d in
                     zip(yr[keep], mo[keep], dy[keep])], dtype="datetime64[D]")
    F = {}
    for k in CONTRACT + OPTIONAL:
        if k in ds:
            x = ds[k].transpose("time", "plev", "lat").values[keep]
            F[k] = np.where(np.abs(x) < C.FILL_SENTINEL, x, np.nan)
    F.update(plev=ds.plev.values.astype("f8"), lat=ds.lat.values.astype("f8"), days=days)
    if "T" not in F:
        F["T"] = F["THzm"] * (F["plev"][None, :, None] / C.P_REF_HPA) ** C.KAPPA
    ds.close()
    return orient(F)


READERS = {"jcm_npz": read_jcm_npz, "netcdf": read_netcdf}


def load_source(name, years):
    """The contract fields of one registered source, levels and latitudes ascending."""
    spec = SOURCES[name]
    return READERS[spec["reader"]](spec, years)
