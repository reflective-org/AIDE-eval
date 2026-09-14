"""
20 - Build the standard series JSON from any zonal-mean source.

This is the general ingest. 07 derives the series from the CESM h6 tape; this
derives the same series from anything else - a reanalysis, another GCM, an
emulator rollout - so that 17 has exactly one scoring path to maintain rather
than one per dataset.

WHAT A SOURCE HAS TO PROVIDE. Zonal-mean u and T on pressure levels, with a
time axis. That is the whole contract. Everything else is optional and the
script reports what it could not compute rather than omitting it silently:

    supplied u, T           -> vortex NH/SH, polar cap T NH/SH, and their
                               12-month climatologies (the seasonal shape check)
    supplied v, omega, v'T' -> additionally the TEM residual circulation: mass
       at daily resolution     flux, w*, the w* profile, the 100 hPa heat flux,
                               and with daily u the SSW count and the DJF
                               percentiles. This is the full CESM set, reduced by
                               the same aide_val_common.series_from_tem_fields
                               that 07 uses on the h6 tape. v'T' must be formed
                               from the 3-D fields before time averaging; a
                               monthly source cannot reach any of these, since
                               SSW and the percentiles need the day the wind
                               reverses.

WHY MONTHLY MEANS ARE EXACT FOR THE FOUR THAT SURVIVE, not approximate. Each is
a chain of linear operations on the field - log-p interpolation to 10 hPa, the
d ln p weighted average across 10-50 hPa, the cos(lat) band mean, the seasonal
average. Linear operators commute with time averaging, so a day-weighted mean of
monthly means equals the mean of the daily series. The day weighting is
recovered by expanding each month across its own calendar days, which also lets
the completeness rules in seasonal_means apply untouched.

INDEPENDENCE IS A PROPERTY OF THE SOURCE, NOT OF ITS PERIOD. A reanalysis whose
years lie inside the anchor's is still an independent product; a window of the
anchor's own model is not, whatever its years. Each source declares this and 17
carries it into the report, so a self-consistency check can never be read as a
validation.

  20_series_from_zonal_mean.py NAME [first_year last_year]
  20_series_from_zonal_mean.py --list

Output: output/20_series__<NAME>.json
"""

import argparse
import calendar
import json
import os
import subprocess
import sys
import zipfile

import glob
import re

import numpy as np
import netCDF4
import xarray as xr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import aide_val_common as C

VORTEX_HPA = 10.0
CAP_HPA = (10.0, 50.0)
CAP_LAT = 60.0
SENTINEL = 1e20            # CLaMS writes +/-1e30 for missing; CESM 1e35 (D3)

# --------------------------------------------------------------------- sources
# A source is: where the fields are, what they are called, and whether it is an
# independent product. Adding a GCM or an emulator rollout is a few lines here
# plus a reader if its layout is new - no change to the scoring path.
CLAMS = dict(
    reader="clams_zip", var_u="U", var_t="TEMP",
    coord_lev="press", coord_lat="lat", time_resolution="monthly",
    calendar="gregorian", independent=True,
    note="ERA5-family fields regridded onto the CLaMS 1-degree zonal-mean grid, "
         "so lightly smoothed relative to the raw reanalysis")

# jcm-strat phase 10: the stripped JCM stratosphere (Dinosaur semi-Lagrangian dycore,
# Polvani-Kushner relaxation in place of radiation, troposphere and the tropical
# zonal-mean wind nudged to ERA5), 30 chained one-year segments 1990-2019 at
# T63 strat63, written 6-hourly. Read from the zonal-mean cache that
# /data/jcm30_validation/build_zm_cache.py reduced from the 2.7 TB of chunk files:
# [u] [v] [T] [omega] and the products [vT] etc., so v'T' = [vT] - [v][T] is
# formed at the 6 h output cadence, before any time averaging (protocol 4.6).
# Pressure is the hybrid coordinate times 1000 hPa: exact where hybrid_b = 0
# (above 80 hPa) and within 0.7 hPa of the true pressure at 100 hPa.
JCM = dict(
    reader="jcm_zmcache", cache="/data/jcm30_validation/zmcache",
    p0_hpa=1000.0, time_resolution="daily", calendar="gregorian",
    independent=True,
    note="stratosphere-only JCM on the Dinosaur dycore; troposphere (p > 150 hPa) "
         "and the tropical zonal-mean wind (QBO) nudged to ERA5, stratospheric "
         "temperature relaxed to a seasonal Polvani-Kushner equilibrium; daily means "
         "of 6-hourly instantaneous output, v'T' formed at 6 h")

SOURCES = {
    "JCM-strat-p10": dict(JCM, long_name="jcm-strat phase 10 (T63 strat63) 1990-2019"),
    "ERA5": dict(CLAMS, archive="/data/CLaMS/CLaMS_v3/clams_v3.1_era5_zm_lat.zip",
                 member="clams_v3.1_era5_zm_lat_press_{year}.nc",
                 long_name="ERA5 (CLaMS v3.1 zonal mean)"),
    "ERA-Interim": dict(CLAMS, archive="/data/CLaMS/CLaMS_v1/clams_v1.0_erain_zm_lat.zip",
                        member="clams_v1.0_erain_zm_lat_press_{year}.nc",
                        long_name="ERA-Interim (CLaMS v1.0 zonal mean)"),
    "JRA-55": dict(CLAMS, archive="/data/CLaMS/CLaMS_v1/clams_v1.0_jra55_zm_lat.zip",
                   member="clams_v1.0_jra55_zm_lat_press_{year}.nc",
                   long_name="JRA-55 (CLaMS v1.0 zonal mean)"),
    "MERRA-2": dict(CLAMS, archive="/data/CLaMS/CLaMS_v1/clams_v1.0_merra2_zm_lat.zip",
                    member="clams_v1.0_merra2_zm_lat_press_{year}.nc",
                    long_name="MERRA-2 (CLaMS v1.0 zonal mean)"),
}


# --------------------------------------------------------------------- readers
def read_clams_zip(spec, years):
    """One netCDF per year inside a zip, dims (month, level, lat).

    Streamed through memory because the pre-extracted directories beside the
    archives are not readable by every uid on this machine.
    """
    U, T, yr, mo, lev, lat = [], [], [], [], None, None
    have = set(zipfile.ZipFile(spec["archive"]).namelist())
    for y in years:
        member = spec["member"].format(year=y)
        if member not in have:
            continue
        raw = subprocess.run(["unzip", "-p", spec["archive"], member],
                             capture_output=True, check=True).stdout
        d = netCDF4.Dataset(f"{y}", memory=raw)
        u = np.asarray(d.variables[spec["var_u"]][:], dtype="f8")
        t = np.asarray(d.variables[spec["var_t"]][:], dtype="f8")
        if lev is None:
            lev = np.asarray(d.variables[spec["coord_lev"]][:], dtype="f8")
            lat = np.asarray(d.variables[spec["coord_lat"]][:], dtype="f8")
        n = u.shape[0]
        U.append(u); T.append(t)
        yr.append(np.full(n, y, int)); mo.append(np.arange(1, n + 1))
        d.close()
    if not U:
        raise SystemExit("no members matched the requested years")
    return dict(U=np.concatenate(U), T=np.concatenate(T), lev=lev, lat=lat,
                yr=np.concatenate(yr), mo=np.concatenate(mo))


def read_netcdf(spec, years):
    """A single netCDF holding the whole record, dims (time, level, lat).

    The path a GCM or an emulator rollout would take. It needs `year` and
    `month` variables, or a CF time axis the caller has already decoded.
    """
    with netCDF4.Dataset(spec["path"]) as d:
        U = np.asarray(d.variables[spec["var_u"]][:], dtype="f8")
        T = np.asarray(d.variables[spec["var_t"]][:], dtype="f8")
        lev = np.asarray(d.variables[spec["coord_lev"]][:], dtype="f8")
        lat = np.asarray(d.variables[spec["coord_lat"]][:], dtype="f8")
        yr = np.asarray(d.variables["year"][:], dtype=int)
        mo = np.asarray(d.variables["month"][:], dtype=int)
    m = np.isin(yr, years)
    return dict(U=U[m], T=T[m], lev=lev, lat=lat, yr=yr[m], mo=mo[m])


def read_jcm_zmcache(spec, years):
    """The jcm30_validation zonal-mean cache: one file per 10-day chunk holding
    6-hourly [u] [v] [T] [omega] [vT] on (time, level, lat), one directory per
    calendar year.

    Every chunk of every requested year is concatenated before the daily mean is
    taken, so the 00:00 frame that ends one year lands in the first day of the
    next together with that day's own frames - the segments are one chained
    integration, not thirty separate runs. Days are labelled by the calendar day
    of their frames; a chunk's last frame is stamped 00:00 of the following day.
    """
    files = []
    for y in years:
        fs = glob.glob(os.path.join(spec["cache"], str(y), "longrun_day*.nc"))
        files += sorted(fs, key=lambda f: int(re.search(r"_day(\d+)\.nc$", f).group(1)))
    if not files:
        raise SystemExit(f"no cache files under {spec['cache']} for the requested years")
    parts = []
    for f in files:
        with xr.open_dataset(f) as d:
            d = d[["u", "v", "T", "omega", "vT"]].load()
        d["vpTp"] = d["vT"] - d["v"] * d["T"]          # eddy flux at the output cadence
        parts.append(d.drop_vars("vT"))
    six = xr.concat(parts, dim="time").sortby("time")
    day = six.resample(time="1D").mean()
    t = day.time.values
    yr = t.astype("datetime64[Y]").astype(int) + 1970
    # empty bins, if any, and the lone 00:00 frame that closes the last year
    have = ~np.isnan(day["u"].isel(level=0, lat=0).values) & np.isin(yr, list(years))
    day = day.isel(time=have); t, yr = t[have], yr[have]
    mo = t.astype("datetime64[M]").astype(int) % 12 + 1
    dy = (t - t.astype("datetime64[M]")).astype("timedelta64[D]").astype(int) + 1
    f64 = lambda v: day[v].values.astype("f8")
    return dict(U=f64("u"), V=f64("v"), T=f64("T"), OMEGA=f64("omega"),
                VPTP=f64("vpTp"), lev=day.level.values.astype("f8") * spec["p0_hpa"],
                lat=day.lat.values.astype("f8"), yr=yr, mo=mo, dy=dy)


READERS = {"clams_zip": read_clams_zip, "netcdf": read_netcdf,
           "jcm_zmcache": read_jcm_zmcache}


# ----------------------------------------------------------------------- core
FIELDS = ("U", "T", "V", "OMEGA", "VPTP")


def orient(F):
    """Put levels and latitudes ascending, as interp_level and band_mean assume."""
    if F["lev"][0] > F["lev"][-1]:
        F["lev"] = F["lev"][::-1]
        F.update({k: F[k][:, ::-1, :] for k in FIELDS if k in F})
    if F["lat"][0] > F["lat"][-1]:
        F["lat"] = F["lat"][::-1]
        F.update({k: F[k][:, :, ::-1] for k in FIELDS if k in F})
    return F


def expand_to_days(x, yr, mo, cal):
    """Repeat each month across its own calendar days, recovering day weighting.

    A `noleap` source gets 28-day Februaries; a Gregorian one gets 29 in leap
    years. Protocol section 4.3: the calendar changes only the day counts that
    the completeness rules see.
    """
    if cal in ("noleap", "365_day"):
        length = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        reps = np.array([length[m - 1] for m in mo])
    else:
        reps = np.array([calendar.monthrange(int(y), int(m))[1]
                         for y, m in zip(yr, mo)])
    return np.repeat(x, reps, axis=0), np.repeat(yr, reps), np.repeat(mo, reps)


def series_from_tem(F, spec):
    """The full CESM series from a daily source carrying u, v, omega, T and v'T'.

    The three derivations of protocol section 4.6, then the same reduction 07
    applies to the h6 tape: w = -H omega / p, theta = T (1000/p)^kappa, and
    v'theta' = v'T' (1000/p)^kappa - exact on a pressure surface because the
    factor is constant along it.
    """
    lev, lat = F["lev"], F["lat"]
    fac = (C.P_REF_HPA / lev)[None, :, None] ** C.KAPPA
    W = -C.H_SCALE * F["OMEGA"] / (lev[None, :, None] * 100.0)
    out = C.series_from_tem_fields(F["U"], F["V"], W, F["T"] * fac, F["VPTP"] * fac,
                                   lev, lat, F["yr"], F["mo"], F["dy"])
    supports = ["mass_flux", "w_star", "vortex_NH", "vortex_SH",
                "polar_cap_T_NH", "polar_cap_T_SH", "heat_flux_100",
                "w_star_profile", "ssw_NH", "daily_distribution", "daily_sigma"]
    return out, supports, {}


def series_from_fields(U, T, lev, lat, yr, mo, spec):
    """The standard series, plus an explicit account of what is missing.

    Every reduction here goes through aide_val_common, so a source and the CESM
    anchor travel the identical code path - protocol section 5, rule 1.
    """
    U = np.ma.filled(np.ma.masked_where(np.abs(U) > SENTINEL, U), np.nan)
    T = np.ma.filled(np.ma.masked_where(np.abs(T) > SENTINEL, T), np.nan)

    u10 = C.interp_level(U, lev, VORTEX_HPA, axis=1)
    u60n = C.interp_lat(u10, lat, CAP_LAT, axis=-1)
    u60s = -C.interp_lat(u10, lat, -CAP_LAT, axis=-1)     # westerly positive

    ksel = (lev >= CAP_HPA[0]) & (lev <= CAP_HPA[1])
    if ksel.sum() < 2:
        raise SystemExit(f"source has {ksel.sum()} level(s) in "
                         f"{CAP_HPA[0]}-{CAP_HPA[1]} hPa; the polar cap layer "
                         f"needs at least two")
    w = np.gradient(np.log(lev[ksel])); w /= w.sum()
    tlev = np.tensordot(T[:, ksel, :], w, axes=([1], [0]))
    tcapN = C.band_mean(tlev, lat, CAP_LAT, 90.0, axis=-1)
    tcapS = C.band_mean(tlev, lat, -90.0, -CAP_LAT, axis=-1)

    # A source carrying v, omega and v'T' at daily resolution takes
    # series_from_tem instead and scores the full set.

    daily = spec["time_resolution"] == "monthly"
    cal = spec.get("calendar", "gregorian")
    if daily:
        u60n_d, yr_d, mo_d = expand_to_days(u60n, yr, mo, cal)
        u60s_d, _, _ = expand_to_days(u60s, yr, mo, cal)
        tcapN_d, _, _ = expand_to_days(tcapN, yr, mo, cal)
        tcapS_d, _, _ = expand_to_days(tcapS, yr, mo, cal)
    else:
        u60n_d, u60s_d, tcapN_d, tcapS_d, yr_d, mo_d = (
            u60n, u60s, tcapN, tcapS, yr, mo)

    out = {}
    sN, out["vortex_NH"] = C.seasonal_means(u60n_d, yr_d, mo_d, [12, 1, 2],
                                            C.season_year)
    _, out["polar_cap_T_NH"] = C.seasonal_means(tcapN_d, yr_d, mo_d, [12, 1, 2],
                                                C.season_year)
    out["_djf_years"] = sN
    sS, out["vortex_SH"] = C.seasonal_means(u60s_d, yr_d, mo_d, [6, 7, 8])
    _, out["polar_cap_T_SH"] = C.seasonal_means(tcapS_d, yr_d, mo_d, [6, 7, 8])
    out["_jja_years"] = sS

    # the 12-month climatology per year, for the seasonal-cycle shape check
    monthly = {"vortex_NH": u60n, "vortex_SH": u60s,
               "polar_cap_T_NH": tcapN, "polar_cap_T_SH": tcapS}
    years = np.array([y for y in np.unique(yr)
                      if set(mo[yr == y]) == set(range(1, 13))])
    out["_monthly_years"] = years
    for k, v in monthly.items():
        out[f"_monthly_{k}"] = np.array(
            [[float(np.nanmean(v[(yr == y) & (mo == m)])) for m in range(1, 13)]
             for y in years])

    supports = ["vortex_NH", "vortex_SH", "polar_cap_T_NH", "polar_cap_T_SH"]
    not_evaluable = {
        "mass_flux": "needs the eddy heat flux v'theta', which cannot be formed "
                     "after time averaging, and omega above 50 hPa",
        "w_star": "same as mass flux",
        "w_star_profile": "same as mass flux",
        "heat_flux_100": "needs v'theta'; the R1 mechanism relation goes with it",
        "ssw_NH": "needs the day the wind reverses; a monthly mean cannot resolve it",
        "daily_distribution": "needs the daily series",
        "daily_sigma": "needs the daily series",
    }
    if spec["time_resolution"] != "monthly":
        for k in ("ssw_NH", "daily_distribution", "daily_sigma"):
            not_evaluable.pop(k, None)
    return out, supports, not_evaluable


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1].strip())
    ap.add_argument("name", nargs="?", help="source name; --list to see them")
    ap.add_argument("years", nargs="*", type=int, help="first_year last_year")
    ap.add_argument("--list", action="store_true")
    a = ap.parse_args()
    if a.list or not a.name:
        for k, v in SOURCES.items():
            print(f"  {k:14s} {v['long_name']}")
        return
    if a.name not in SOURCES:
        raise SystemExit(f"unknown source {a.name!r}; --list to see them")

    spec = SOURCES[a.name]
    years = (range(a.years[0], a.years[1] + 1) if len(a.years) > 1
             else range(1900, 2101))
    F = orient(READERS[spec["reader"]](spec, years))
    U, T, lev, lat, yr, mo = (F[k] for k in ("U", "T", "lev", "lat", "yr", "mo"))
    print(f"{a.name}: {U.shape[0]} times {yr.min()}-{yr.max()}, "
          f"{lev.size} levels {lev.min():g}-{lev.max():g} hPa, {lat.size} lats")

    full = all(k in F for k in ("V", "OMEGA", "VPTP", "dy")) \
        and spec["time_resolution"] == "daily"
    series, supports, not_evaluable = (series_from_tem(F, spec) if full else
                                       series_from_fields(U, T, lev, lat, yr, mo, spec))
    res = {
        "name": a.name, "long_name": spec["long_name"],
        "independent": bool(spec["independent"]),
        "time_resolution": spec["time_resolution"],
        "period": [int(yr.min()), int(yr.max())],
        "levels_hPa": [float(x) for x in lev],
        "note": spec.get("note", ""),
        "supports": supports, "not_evaluable": not_evaluable,
        "series": {k: (np.asarray(v).tolist() if not np.isscalar(v) else v)
                   for k, v in series.items()},
    }
    out = os.path.join(C.OUTDIR, f"20_series__{a.name}.json")
    with open(out, "w") as f:
        json.dump(res, f, indent=2)
    for k in supports:
        if k in series:
            v = np.asarray(series[k], float)
            print(f"  {k:16s} n={v.size:3d}  mean {v.mean():9.3f}")
    if "_ssw_count" in series:
        print(f"  {'ssw_NH':16s} {series['_ssw_count']} events in "
              f"{series['_ssw_winters']} winters")
    print(f"  not evaluable: {', '.join(sorted(not_evaluable)) or 'nothing'}")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
