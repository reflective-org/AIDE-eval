"""
21 - Fetch ERA5 monthly means from CDS and reduce them to a zonal-mean tape.

RESTRICTED SCOPE. Monthly means can support only part of the protocol, because
the eddy heat flux v'theta' has to be formed before any time average. From
monthly fields it is not recoverable, so mass flux and w* - and the R1
mechanism relation - cannot be computed at all from this tape. The SSW count
and the daily DJF sigma need the day the wind reverses and are equally out of
reach. What remains, and what 22/23 score, is:

    vortex_NH, vortex_SH, polar_cap_T_NH, polar_cap_T_SH

Those four are exact from monthly data, not approximate. Every step that
produces them - log-p interpolation to 10 hPa, the d ln p weighted average
over the 10-50 hPa layer, the cos(lat) band mean - is linear in the field, so
a day-weighted mean of monthly means equals the mean of the daily series. 22
recovers that day weighting by expanding each month across its own days, which
also lets the existing seasonal_means completeness rules apply unchanged.

Only two variables are needed: u for the vortex, t for the polar cap.

The monthly dataset's cost limit is 120000 against the daily dataset's 400, so
a year fits in one request with room to spare. Requests are still issued one
per year and run concurrently, which is how this script measures the CDS
concurrency limit - the number that decides whether a future daily pull takes
days or hours.

Run with the download environment:
  ../era5_env/bin/python 21_era5_monthly_ingest.py [first_year last_year]

Output: output/era5_monthly_tape.nc, logs/21_era5_monthly_ingest.txt
"""
import os
import sys
import time
import zipfile
import threading
import collections
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import netCDF4

DATASET = "reanalysis-era5-pressure-levels-monthly-means"
# 10 hPa for the vortex; 10-50 hPa for the polar cap. 7 and 70 bracket the
# range so the log-p interpolation to 10 hPa never extrapolates.
LEVELS = ["7", "10", "20", "30", "50", "70"]
VARIABLES = ["u_component_of_wind", "temperature"]
WORKERS = 5                      # concurrent requests; also the concurrency probe

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTDIR = os.path.join(_ROOT, "output")
LOGDIR = os.path.join(_ROOT, "logs")
RAW = os.path.join(OUTDIR, "era5_monthly_raw")
TAPE = os.path.join(OUTDIR, "era5_monthly_tape.nc")

_print_lock = threading.Lock()
_log_fh = None


def w(line=""):
    with _print_lock:
        print(line, flush=True)
        if _log_fh:
            _log_fh.write(line + "\n")
            _log_fh.flush()


def request_for(year):
    return {
        "product_type": ["monthly_averaged_reanalysis"],
        "variable": VARIABLES,
        "pressure_level": LEVELS,
        "year": [str(year)],
        "month": [f"{m:02d}" for m in range(1, 13)],
        "time": ["00:00"],
        "data_format": "netcdf",
        "download_format": "unarchived",
    }


def fetch_year(year, marks):
    """Download one year. Records timings so concurrency can be measured."""
    import cdsapi
    path = os.path.join(RAW, f"era5_monthly_{year}.nc")
    if os.path.exists(path) and os.path.getsize(path) > 0:
        w(f"  {year}  already present, skipping")
        marks.append((year, None, None))
        return path
    t0 = time.time()
    w(f"  {year}  submitting")
    cdsapi.Client(quiet=True, progress=False).retrieve(DATASET, request_for(year), path)
    t1 = time.time()
    w(f"  {year}  done in {(t1 - t0) / 60:5.1f} min  "
      f"{os.path.getsize(path) / 1e6:7.1f} MB")
    marks.append((year, t0, t1))
    return path


def open_members(path):
    """Yield netCDF datasets from a file that may be a bare .nc or a zip."""
    if zipfile.is_zipfile(path):
        z = zipfile.ZipFile(path)
        for name in z.namelist():
            yield netCDF4.Dataset(name, memory=z.read(name))
    else:
        yield netCDF4.Dataset(path)


def reduce_year(path):
    """Zonal-mean one year's file to (12, nlev, nlat), axes ascending."""
    found, lat, lev, ntime = {}, None, None, None
    for ds in open_members(path):
        for short, key in (("u", "Uzm"), ("t", "Tzm")):
            if short not in ds.variables:
                continue
            a = ds.variables[short]
            x = np.asarray(a[:], dtype="f8")          # (time, level, lat, lon)
            if x.ndim != 4:
                raise SystemExit(f"unexpected shape {x.shape} for {short} in {path}")
            found[key] = x.mean(axis=3)               # zonal mean over longitude
            lat = np.asarray(ds.variables["latitude"][:], dtype="f8")
            lev = np.asarray(ds.variables["pressure_level"][:], dtype="f8")
            ntime = x.shape[0]
        ds.close()
    missing = {"Uzm", "Tzm"} - set(found)
    if missing:
        raise SystemExit(f"{path}: missing {missing}")
    if ntime != 12:
        raise SystemExit(f"{path}: expected 12 months, got {ntime}")

    # ERA5 ships latitude +90 -> -90 and pressure high -> low. Both helpers in
    # aide_val_common assume ascending: interp_level uses searchsorted on
    # log(p), and interp_lat on latitude. Flipping here is the whole reason the
    # rest of the pipeline can be reused untouched.
    if lat[0] > lat[-1]:
        lat = lat[::-1]
        found = {k: v[:, :, ::-1] for k, v in found.items()}
    if lev[0] > lev[-1]:
        lev = lev[::-1]
        found = {k: v[:, ::-1, :] for k, v in found.items()}
    return found, lev, lat


def main():
    os.makedirs(RAW, exist_ok=True)
    os.makedirs(LOGDIR, exist_ok=True)
    global _log_fh
    _log_fh = open(os.path.join(LOGDIR, "21_era5_monthly_ingest.txt"), "a")

    lo, hi = (int(sys.argv[1]), int(sys.argv[2])) if len(sys.argv) > 2 else (1990, 1994)
    years = list(range(lo, hi + 1))

    w("")
    w("=" * 72)
    w(f"ERA5 monthly ingest {lo}-{hi}  -  {time.strftime('%Y-%m-%d %H:%M')}")
    w(f"dataset   {DATASET}")
    w(f"variables {', '.join(VARIABLES)}   levels {', '.join(LEVELS)} hPa")
    w(f"workers   {WORKERS} concurrent requests")
    w("scope     vortex_NH, vortex_SH, polar_cap_T_NH, polar_cap_T_SH only")
    w("          mass flux, w*, R1, SSW and daily sigma need sub-monthly data")
    w("")

    marks = []
    t_all = time.time()
    with ThreadPoolExecutor(max_workers=WORKERS) as ex:
        paths = list(ex.map(lambda y: fetch_year(y, marks), years))
    wall = time.time() - t_all

    # -- concurrency actually achieved, from the request intervals
    timed = [(a, b) for _, a, b in marks if a is not None]
    if timed:
        serial = sum(b - a for a, b in timed)
        edges = sorted([(a, 1) for a, _ in timed] + [(b, -1) for _, b in timed])
        cur = peak = 0
        for _, d in edges:
            cur += d
            peak = max(peak, cur)
        w("")
        w(f"CONCURRENCY  {len(timed)} requests, {WORKERS} workers")
        w(f"  peak simultaneous in flight   {peak}")
        w(f"  sum of request durations      {serial / 60:6.1f} min")
        w(f"  wall clock                    {wall / 60:6.1f} min")
        w(f"  effective speed-up            {serial / max(wall, 1e-9):5.2f}x")
        w(f"  mean per request              {serial / len(timed) / 60:6.1f} min")

    # -- reduce and concatenate
    w("")
    w("REDUCING")
    U, T, yy, mm = [], [], [], []
    lev = lat = None
    for year, path in zip(years, paths):
        f, lv, la = reduce_year(path)
        if lev is None:
            lev, lat = lv, la
        elif not (np.array_equal(lev, lv) and np.array_equal(lat, la)):
            raise SystemExit(f"{year}: grid differs from the first year")
        U.append(f["Uzm"])
        T.append(f["Tzm"])
        yy += [year] * 12
        mm += list(range(1, 13))
        w(f"  {year}  {f['Uzm'].shape} zonal-meaned")

    U = np.concatenate(U, axis=0)
    T = np.concatenate(T, axis=0)

    with netCDF4.Dataset(TAPE, "w") as ds:
        ds.createDimension("time", U.shape[0])
        ds.createDimension("level", len(lev))
        ds.createDimension("lat", len(lat))
        for name, dims, val, units in (
                ("Uzm", ("time", "level", "lat"), U, "m s-1"),
                ("Tzm", ("time", "level", "lat"), T, "K"),
                ("level", ("level",), lev, "hPa"),
                ("lat", ("lat",), lat, "degrees_north"),
                ("year", ("time",), np.array(yy), "1"),
                ("month", ("time",), np.array(mm), "1")):
            v = ds.createVariable(name, "f8" if val.dtype.kind == "f" else "i4", dims)
            v[:] = val
            v.units = units
        ds.source = DATASET
        ds.note = ("Zonal-mean ERA5 monthly means. Latitude and pressure both "
                   "ascending. Supports vortex and polar-cap diagnostics only; "
                   "v'theta' cannot be formed from monthly fields.")
        ds.years = f"{lo}-{hi}"

    w("")
    w(f"wrote {TAPE}")
    w(f"  Uzm/Tzm {U.shape}  levels {lev.tolist()}  lat {lat[0]:+.2f}..{lat[-1]:+.2f}")
    _log_fh.close()


if __name__ == "__main__":
    main()
