"""
19 - Probe the CDS before committing to a 36-year ERA5 pull.

Three unknowns decide the shape of the bulk download, and guessing any of them
wrong costs days of transfer:

  1. Does derived-era5-pressure-levels-daily-statistics honour a `grid` key?
     At 1.0 degree the full 1990-2025 pull is ~164 GB; at native 0.25 degree it
     is ~2.6 TB. A factor of 16 on the schedule.
  2. What throughput and queue latency does CDS actually give? Daily statistics
     are computed at retrieval rather than served from an archive, so queueing
     dominates and cannot be estimated from the byte count.
  3. What exactly comes back - level list, latitude ordering (ERA5 runs
     +90 to -90, CESM ascending), variable names and units.

One day of data settles all three. Nothing here is part of the pipeline; it
writes a log and is read by a human.

Run with the download environment. Stages are a (one day native),
b (one day regridded) and c (one month, batching efficiency):
  ../era5_env/bin/python 19_era5_probe.py          # all three
  ../era5_env/bin/python 19_era5_probe.py b,c      # just these

Requests take minutes each and the log is appended, so stages can be run
separately. Do not wrap this in `timeout` - a killed request still occupies
a CDS queue slot.

Output: logs/19_probe.txt
"""
import os
import sys
import time
import traceback

DATASET = "derived-era5-pressure-levels-daily-statistics"
LEVELS = ["3", "5", "7", "10", "20", "30", "50", "70", "100", "125", "150", "175"]
VARIABLES = ["u_component_of_wind", "v_component_of_wind",
             "vertical_velocity", "temperature"]

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LOGDIR = os.path.join(_ROOT, "logs")
SCRATCH = os.path.join(_ROOT, "output", "era5_probe")


def base_request():
    return {
        "product_type": "reanalysis",
        "variable": VARIABLES,
        "year": "2000",
        "month": "01",
        "day": ["01"],
        "pressure_level": LEVELS,
        "daily_statistic": "daily_mean",
        "time_zone": "utc+00:00",
        "frequency": "1_hourly",
        "data_format": "netcdf",
        "download_format": "unarchived",
    }


def fetch(client, request, target, w):
    """Retrieve one request, timing it. Returns bytes written, or None on error."""
    t0 = time.time()
    try:
        client.retrieve(DATASET, request, target)
    except Exception as exc:
        w(f"  FAILED after {time.time() - t0:6.1f} s")
        w(f"  {type(exc).__name__}: {exc}")
        return None
    dt = time.time() - t0
    n = os.path.getsize(target)
    w(f"  ok  {n / 1e6:8.1f} MB in {dt:6.1f} s  ({n / 1e6 / max(dt, 1e-9):.2f} MB/s)")
    return n


def describe(path, w):
    """Report what actually came back: grid, levels, latitude order, units."""
    try:
        import netCDF4
    except ImportError:
        w("  netCDF4 not available; skipping inspection")
        return
    import zipfile
    if zipfile.is_zipfile(path):
        names = zipfile.ZipFile(path).namelist()
        w(f"  ZIP archive, {len(names)} members: {', '.join(names[:6])}")
        w("  -> ingest must unzip; one netCDF per variable")
        return
    try:
        ds = netCDF4.Dataset(path)
    except Exception as exc:
        w(f"  could not open as netCDF: {exc}")
        return
    w(f"  dims: {dict((k, len(v)) for k, v in ds.dimensions.items())}")
    for name in ds.variables:
        v = ds.variables[name]
        units = getattr(v, "units", "?")
        w(f"    {name:28s} {str(v.dimensions):38s} {units}")
    for cand in ("latitude", "lat"):
        if cand in ds.variables:
            lat = ds.variables[cand][:]
            order = "descending +90 to -90" if lat[0] > lat[-1] else "ascending -90 to +90"
            w(f"  latitude: n={len(lat)} {lat[0]:+.2f} to {lat[-1]:+.2f}, {order}")
            if len(lat) > 1:
                w(f"            spacing {abs(float(lat[1] - lat[0])):.3f} deg")
            break
    for cand in ("pressure_level", "level", "plev"):
        if cand in ds.variables:
            w(f"  levels: {list(ds.variables[cand][:])}")
            break
    ds.close()


def main():
    os.makedirs(LOGDIR, exist_ok=True)
    os.makedirs(SCRATCH, exist_ok=True)
    stages = sys.argv[1].split(",") if len(sys.argv) > 1 else ["a", "b", "c"]
    log = open(os.path.join(LOGDIR, "19_probe.txt"), "a")

    def w(line=""):
        print(line)
        log.write(line + "\n")
        log.flush()

    w("")
    w("=" * 72)
    w(f"ERA5 CDS probe - stages {','.join(stages)} - {time.strftime('%Y-%m-%d %H:%M')}")
    w(f"dataset   {DATASET}")
    w(f"levels    {len(LEVELS)}: {' '.join(LEVELS)} hPa")
    w(f"variables {', '.join(VARIABLES)}")
    w()

    try:
        import cdsapi
    except ImportError:
        w("cdsapi not installed. Run with ../era5_env/bin/python")
        return 1

    if not os.path.exists(os.path.expanduser("~/.cdsapirc")) \
            and not os.environ.get("CDSAPI_KEY"):
        w("NO CREDENTIALS. Write a CDS key to ~/.cdsapirc and re-run.")
        return 2

    try:
        client = cdsapi.Client()
    except Exception as exc:
        w(f"could not construct client: {exc}")
        return 2

    got = {}

    # -- A: one day, native grid. Per-day byte cost and the true layout.
    if "a" in stages:
        w("A. one day, native grid (no grid key)")
        path = os.path.join(SCRATCH, "probe_native.nc")
        got["a"] = fetch(client, base_request(), path, w)
        if got["a"]:
            describe(path, w)
        w()

    # -- B: does the server regrid? This is the factor of 16 on volume.
    if "b" in stages:
        w("B. one day, grid = [1.0, 1.0]")
        req = base_request()
        req["grid"] = [1.0, 1.0]
        path = os.path.join(SCRATCH, "probe_1deg.nc")
        got["b"] = fetch(client, req, path, w)
        if got["b"]:
            describe(path, w)
        w()

    # -- C: a whole month. A costs ~6.5 min of which ~7 s is transfer, so the
    #       schedule is set by per-request processing, not by bytes. If a month
    #       costs about what a day costs, 432 month-requests is the plan; if it
    #       scales with days, the pull is 30x longer and needs rethinking.
    if "c" in stages:
        w("C. one month (Jan 2000), native grid - batching efficiency")
        req = base_request()
        req["day"] = [f"{d:02d}" for d in range(1, 32)]
        path = os.path.join(SCRATCH, "probe_month.nc")
        got["c"] = fetch(client, req, path, w)
        if got["c"]:
            describe(path, w)
        w()

    w("IMPLICATIONS for 1990-2025 (36 years, 13149 days, 432 months)")
    if got.get("a"):
        w(f"  native 0.25 deg   {got['a'] / 1e6:7.1f} MB/day"
          f" -> {got['a'] * 13149 / 1e12:5.2f} TB total")
    if got.get("b"):
        w(f"  regridded 1.0 deg {got['b'] / 1e6:7.1f} MB/day"
          f" -> {got['b'] * 13149 / 1e12:5.2f} TB total")
        if got.get("a"):
            w(f"  grid key HONOURED - {got['a'] / got['b']:.1f}x less transfer at 1.0 deg")
    elif "b" in stages:
        w("  grid key REJECTED - stay native 0.25 deg, or switch to")
        w("  reanalysis-era5-complete, which is MARS-backed and supports grid/levelist")
    if got.get("c") and got.get("a"):
        w(f"  month request returned {got['c'] / got['a']:.1f}x one day's bytes"
          f" ({got['c'] / 1e6:.1f} MB)")
    w()
    w("Raw probe files are under output/era5_probe/ and can be deleted.")
    log.close()
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception:
        traceback.print_exc()
        sys.exit(1)
