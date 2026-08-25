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

Run with the download environment:
  ../era5_env/bin/python 19_era5_probe.py

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
    log = open(os.path.join(LOGDIR, "19_probe.txt"), "w")

    def w(line=""):
        print(line)
        log.write(line + "\n")
        log.flush()

    w("ERA5 CDS probe")
    w(f"dataset   {DATASET}")
    w(f"levels    {len(LEVELS)}: {' '.join(LEVELS)} hPa")
    w(f"variables {', '.join(VARIABLES)}")
    w("probe day 2000-01-01")
    w()

    try:
        import cdsapi
    except ImportError:
        w("cdsapi not installed. Run with ../era5_env/bin/python")
        return 1

    rc = os.path.expanduser("~/.cdsapirc")
    if not os.path.exists(rc) and not os.environ.get("CDSAPI_KEY"):
        w("NO CREDENTIALS.")
        w("  Register at https://cds.climate.copernicus.eu, accept the ERA5 licence,")
        w("  and write the key to ~/.cdsapirc. Then re-run this probe.")
        return 2

    try:
        client = cdsapi.Client()
    except Exception as exc:
        w(f"could not construct client: {exc}")
        return 2

    # -- A: native grid, to establish the per-day byte cost and the true layout
    w("A. native grid (no grid key)")
    a_path = os.path.join(SCRATCH, "probe_native.nc")
    a_bytes = fetch(client, base_request(), a_path, w)
    if a_bytes:
        describe(a_path, w)
    w()

    # -- B: does the server regrid? This is the factor of 16.
    w("B. grid = [1.0, 1.0]")
    b_req = base_request()
    b_req["grid"] = [1.0, 1.0]
    b_path = os.path.join(SCRATCH, "probe_1deg.nc")
    b_bytes = fetch(client, b_req, b_path, w)
    if b_bytes:
        describe(b_path, w)
    w()

    # -- what this implies for the real pull
    w("IMPLICATIONS for 1990-2025 (36 years, 13149 days)")
    for tag, n in (("native 0.25 deg", a_bytes), ("regridded 1.0 deg", b_bytes)):
        if n:
            w(f"  {tag:20s} {n / 1e6:7.1f} MB/day -> {n * 13149 / 1e12:6.2f} TB total")
    if b_bytes and a_bytes:
        w(f"  grid key is HONOURED - use 1.0 deg, {a_bytes / b_bytes:.1f}x less transfer")
    elif a_bytes and not b_bytes:
        w("  grid key REJECTED - fall back to native 0.25 deg, or to")
        w("  reanalysis-era5-complete, which is MARS-backed and supports grid/levelist")
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
