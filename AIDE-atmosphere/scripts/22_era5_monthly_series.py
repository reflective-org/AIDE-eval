"""
22 - Turn the ERA5 monthly tape into the four per-season series tier 1 can score.

Reads output/era5_monthly_tape.nc (written by 21) and produces the same
quantities 07_period_split.py derives from the CESM daily tape, for the four
diagnostics that survive a monthly archive:

    vortex_NH        DJF mean of u at 10 hPa, 60N
    vortex_SH        JJA mean of -u at 10 hPa, 60S   (westerly positive)
    polar_cap_T_NH   DJF mean of T over 10-50 hPa, 60-90N
    polar_cap_T_SH   JJA mean of T over 10-50 hPa, 90-60S

WHY MONTHLY IS EXACT HERE, NOT APPROXIMATE. Each of those is a chain of
linear operations on the field - log-p interpolation to 10 hPa, the d ln p
weighted average across the 10-50 hPa layer, the cos(lat) band mean, and the
seasonal average. Linear operators commute with time averaging, so a
day-weighted mean of monthly means is identical to the mean of the daily
series. The day weighting is recovered by expanding each month across its own
calendar days, which has the further benefit that seasonal_means applies its
normal completeness rules (>= 81 days per season) untouched.

The expansion produces a piecewise-constant series. That is harmless for these
four seasonal means and would be meaningless for anything reading day-to-day
variability, so this script deliberately does not emit _u60n_djf_daily,
_ssw_seasons or any annual series - the diagnostics that need them cannot be
computed from monthly data at all (see 21's header).

Run with the pinned analysis environment:
  ../../.AIDE-eval_env/bin/python 22_era5_monthly_series.py

Output: output/22_era5_monthly_series.json
"""
import os
import sys
import json
import calendar

import numpy as np
import netCDF4

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import aide_val_common as C

TAPE = os.path.join(C.OUTDIR, "era5_monthly_tape.nc")
OUT = os.path.join(C.OUTDIR, "22_era5_monthly_series.json")


def load_tape():
    with netCDF4.Dataset(TAPE) as ds:
        U = np.asarray(ds.variables["Uzm"][:], dtype="f8")
        T = np.asarray(ds.variables["Tzm"][:], dtype="f8")
        lev = np.asarray(ds.variables["level"][:], dtype="f8")
        lat = np.asarray(ds.variables["lat"][:], dtype="f8")
        yr = np.asarray(ds.variables["year"][:], dtype=int)
        mo = np.asarray(ds.variables["month"][:], dtype=int)
    if not np.all(np.diff(lev) > 0):
        raise SystemExit("tape levels are not ascending; interp_level would misread them")
    if not np.all(np.diff(lat) > 0):
        raise SystemExit("tape latitudes are not ascending")
    return U, T, lev, lat, yr, mo


def expand_to_days(x, yr, mo):
    """Repeat each month across its calendar days, recovering the day weighting.

    Returns (x_daily, yr_daily, mo_daily). Gregorian, so leap years get 29 days
    in February - ERA5's calendar, not CESM's noleap. Protocol section 4.3 says
    a different calendar changes only the day counts in the completeness rules.
    """
    reps = np.array([calendar.monthrange(int(y), int(m))[1] for y, m in zip(yr, mo)])
    return (np.repeat(x, reps, axis=0),
            np.repeat(yr, reps),
            np.repeat(mo, reps))


def main():
    U, T, lev, lat, yr, mo = load_tape()
    print(f"tape {U.shape}  levels {lev.tolist()}  "
          f"{yr.min()}-{yr.max()}  {len(yr)} months")

    # -- the same reductions 07_period_split.py applies to the CESM tape
    u10 = C.interp_level(U, lev, 10.0, axis=1)
    u60n = C.interp_lat(u10, lat, 60.0, axis=-1)
    u60s = -C.interp_lat(u10, lat, -60.0, axis=-1)          # westerly positive

    ksel = (lev >= 10.0) & (lev <= 50.0)
    wgt = np.gradient(np.log(lev[ksel]))
    wgt /= wgt.sum()
    tlev = np.tensordot(T[:, ksel, :], wgt, axes=([1], [0]))
    tcapN = C.band_mean(tlev, lat, 60.0, 90.0, axis=-1)
    tcapS = C.band_mean(tlev, lat, -90.0, -60.0, axis=-1)
    print(f"polar cap layer uses {ksel.sum()} levels inside 10-50 hPa "
          f"({lev[ksel].tolist()}); CESM has 8")

    # -- day-weight, then reduce exactly as the CESM path does
    out, series = {}, {}
    for name, x, months, label in (
            ("vortex_NH", u60n, [12, 1, 2], C.season_year),
            ("polar_cap_T_NH", tcapN, [12, 1, 2], C.season_year),
            ("vortex_SH", u60s, [6, 7, 8], None),
            ("polar_cap_T_SH", tcapS, [6, 7, 8], None)):
        xd, yd, md = expand_to_days(x, yr, mo)
        years, vals = C.seasonal_means(xd, yd, md, months, label)
        series[name] = (years, vals)
        out[name] = [float(v) for v in vals]
        out[("_djf_years" if 12 in months else "_jja_years")] = \
            [int(t) for t in years]
        print(f"  {name:16s} {len(vals)} seasons "
              f"{years.min()}-{years.max()}  "
              f"mean {vals.mean():8.3f}")

    doc = {
        "source": "ERA5 monthly means via CDS "
                  "(reanalysis-era5-pressure-levels-monthly-means)",
        "tape": os.path.basename(TAPE),
        "years": [int(yr.min()), int(yr.max())],
        "levels_hPa": lev.tolist(),
        "n_levels_in_polar_cap_layer": int(ksel.sum()),
        "scope": "monthly-recoverable diagnostics only",
        "not_computable_from_monthly": {
            "mass_flux": "needs v'theta', which cannot be formed after time averaging",
            "w_star": "needs v'theta'",
            "heat_flux_100": "is v'theta'",
            "ssw_count": "needs the day the wind reverses",
            "daily_DJF_sigma": "needs the daily series",
        },
        "exactness": "linear reductions commute with time averaging, so these "
                     "four seasonal means are exact given day weighting",
        "series": out,
    }
    with open(OUT, "w") as f:
        json.dump(doc, f, indent=2)
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
