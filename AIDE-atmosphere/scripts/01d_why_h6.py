"""
01d - Why the eddy fluxes are taken from h6 rather than rebuilt from the 3-D files.

The h6 tape carries VTHzm, the eddy heat flux v'theta', accumulated by CAM at
every 1800 s model timestep and written as a daily mean. The 3-D tapes carry V
and T, but both are TIME AVERAGES (h1 hourly means, h7 daily means), so any eddy
flux rebuilt from them loses the covariance carried by fluctuations inside the
averaging window.

This script measures how much is lost, by computing the same quantity three ways
for a set of sample days:

  A  model      VTHzm from h6            - accumulated every 1800 s
  B  hourly     eddy flux from the 24 hourly-mean h1 snapshots, then time-averaged
  C  daily      eddy flux from the daily-mean h7 fields

B vs A isolates the loss from hourly averaging; C vs B isolates the further loss
from daily averaging. If C is far from A, rebuilding from the daily 3-D files is
not a viable substitute - which is the justification for using h6.

Only levels with hybm == 0 are used (p >= ~182 hPa upward), where the pressure of
a model level is exactly hyam*P0 and no surface pressure is needed.

Output: output/01d_why_h6.json
"""

import json
import os
import numpy as np
import xarray as xr

import aide_val_common as C

SEG = "1996-2014"
CASE = C.SEGMENTS[SEG]["case"]
ROOT = f"{C.ROOT}/{CASE}/{C.SEGMENTS[SEG]['subdir']}"
OUT = os.path.join(C.OUTDIR, "01d_why_h6.json")

DAYS = [10, 25, 340, 355]          # two January days, two December days (1996)
H1_OFFSET = -23                    # see the note in main(); verified in 01d log
TEST_P = [10.0, 30.0, 70.0, 100.0]
BAND = (45.0, 75.0)                # the R1 predictor band


def main():
    h1v = xr.open_dataset(f"{ROOT}/../hour_1/{CASE}.cam.h1.V.1996010100-2014123100.nc")
    h1t = xr.open_dataset(f"{ROOT}/../hour_1/{CASE}.cam.h1.T.1996010100-2014123100.nc")
    h7v = xr.open_dataset(f"{ROOT}/{CASE}.cam.h7.V.19960101-20141231.nc")
    h7t = xr.open_dataset(f"{ROOT}/{CASE}.cam.h7.T.19960101-20141231.nc")

    lev = h1v["lev"].values.astype(np.float64)
    lat = h1v["lat"].values.astype(np.float64)
    hyam = h1v["hyam"].values.astype(np.float64)
    hybm = h1v["hybm"].values.astype(np.float64)
    P0 = float(h1v["P0"].values)

    pure = hybm == 0.0
    kidx = [int(np.argmin(np.abs(lev - p))) for p in TEST_P]
    assert all(pure[k] for k in kidx), "a test level is not in the pure-pressure region"
    p_lev = hyam[kidx] * P0                       # Pa, exact
    theta_fac = (100000.0 / p_lev) ** C.KAPPA     # (nk,)

    grid = C.load_grid(SEG)
    ilev = grid["ilev"]
    vthzm = C.load_h6(SEG, "VTHzm")

    jsel = (lat >= BAND[0]) & (lat <= BAND[1])
    wj = np.cos(np.deg2rad(lat[jsel])); wj /= wj.sum()

    def band_of(field_klat):
        return float(np.dot(field_klat[jsel], wj))

    def eddy(V, T):
        """Zonal eddy heat flux from (nk, nlat, nlon) V and T on the test levels."""
        TH = T * theta_fac[:, None, None]
        Vb = V.mean(axis=2, keepdims=True)
        Tb = TH.mean(axis=2, keepdims=True)
        return ((V - Vb) * (TH - Tb)).mean(axis=2)      # (nk, nlat)

    rows = []
    print(f"{'day':>4} {'p_hPa':>7} {'A model':>10} {'B hourly':>10} {'C daily':>10} "
          f"{'B/A':>7} {'C/A':>7}   (K m/s, 45-75N)")
    print("-" * 74)

    for d in DAYS:
        # --- B: 24 hourly snapshots, eddy flux each hour, then averaged
        # TIME ALIGNMENT: both tapes stamp at the END of the averaging interval and
        # carry a degenerate first record (time_bnds = [0,0]). Day index d on h6/h7
        # therefore corresponds to h1 indices d*24-23 .. d*24 - the 24 hours ENDING
        # at time d. Verified to float32 roundoff (5e-8) against the h7 daily mean;
        # using the naive d*24 .. d*24+23 shifts the comparison by a full day.
        acc = None
        for h in range(-23, 1):
            t = d * 24 + h
            V = h1v["V"].isel(time=t, lev=kidx).values.astype(np.float64)
            T = h1t["T"].isel(time=t, lev=kidx).values.astype(np.float64)
            e = eddy(V, T)
            acc = e if acc is None else acc + e
        B = acc / 24.0

        # --- C: daily-mean fields, single eddy flux
        Vd = h7v["V"].isel(time=d, lev=kidx).values.astype(np.float64)
        Td = h7t["T"].isel(time=d, lev=kidx).values.astype(np.float64)
        Cc = eddy(Vd, Td)

        vz_day = vthzm[d].values.astype(np.float64)[None]      # (1, ilev, lat)
        for i, pl in enumerate(TEST_P):
            # Interpolate the h6 field onto the EXACT h1/h7 mid-layer pressure, so
            # the comparison isolates the time-averaging effect rather than the
            # interface-vs-midpoint level offset.
            A = float(np.dot(
                C.interp_level(vz_day, ilev, p_lev[i] / 100.0, axis=1)[0][jsel], wj))
            b, c = band_of(B[i]), band_of(Cc[i])
            rows.append(dict(day=int(d), p_hPa=float(pl), model_h6=A,
                             hourly_h1=b, daily_h7=c,
                             ratio_hourly=b / A if A else np.nan,
                             ratio_daily=c / A if A else np.nan))
            print(f"{d:4d} {pl:7.1f} {A:10.2f} {b:10.2f} {c:10.2f} "
                  f"{b/A:7.3f} {c/A:7.3f}")

    for ds in (h1v, h1t, h7v, h7t):
        ds.close()

    rb = np.array([r["ratio_hourly"] for r in rows])
    rc = np.array([r["ratio_daily"] for r in rows])
    res = dict(
        days=DAYS, levels_hPa=TEST_P, band_deg=list(BAND),
        model_timestep_s=1800,
        rows=rows,
        median_ratio_hourly_to_model=float(np.median(rb)),
        median_ratio_daily_to_model=float(np.median(rc)),
        median_abs_error_hourly_pct=float(np.median(np.abs(rb - 1)) * 100),
        median_abs_error_daily_pct=float(np.median(np.abs(rc - 1)) * 100),
    )
    print("-" * 74)
    print(f"median |error| vs the model's own accumulation:")
    print(f"   B  rebuilt from 24 HOURLY means : "
          f"{res['median_abs_error_hourly_pct']:6.1f} %")
    print(f"   C  rebuilt from the DAILY mean  : "
          f"{res['median_abs_error_daily_pct']:6.1f} %")
    print(f"\nmedian ratio to model: hourly {np.median(rb):.3f}, "
          f"daily {np.median(rc):.3f}")
    print("A ratio below 1 means the rebuilt flux is too WEAK - covariance lost "
          "inside the averaging window.")

    with open(OUT, "w") as f:
        json.dump(res, f, indent=2)
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
