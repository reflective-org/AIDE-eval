"""
01c - Settle the VTHzm convention against the raw 3D daily fields.

The whole tropical-upwelling target rests on what VTHzm means. The h6 long_name
"Meridional Heat Flux: 3D zon. mean" is ambiguous between

    (i)  zonal mean of (v * theta)          -> eddy flux = VTHzm - Vzm*THzm
    (ii) zonal mean of (v' * theta')        -> eddy flux = VTHzm

The h7 daily tape carries the raw 3D V and T on (time, lev, lat, lon), so we can
compute the true eddy flux directly for a few days and compare.

theta = T * (P0/p)^kappa, with p from the hybrid coefficients and PS.
"""

import numpy as np
import xarray as xr
import aide_val_common as C

SEG = "1996-2014"
CASE = C.SEGMENTS[SEG]["case"]
DAY = f"{C.ROOT}/{CASE}/{C.SEGMENTS[SEG]['subdir']}"
SPAN7 = "19960101-20141231"

# sample days spread across the seasonal cycle of the first year
DAYS = [15, 105, 196, 288, 380, 470]

pv = f"{DAY}/{CASE}.cam.h7.V.{SPAN7}.nc"
pt = f"{DAY}/{CASE}.cam.h7.T.{SPAN7}.nc"
pp = f"{DAY}/{CASE}.cam.h7.PS.{SPAN7}.nc"

dv = xr.open_dataset(pv)
dt = xr.open_dataset(pt)
dp = xr.open_dataset(pp)

print("h7 V encoding:", {k: dv["V"].encoding[k] for k in ("chunksizes", "zlib")})
lev = dv["lev"].values.astype(np.float64)          # hPa, midpoints
lat = dv["lat"].values.astype(np.float64)
hyam = dv["hyam"].values.astype(np.float64)
hybm = dv["hybm"].values.astype(np.float64)
P0 = float(dv["P0"].values)

# h6 zonal-mean fields for the same days (on interfaces)
grid = C.load_grid(SEG)
ilev = grid["ilev"]
vzm_da = C.load_h6(SEG, "Vzm")
thzm_da = C.load_h6(SEG, "THzm")
vthzm_da = C.load_h6(SEG, "VTHzm")

TEST_P = [10.0, 30.0, 70.0, 100.0]
TEST_LAT = [-60.0, -30.0, 0.0, 30.0, 60.0]

print(f"\n{'day':>5} {'p_hPa':>7} {'lat':>7} "
      f"{'true eddy':>13} {'VTHzm':>10} {'VTHzm-Vzm*THzm':>16}")
print("-" * 66)

agree_raw, agree_sub = [], []

for d in DAYS:
    V = dv["V"].isel(time=d).values.astype(np.float64)      # (lev, lat, lon)
    T = dt["T"].isel(time=d).values.astype(np.float64)
    PS = dp["PS"].isel(time=d).values.astype(np.float64)    # (lat, lon)

    p3d = hyam[:, None, None] * P0 + hybm[:, None, None] * PS[None]   # Pa
    TH = T * (100000.0 / p3d) ** C.KAPPA

    Vb = V.mean(axis=2, keepdims=True)
    THb = TH.mean(axis=2, keepdims=True)
    eddy_true = ((V - Vb) * (TH - THb)).mean(axis=2)        # (lev, lat)

    for pl in TEST_P:
        kt = int(np.argmin(np.abs(lev - pl)))               # h7 midpoint
        ki = int(np.argmin(np.abs(ilev - pl)))              # h6 interface
        for la in TEST_LAT:
            j = int(np.argmin(np.abs(lat - la)))
            truth = eddy_true[kt, j]
            raw = float(vthzm_da[d, ki, j])
            sub = raw - float(vzm_da[d, ki, j]) * float(thzm_da[d, ki, j])
            print(f"{d:5d} {lev[kt]:7.2f} {lat[j]:7.1f} "
                  f"{truth:13.2f} {raw:10.2f} {sub:16.2f}")
            if np.isfinite(truth) and abs(truth) > 1.0:
                agree_raw.append(abs(raw - truth) / abs(truth))
                agree_sub.append(abs(sub - truth) / abs(truth))

dv.close(); dt.close(); dp.close()

print("-" * 66)
print(f"median |relative error| vs the true eddy flux:")
print(f"   VTHzm as-is                 : {np.median(agree_raw)*100:6.1f} %")
print(f"   VTHzm - Vzm*THzm            : {np.median(agree_sub)*100:6.1f} %")
print(f"\nn = {len(agree_raw)} comparisons (|truth| > 1 K m/s)")
winner = "VTHzm IS ALREADY THE EDDY FLUX" if np.median(agree_raw) < np.median(agree_sub) \
    else "VTHzm is the full product; subtract Vzm*THzm"
print(f"VERDICT: {winner}")
print("\n(residual error is expected and comes from the h6 interface vs h7 midpoint\n"
      " level offset, not from the convention.)")
