"""
01b - Why is the tropical w* coming out downward? Diagnose term by term.

The suspects, in order:
  S1  VTHzm is ALREADY the eddy flux v'theta', not the full product v*theta.
      If so, subtracting Vzm*THzm corrupts it.
  S2  d(theta)/dz sign/stencil on a coordinate that decreases with index.
  S3  the d/dphi stencil across the equator.

Prints latitude profiles at 70 hPa and vertical profiles in the deep tropics
for every term, under both interpretations of VTHzm.
"""

import numpy as np
import aide_val_common as C

SEG = "1996-2014"

grid = C.load_grid(SEG)
p, lat = grid["ilev"], grid["lat"]
z = C.log_p_height(p)
phi = np.deg2rad(lat)
cosphi = np.cos(phi)

vzm = np.nanmean(C.load_h6(SEG, "Vzm").values.astype(np.float64), axis=0)
wzm = np.nanmean(C.load_h6(SEG, "Wzm").values.astype(np.float64), axis=0)
thzm = np.nanmean(C.load_h6(SEG, "THzm").values.astype(np.float64), axis=0)
vthzm = np.nanmean(C.load_h6(SEG, "VTHzm").values.astype(np.float64), axis=0)

k70 = int(np.argmin(np.abs(p - 70.0)))
print(f"level used: index {k70}, p = {p[k70]:.3f} hPa\n")

# ---------------------------------------------------------------- S2 first
dthdz = np.gradient(thzm, z, axis=0)
print("S2  d(theta)/dz sanity  (must be POSITIVE everywhere in the stratosphere)")
print(f"    z is {'decreasing' if z[1] < z[0] else 'increasing'} with index; "
      f"z[0]={z[0]:.0f} m  z[-1]={z[-1]:.0f} m")
for pl in [10.0, 30.0, 70.0, 100.0]:
    k = int(np.argmin(np.abs(p - pl)))
    eq = C.band_mean(dthdz[k][None], lat, -10, 10, axis=-1)[0]
    print(f"    {p[k]:7.2f} hPa   theta={C.band_mean(thzm[k][None],lat,-10,10,axis=-1)[0]:7.1f} K"
          f"   dtheta/dz(10S-10N) = {eq*1000:+8.3f} K/km")
print()

# ---------------------------------------------------------------- S1
print("S1  is VTHzm the full product or already the eddy flux?")
print("    latitude profile at 70 hPa, annual mean")
print(f"    {'lat':>7} {'VTHzm':>10} {'Vzm*THzm':>10} {'difference':>11}   (K m/s)")
for la in [-70, -50, -30, -10, 0, 10, 30, 50, 70]:
    j = int(np.argmin(np.abs(lat - la)))
    full = vthzm[k70, j]
    mean = vzm[k70, j] * thzm[k70, j]
    print(f"    {lat[j]:7.1f} {full:10.2f} {mean:10.2f} {full-mean:11.2f}")
print("    -> if column 1 ~ column 2, VTHzm is the FULL product (subtract).")
print("       if column 1 is small and column 2 is large, VTHzm is already eddy.\n")


def wstar_from(vpt, tag):
    psi = vpt / dthdz
    num = psi * cosphi[None, :]
    dnum = np.gradient(num, phi, axis=1) / (C.A_EARTH * cosphi[None, :])
    ws = wzm + dnum
    trop_w = C.band_mean(ws[k70][None], lat, -10, 10, axis=-1)[0]
    trop_eddy = C.band_mean(dnum[k70][None], lat, -10, 10, axis=-1)[0]
    print(f"  {tag}")
    print(f"    10S-10N @ {p[k70]:.1f} hPa:  w_bar = {wzm[k70].mean()*0 + C.band_mean(wzm[k70][None],lat,-10,10,axis=-1)[0]*1e3:+7.4f}"
          f"  eddy = {trop_eddy*1e3:+7.4f}  w* = {trop_w*1e3:+7.4f} mm/s")
    print(f"    {'lat':>7} {'psi(m2/s)':>12} {'eddy term':>12} {'w*':>10}  (mm/s)")
    for la in [-60, -40, -20, -10, 0, 10, 20, 40, 60]:
        j = int(np.argmin(np.abs(lat - la)))
        print(f"    {lat[j]:7.1f} {psi[k70,j]:12.1f} {dnum[k70,j]*1e3:12.4f} {ws[k70,j]*1e3:10.4f}")
    print()
    return ws


print("S3  w* under each interpretation")
ws_sub = wstar_from(vthzm - vzm * thzm, "(i)  v'theta' = VTHzm - Vzm*THzm   [current assumption]")
ws_raw = wstar_from(vthzm, "(ii) v'theta' = VTHzm directly      [alternative]")

print("Expected from the literature (Butchart 2014 review, CCMVal-2 multi-model):")
print("  tropical w* at 70 hPa, annual mean, ~0.2-0.4 mm/s UPWARD;")
print("  turnaround latitudes near +/-30 deg; midlatitude w* downward.")
