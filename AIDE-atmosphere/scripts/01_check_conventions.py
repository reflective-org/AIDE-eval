"""
01 - Establish what the h6 TEM fields actually are, before trusting any target.

Three things are checked, because each one silently changes the tropical
upwelling number by 10-40% if assumed wrong:

  (a) Is Wzm the log-pressure vertical velocity w = -H*omega/p (H = 7 km), or a
      geometric one w = -omega*R*T/(g*p) (local scale height ~6.1 km at 210 K)?
      Determined by mass continuity: in log-pressure coordinates
          (1/(a cos.phi)) d(v_bar cos.phi)/d.phi + (1/rho0) d(rho0 w_bar)/dz = 0
      holds exactly with rho0 ~ exp(-z/H) only for the log-pressure w. We
      integrate the divergence downward from the model top and regress the
      implied w against Wzm. Slope ~1 => log-pressure. Slope ~H/H_local => not.

  (b) Does the eddy heat flux reconstruction VTHzm - Vzm*THzm behave like a
      physical v'theta' (poleward in the winter stratosphere, right magnitude)?

  (c) Do the two independent routes to w* agree - the direct AHL 3.5.1 form and
      the one obtained by integrating v* through continuity? Disagreement means
      the derivative stencils or the level spacing are the limiting error, and
      that error has to be carried into the tolerance budget.

Output: output/01_conventions.json and a printed report.
"""

import json
import os
import numpy as np

import aide_val_common as C

SEG = "1996-2014"
OUT = os.path.join(C.OUTDIR, "01_conventions.json")


def main():
    print(f"loading h6 TEM fields, segment {SEG} ...", flush=True)
    grid = C.load_grid(SEG)
    p = grid["ilev"]                 # hPa, increasing downward
    lat = grid["lat"]
    z = C.log_p_height(p)

    vzm = C.load_h6(SEG, "Vzm").values.astype(np.float64)
    wzm = C.load_h6(SEG, "Wzm").values.astype(np.float64)
    thzm = C.load_h6(SEG, "THzm").values.astype(np.float64)
    vthzm = C.load_h6(SEG, "VTHzm").values.astype(np.float64)
    print(f"  shapes {vzm.shape}, levels {p[0]:.3f}-{p[-1]:.1f} hPa", flush=True)

    report = {"segment": SEG, "n_days": int(vzm.shape[0]),
              "p_hpa_range": [float(p[0]), float(p[-1])]}

    # ---------------------------------------------------------------- (a)
    # Time-mean fields; continuity is a statement about the mean circulation and
    # the daily fields are far noisier than the stencil error we are probing.
    v_bar = np.nanmean(vzm, axis=0)          # (lev, lat)
    w_bar = np.nanmean(wzm, axis=0)

    phi = np.deg2rad(lat)
    cosphi = np.cos(phi)
    rho0 = np.exp(-z / C.H_SCALE)

    # horizontal mass divergence, 1/(a cos) d(v cos)/dphi
    div = np.gradient(v_bar * cosphi[None, :], phi, axis=1) / (C.A_EARTH * cosphi[None, :])

    # integrate rho0*w downward from the top: d(rho0 w)/dz = -rho0 * div
    integrand = -rho0[:, None] * div
    rho0w = np.zeros_like(integrand)
    # z decreases with index (p increases), so integrate with -dz
    for k in range(1, len(z)):
        dz = z[k] - z[k - 1]
        rho0w[k] = rho0w[k - 1] + 0.5 * (integrand[k] + integrand[k - 1]) * dz
    w_cont = rho0w / rho0[:, None]

    # regress w_cont on w_bar where the signal is meaningful: 5-150 hPa, |lat|<60
    ksel = (p >= 5.0) & (p <= 150.0)
    jsel = np.abs(lat) < 60.0
    x = w_bar[np.ix_(ksel, jsel)].ravel()
    y = w_cont[np.ix_(ksel, jsel)].ravel()
    good = np.isfinite(x) & np.isfinite(y)
    slope = float(np.sum(x[good] * y[good]) / np.sum(x[good] ** 2))
    corr = float(np.corrcoef(x[good], y[good])[0, 1])

    # what slope would each convention predict?
    t_typ = 215.0                                   # K, lower stratosphere
    h_local = C.RGAS * t_typ / C.G0
    report["continuity_check"] = dict(
        slope_wcont_on_Wzm=slope,
        correlation=corr,
        n_points=int(good.sum()),
        expected_slope_if_logp=1.0,
        expected_slope_if_geometric=float(C.H_SCALE / h_local),
        h_local_m_at_215K=float(h_local),
    )
    print("\n(a) Wzm convention via mass continuity")
    print(f"    slope(w_continuity on Wzm) = {slope:.3f}   r = {corr:.3f}   "
          f"(n={good.sum()})")
    print(f"    expected 1.000 if log-pressure w = -H.omega/p  (H = 7000 m)")
    print(f"    expected {C.H_SCALE/h_local:.3f} if geometric w = -omega.R.T/(g.p) "
          f"(H_local = {h_local:.0f} m at 215 K)")
    verdict = "log-pressure" if abs(slope - 1.0) < abs(slope - C.H_SCALE / h_local) \
        else "geometric"
    report["continuity_check"]["verdict"] = verdict
    print(f"    -> Wzm is closer to the {verdict.upper()} convention")

    # ---------------------------------------------------------------- (b)
    # VTHzm is already the eddy flux v'theta' - proved in 01c against the raw
    # 3D h7 fields (15% error as-is vs 1087% if Vzm*THzm is subtracted).
    vpt = vthzm
    vpt_bar = np.nanmean(vpt, axis=0)
    # DJF, 45-75N, 10-100 hPa: should be strongly poleward (positive) in NH winter
    yr, mo, dy = C.time_index(SEG)
    djf = C.djf_mask(mo)
    jja = np.isin(mo, [6, 7, 8])
    ksel2 = (p >= 10.0) & (p <= 100.0)

    def band(field_t, j0, j1):
        jj = (lat >= j0) & (lat <= j1)
        return float(np.nanmean(field_t[np.ix_(ksel2, jj)]))

    vpt_djf = np.nanmean(vpt[djf], axis=0)
    vpt_jja = np.nanmean(vpt[jja], axis=0)
    report["eddy_heat_flux"] = dict(
        djf_45_75N_Kms=band(vpt_djf, 45, 75),
        jja_45_75N_Kms=band(vpt_djf * 0 + vpt_jja, 45, 75),
        djf_75_45S_Kms=band(vpt_djf, -75, -45),
        jja_75_45S_Kms=band(vpt_djf * 0 + vpt_jja, -75, -45),
        annual_45_75N_Kms=band(vpt_bar, 45, 75),
    )
    print("\n(b) eddy heat flux v-theta = VTHzm (already eddy, see 01c), 10-100 hPa mean")
    print(f"    DJF 45-75N  {band(vpt_djf,45,75):+8.2f} K m/s   "
          f"(expect strongly POSITIVE, poleward, NH winter waves)")
    print(f"    JJA 45-75N  {band(vpt_jja,45,75):+8.2f} K m/s   (expect weak)")
    print(f"    JJA 45-75S  {band(vpt_jja,-75,-45):+8.2f} K m/s  "
          f"(expect NEGATIVE = poleward in SH, weaker than NH DJF)")
    print(f"    DJF 45-75S  {band(vpt_djf,-75,-45):+8.2f} K m/s  (expect weak)")

    # ---------------------------------------------------------------- (c)
    vstar, wstar, psi = C.tem_residual(vzm, wzm, thzm, vthzm, p, lat)
    vstar_bar = np.nanmean(vstar, axis=0)
    wstar_bar = np.nanmean(wstar, axis=0)

    # w* from v* through continuity
    divs = np.gradient(vstar_bar * cosphi[None, :], phi, axis=1) / (C.A_EARTH * cosphi[None, :])
    integ = -rho0[:, None] * divs
    rho0ws = np.zeros_like(integ)
    for k in range(1, len(z)):
        dz = z[k] - z[k - 1]
        rho0ws[k] = rho0ws[k - 1] + 0.5 * (integ[k] + integ[k - 1]) * dz
    wstar_cont = rho0ws / rho0[:, None]

    w70_direct = float(C.band_mean(C.interp_level(wstar_bar[None], p, 70.0, axis=1)[0],
                                   lat, -10, 10, axis=-1))
    w70_cont = float(C.band_mean(C.interp_level(wstar_cont[None], p, 70.0, axis=1)[0],
                                 lat, -10, 10, axis=-1))
    w70_eulerian = float(C.band_mean(C.interp_level(w_bar[None], p, 70.0, axis=1)[0],
                                     lat, -10, 10, axis=-1))

    report["wstar_cross_check"] = dict(
        wstar_70hPa_10S10N_direct_mm_s=w70_direct * 1e3,
        wstar_70hPa_10S10N_continuity_mm_s=w70_cont * 1e3,
        eulerian_wbar_70hPa_10S10N_mm_s=w70_eulerian * 1e3,
        relative_difference_pct=100.0 * abs(w70_direct - w70_cont) / abs(w70_direct),
    )
    print("\n(c) tropical mean vertical velocity at 70 hPa, 10S-10N, 19-yr mean")
    print(f"    Eulerian  w_bar     {w70_eulerian*1e3:8.4f} mm/s")
    print(f"    residual  w* direct {w70_direct*1e3:8.4f} mm/s   (AHL 3.5.1)")
    print(f"    residual  w* contin {w70_cont*1e3:8.4f} mm/s   (from v*)")
    print(f"    two w* routes differ by "
          f"{100*abs(w70_direct-w70_cont)/abs(w70_direct):.1f} %  "
          f"<- this is the METHOD uncertainty floor for any w* target")

    ts, tn = C.turnaround_latitudes(
        C.interp_level(wstar_bar[None], p, 70.0, axis=1)[0], lat)
    report["turnaround_latitudes_70hPa"] = [float(ts), float(tn)]
    print(f"\n    turnaround latitudes at 70 hPa: {ts:.1f} to {tn:.1f} deg "
          f"(w* = 0 crossings)")

    os.makedirs(C.OUTDIR, exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(report, f, indent=2)
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
