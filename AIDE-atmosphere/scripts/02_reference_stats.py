"""
02 - Reference statistics from CESM2-WACCM6 for the two validation targets.

Computes, per segment and for the pooled 45-year record:

  TROPICAL UPWELLING
    w* at 70 hPa, cos-weighted 10S-10N            (TEM residual, AHL 3.5.1)
    w* at 70 hPa averaged between turnaround latitudes
    tropical upward mass flux across 70 hPa       (Butchart/CCMVal convention)
    - annual mean, 12-month climatology, daily/monthly/interannual sigma

  POLAR VORTEX
    u at 60N and 60S, 10 hPa                      (interpolated to exact lat/p)
    - DJF (NH) / JJA (SH) mean, monthly climatology
    - daily sigma within season, interannual sigma of the seasonal mean
    - 5th/95th percentiles of daily values in season
    - SSW frequency, Charlton-Polvani (2007) major midwinter warmings
    polar cap temperature 60-90N, 10-50 hPa, from THzm

Every sigma is reported with its own sampling uncertainty, because that sets the
floor on how tight a target can honestly be.

Output: output/02_reference_stats.json  and a printed report.
"""

import json
import os
import numpy as np

import aide_val_common as C

OUT = os.path.join(C.OUTDIR, "02_reference_stats.json")


def process_segment(seg):
    print(f"\n{'='*78}\nSEGMENT {seg}\n{'='*78}", flush=True)
    grid = C.load_grid(seg)
    p, lat = grid["ilev"], grid["lat"]
    yr, mo, dy = C.time_index(seg)

    vzm = C.load_h6(seg, "Vzm").values.astype(np.float64)
    wzm = C.load_h6(seg, "Wzm").values.astype(np.float64)
    thzm = C.load_h6(seg, "THzm").values.astype(np.float64)
    vthzm = C.load_h6(seg, "VTHzm").values.astype(np.float64)
    uzm = C.load_h6(seg, "Uzm").values.astype(np.float64)
    print(f"loaded {vzm.shape[0]} days x {vzm.shape[1]} levels x {vzm.shape[2]} lats",
          flush=True)

    vstar, wstar, psi = C.tem_residual(vzm, wzm, thzm, vthzm, p, lat)

    r = {"segment": seg, "n_days": int(len(yr)),
         "years": [int(yr.min()), int(yr.max())]}

    # ---------------------------------------------------------- upwelling
    w70 = C.interp_level(wstar, p, 70.0, axis=1)             # (time, lat)
    w70_trop = C.band_mean(w70, lat, -10.0, 10.0, axis=-1)   # (time,)

    w70_clim = np.nanmean(w70, axis=0)
    ts, tn = C.turnaround_latitudes(w70_clim, lat)
    w70_turn = C.band_mean(w70, lat, ts, tn, axis=-1)

    mflux = C.tropical_upward_mass_flux(wstar, p, lat, 70.0)

    # same at 100 hPa and 10 hPa, for the profile
    extra = {}
    for pl in (100.0, 50.0, 30.0, 10.0):
        wl = C.interp_level(wstar, p, pl, axis=1)
        extra[f"w_star_{pl:g}hPa_10S10N_mm_s"] = {
            "mean": float(np.nanmean(C.band_mean(wl, lat, -10, 10, axis=-1)) * 1e3)}

    def upwelling_block(x, name, scale, unit):
        yrs_a, ann = C.annual_means(x, yr)
        clim = C.monthly_climatology(x, mo)
        sig_ia = float(np.std(ann, ddof=1))
        d = dict(
            units=unit,
            annual_mean=float(np.nanmean(x)) * scale,
            monthly_climatology=[float(v) * scale for v in clim],
            seasonal_amplitude=float(clim.max() - clim.min()) * scale,
            month_of_max=int(np.argmax(clim) + 1),
            month_of_min=int(np.argmin(clim) + 1),
            sigma_daily=float(np.nanstd(x, ddof=1)) * scale,
            sigma_interannual=sig_ia * scale,
            sigma_interannual_rel_pct=100.0 * sig_ia / abs(np.nanmean(x)),
            sigma_interannual_sampling_err=float(
                C.sigma_sampling_error(sig_ia, len(ann))) * scale,
            n_years=int(len(ann)),
            annual_means=[float(v) * scale for v in ann],
            annual_years=[int(v) for v in yrs_a],
        )
        lo, hi = C.bootstrap_ci(ann, lambda a: np.std(a, ddof=1))
        d["sigma_interannual_ci95"] = [float(lo) * scale, float(hi) * scale]
        lo, hi = C.bootstrap_ci(ann, np.mean)
        d["annual_mean_ci95"] = [float(lo) * scale, float(hi) * scale]
        return d

    r["tropical_upwelling"] = {
        "turnaround_latitudes_70hPa": [float(ts), float(tn)],
        "w_star_70hPa_10S10N": upwelling_block(w70_trop, "w70t", 1e3, "mm/s"),
        "w_star_70hPa_turnaround": upwelling_block(w70_turn, "w70n", 1e3, "mm/s"),
        "upward_mass_flux_70hPa": upwelling_block(mflux, "mf", 1e-9, "1e9 kg/s"),
        "profile": extra,
    }

    u = r["tropical_upwelling"]["w_star_70hPa_10S10N"]
    print(f"\nTROPICAL UPWELLING  w* 70 hPa, 10S-10N")
    print(f"  annual mean            {u['annual_mean']:.4f} mm/s  "
          f"(95% CI {u['annual_mean_ci95'][0]:.4f} - {u['annual_mean_ci95'][1]:.4f})")
    print(f"  seasonal amplitude     {u['seasonal_amplitude']:.4f} mm/s  "
          f"(max month {u['month_of_max']}, min month {u['month_of_min']})")
    print(f"  sigma interannual      {u['sigma_interannual']:.4f} mm/s  "
          f"= {u['sigma_interannual_rel_pct']:.1f} % of the mean  "
          f"(+/- {u['sigma_interannual_sampling_err']:.4f} sampling)")
    print(f"  sigma daily            {u['sigma_daily']:.4f} mm/s")
    print(f"  turnaround latitudes   {ts:.1f} to {tn:.1f} deg")
    m = r["tropical_upwelling"]["upward_mass_flux_70hPa"]
    print(f"  upward mass flux       {m['annual_mean']:.3f} x1e9 kg/s  "
          f"(sigma_ia {m['sigma_interannual']:.3f} = {m['sigma_interannual_rel_pct']:.1f} %)")
    print(f"  w* profile 10S-10N     " + "  ".join(
        f"{k.split('_')[2]}={v['mean']:.3f}" for k, v in extra.items()))

    # ---------------------------------------------------------- vortex
    r["polar_vortex"] = {}
    for hemi, latv, months, seasname in (("NH", 60.0, [12, 1, 2], "DJF"),
                                         ("SH", -60.0, [6, 7, 8], "JJA")):
        u10 = C.interp_level(uzm, p, 10.0, axis=1)              # (time, lat)
        u60 = C.interp_lat(u10, lat, latv, axis=-1)             # (time,)
        sgn = 1.0 if hemi == "NH" else -1.0
        u60s = sgn * u60          # sign-flipped so westerly is positive in both

        labfn = C.season_year if hemi == "NH" else None
        yrs_s, seas = C.seasonal_means(u60s, yr, mo, months, labfn)
        sel = np.isin(mo, months)
        clim = C.monthly_climatology(u60s, mo)

        sig_ia = float(np.std(seas, ddof=1))
        sig_daily = float(np.nanstd(u60s[sel], ddof=1))
        lo_ia, hi_ia = C.bootstrap_ci(seas, lambda a: np.std(a, ddof=1))
        lo_m, hi_m = C.bootstrap_ci(seas, np.mean)

        ssw = C.detect_ssw_cp07(u60s, yr, mo, dy, hemi)
        n_seasons = len(np.unique(
            (C.season_year(yr, mo) if hemi == "NH" else yr)[
                np.isin(mo, [11, 12, 1, 2, 3] if hemi == "NH" else [5, 6, 7, 8, 9])]))
        # drop edge seasons that the record cannot fully cover
        n_seasons = max(n_seasons - 1, 1)
        freq = len(ssw) / n_seasons

        blk = dict(
            latitude=latv, level_hPa=10.0, season=seasname,
            units="m/s",
            sign_convention=("as stored" if hemi == "NH"
                             else "sign-flipped so westerly > 0"),
            season_mean=float(np.mean(seas)),
            season_mean_ci95=[float(lo_m), float(hi_m)],
            monthly_climatology=[float(v) for v in clim],
            sigma_daily_in_season=sig_daily,
            sigma_interannual_of_season_mean=sig_ia,
            sigma_interannual_ci95=[float(lo_ia), float(hi_ia)],
            sigma_interannual_sampling_err=float(
                C.sigma_sampling_error(sig_ia, len(seas))),
            pct05_daily_in_season=float(np.nanpercentile(u60s[sel], 5)),
            pct95_daily_in_season=float(np.nanpercentile(u60s[sel], 95)),
            n_seasons=int(len(seas)),
            season_means=[float(v) for v in seas],
            season_years=[int(v) for v in yrs_s],
            ssw_count=len(ssw),
            ssw_seasons_analysed=int(n_seasons),
            ssw_freq_per_winter=float(freq),
            ssw_events=ssw,
        )
        # Poisson 95% interval on the count, then converted to a rate
        k = len(ssw)
        from scipy.stats import chi2
        lo_k = 0.0 if k == 0 else chi2.ppf(0.025, 2 * k) / 2.0
        hi_k = chi2.ppf(0.975, 2 * (k + 1)) / 2.0
        blk["ssw_freq_ci95_poisson"] = [lo_k / n_seasons, hi_k / n_seasons]
        r["polar_vortex"][hemi] = blk

        print(f"\nPOLAR VORTEX {hemi}  u at {latv:.0f} deg, 10 hPa  ({seasname})")
        print(f"  {seasname} mean            {blk['season_mean']:.2f} m/s  "
              f"(95% CI {lo_m:.2f} - {hi_m:.2f})")
        print(f"  sigma interannual      {sig_ia:.2f} m/s  "
              f"(+/- {blk['sigma_interannual_sampling_err']:.2f} sampling, "
              f"95% CI {lo_ia:.2f} - {hi_ia:.2f})")
        print(f"  sigma daily in season  {sig_daily:.2f} m/s")
        print(f"  5th / 95th pct daily   {blk['pct05_daily_in_season']:.2f} / "
              f"{blk['pct95_daily_in_season']:.2f} m/s")
        print(f"  SSW (CP07 major)       {k} events in {n_seasons} winters "
              f"= {freq:.3f}/winter  "
              f"(Poisson 95% CI {lo_k/n_seasons:.3f} - {hi_k/n_seasons:.3f})")

    # ---------------------------------------------------------- polar cap T
    # T = theta * (p/p0)^kappa, layer-mean over 10-50 hPa, 60-90 deg
    ksel = (p >= 10.0) & (p <= 50.0)
    tzm = thzm[:, ksel, :] * (p[ksel][None, :, None] / C.P_REF_HPA) ** C.KAPPA
    logp = np.log(p[ksel])
    wts = np.gradient(logp)
    wts = wts / wts.sum()
    t_lev = np.tensordot(tzm, wts, axes=([1], [0]))          # (time, lat)

    r["polar_cap_T"] = {}
    for hemi, l0, l1, months, seasname in (("NH", 60.0, 90.0, [12, 1, 2], "DJF"),
                                           ("SH", -90.0, -60.0, [6, 7, 8], "JJA")):
        tcap = C.band_mean(t_lev, lat, l0, l1, axis=-1)
        labfn = C.season_year if hemi == "NH" else None
        yrs_s, seas = C.seasonal_means(tcap, yr, mo, months, labfn)
        sel = np.isin(mo, months)
        sig_ia = float(np.std(seas, ddof=1))
        lo_m, hi_m = C.bootstrap_ci(seas, np.mean)
        r["polar_cap_T"][hemi] = dict(
            lat_band=[l0, l1], p_band_hPa=[10.0, 50.0], season=seasname, units="K",
            annual_mean=float(np.nanmean(tcap)),
            season_mean=float(np.mean(seas)),
            season_mean_ci95=[float(lo_m), float(hi_m)],
            monthly_climatology=[float(v) for v in C.monthly_climatology(tcap, mo)],
            sigma_daily_in_season=float(np.nanstd(tcap[sel], ddof=1)),
            sigma_interannual_of_season_mean=sig_ia,
            sigma_interannual_sampling_err=float(
                C.sigma_sampling_error(sig_ia, len(seas))),
            n_seasons=int(len(seas)),
            season_means=[float(v) for v in seas],
        )
        b = r["polar_cap_T"][hemi]
        print(f"\nPOLAR CAP T {hemi}  {abs(l0):.0f}-{abs(l1):.0f} deg, 10-50 hPa")
        print(f"  {seasname} mean            {b['season_mean']:.2f} K  "
              f"(95% CI {lo_m:.2f} - {hi_m:.2f})")
        print(f"  sigma interannual      {sig_ia:.2f} K  "
              f"(+/- {b['sigma_interannual_sampling_err']:.2f} sampling)")
        print(f"  sigma daily in season  {b['sigma_daily_in_season']:.2f} K")

    return r


def main():
    os.makedirs(C.OUTDIR, exist_ok=True)
    allr = {}
    for seg in ("1996-2014", "1970-1995"):
        allr[seg] = process_segment(seg)

    # ------------------------------------------------------------ pooled
    print(f"\n{'='*78}\nPOOLED 1970-2014 (both segments, 45 years)\n{'='*78}")
    pooled = {}

    def pool(path, key):
        a = allr["1970-1995"]
        b = allr["1996-2014"]
        for k in path:
            a, b = a[k], b[k]
        vals = np.array(a[key] + b[key])
        return vals

    for name, path, unit in (
        ("w_star_70hPa_10S10N", ["tropical_upwelling", "w_star_70hPa_10S10N"], "mm/s"),
        ("upward_mass_flux_70hPa", ["tropical_upwelling", "upward_mass_flux_70hPa"],
         "1e9 kg/s"),
    ):
        vals = pool(path, "annual_means")
        sig = float(np.std(vals, ddof=1))
        lo, hi = C.bootstrap_ci(vals, lambda a: np.std(a, ddof=1))
        pooled[name] = dict(units=unit, n_years=len(vals),
                            mean=float(vals.mean()), sigma_interannual=sig,
                            sigma_interannual_ci95=[float(lo), float(hi)],
                            sigma_interannual_sampling_err=float(
                                C.sigma_sampling_error(sig, len(vals))),
                            sigma_interannual_rel_pct=100.0 * sig / abs(vals.mean()))
        print(f"  {name:26s} mean {vals.mean():8.4f} {unit}   "
              f"sigma_ia {sig:.4f} +/- "
              f"{C.sigma_sampling_error(sig, len(vals)):.4f}  (n={len(vals)})")

    for hemi in ("NH", "SH"):
        vals = pool(["polar_vortex", hemi], "season_means")
        sig = float(np.std(vals, ddof=1))
        lo, hi = C.bootstrap_ci(vals, lambda a: np.std(a, ddof=1))
        k = (allr["1970-1995"]["polar_vortex"][hemi]["ssw_count"]
             + allr["1996-2014"]["polar_vortex"][hemi]["ssw_count"])
        n = (allr["1970-1995"]["polar_vortex"][hemi]["ssw_seasons_analysed"]
             + allr["1996-2014"]["polar_vortex"][hemi]["ssw_seasons_analysed"])
        from scipy.stats import chi2
        lo_k = 0.0 if k == 0 else chi2.ppf(0.025, 2 * k) / 2.0
        hi_k = chi2.ppf(0.975, 2 * (k + 1)) / 2.0
        pooled[f"vortex_{hemi}"] = dict(
            units="m/s", n_seasons=len(vals), mean=float(vals.mean()),
            sigma_interannual=sig, sigma_interannual_ci95=[float(lo), float(hi)],
            sigma_interannual_sampling_err=float(C.sigma_sampling_error(sig, len(vals))),
            ssw_count=int(k), ssw_seasons=int(n), ssw_freq_per_winter=k / n,
            ssw_freq_ci95_poisson=[lo_k / n, hi_k / n])
        print(f"  vortex_{hemi:22s} mean {vals.mean():8.2f} m/s      "
              f"sigma_ia {sig:.2f} +/- {C.sigma_sampling_error(sig,len(vals)):.2f}  "
              f"(n={len(vals)})")
        print(f"  {'':26s} SSW {k}/{n} = {k/n:.3f}/winter "
              f"(Poisson 95% CI {lo_k/n:.3f}-{hi_k/n:.3f})")

    for hemi in ("NH", "SH"):
        vals = pool(["polar_cap_T", hemi], "season_means")
        sig = float(np.std(vals, ddof=1))
        pooled[f"polar_cap_T_{hemi}"] = dict(
            units="K", n_seasons=len(vals), mean=float(vals.mean()),
            sigma_interannual=sig,
            sigma_interannual_sampling_err=float(C.sigma_sampling_error(sig, len(vals))))
        print(f"  polar_cap_T_{hemi:17s} mean {vals.mean():8.2f} K        "
              f"sigma_ia {sig:.2f}  (n={len(vals)})")

    allr["pooled_1970_2014"] = pooled

    with open(OUT, "w") as f:
        json.dump(allr, f, indent=2)
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
