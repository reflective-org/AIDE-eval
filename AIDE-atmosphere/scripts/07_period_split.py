"""
07 - Set the targets on 1996-2014, then verify them out-of-sample on 1970-1995.

The targets so far were anchored on the pooled 45-year record, which is circular:
the data that sets a threshold cannot also test it. This script does the honest
version.

  TRAIN  1996-2014 (19 yr)  -> sigma, mean, slopes -> thresholds
  TEST   1970-1995 (26 yr)  -> treated as if it were a model rollout and scored
                              against those thresholds

If the CESM 1970-1995 period fails a target derived from CESM 1996-2014, the
threshold is too tight for a *different sample of the same model* to meet, and no
emulator should be held to it. A target that CESM itself cannot pass is not a
target, it is a trap.

Everything is computed twice - raw and detrended within each period - because the
forced BDC trend is exactly the thing that makes the two periods differ.

Output: output/07_period_split.json
"""

import json
import os
import numpy as np
from scipy import stats

import aide_val_common as C

OUT = os.path.join(C.OUTDIR, "07_period_split.json")
TRAIN, TEST = "1996-2014", "1970-1995"

PROFILE_LEVELS, DJF_PCTL = C.PROFILE_LEVELS, C.DJF_PCTL


def monthly_by_year(daily, yr, mo, min_days=20):
    """Per-year 12-month means of each daily series.

    A threshold on the seasonal cycle needs an interannual sigma, so the
    climatology has to be resolved per year rather than pooled across the record.
    Years with a month thinner than min_days are dropped whole.
    """
    years = np.array([int(y) for y in np.unique(yr)
                      if min(np.sum((yr == y) & (mo == k)) for k in range(1, 13))
                      >= min_days])
    rows = {k: np.array([[np.nanmean(x[(yr == y) & (mo == m)]) for m in range(1, 13)]
                         for y in years])
            for k, x in daily.items()}
    return years, rows


def profile_by_year(wstar, p, lat, yr):
    """Tropical 10S-10N w* at each profile level, calendar-year means, mm/s."""
    years = np.array([int(y) for y in np.unique(yr) if np.sum(yr == y) >= 350])
    out = {}
    for pl in PROFILE_LEVELS:
        b = C.band_mean(C.interp_level(wstar, p, pl, axis=1),
                        lat, -10, 10, axis=-1) * 1e3
        out[f"{pl:g}"] = np.array([np.nanmean(b[yr == y]) for y in years])
    return years, out


def winter_percentiles(x, yr, mo, min_days=81):
    """Per-winter DJF percentiles of a daily series.

    Computed within each winter rather than over the pooled days, because winters
    are independent: the D5 rule then applies with n = winters and the 14-day
    decorrelation time of the daily series never enters.
    """
    sel = C.djf_mask(mo)
    lab = C.season_year(yr, mo)
    years, rows = [], []
    for y in np.unique(lab[sel]):
        s = sel & (lab == y)
        if s.sum() >= min_days:
            years.append(int(y))
            rows.append([float(np.nanpercentile(x[s], q)) for q in DJF_PCTL])
    return np.array(years), np.array(rows)


def diagnostics_for(seg):
    """Every scalar series this validation suite needs, for one period."""
    grid = C.load_grid(seg)
    p, lat = grid["ilev"], grid["lat"]
    yr, mo, dy = C.time_index(seg)

    f = {v: C.load_h6(seg, v).values.astype(np.float64)
         for v in ("Uzm", "Vzm", "Wzm", "THzm", "VTHzm")}
    _, wstar, _ = C.tem_residual(f["Vzm"], f["Wzm"], f["THzm"], f["VTHzm"], p, lat)

    w70 = C.band_mean(C.interp_level(wstar, p, 70.0, axis=1), lat, -10, 10, axis=-1)
    mflux = C.tropical_upward_mass_flux(wstar, p, lat, 70.0) * 1e-9

    u10 = C.interp_level(f["Uzm"], p, 10.0, axis=1)
    u60n = C.interp_lat(u10, lat, 60.0, axis=-1)
    u60s = -C.interp_lat(u10, lat, -60.0, axis=-1)      # westerly positive

    ksel = (p >= 10.0) & (p <= 50.0)
    tz = f["THzm"][:, ksel, :] * (p[ksel][None, :, None] / C.P_REF_HPA) ** C.KAPPA
    wgt = np.gradient(np.log(p[ksel])); wgt /= wgt.sum()
    tlev = np.tensordot(tz, wgt, axes=([1], [0]))
    tcapN = C.band_mean(tlev, lat, 60.0, 90.0, axis=-1)
    tcapS = C.band_mean(tlev, lat, -90.0, -60.0, axis=-1)

    vpt100 = C.band_mean(C.interp_level(f["VTHzm"], p, 100.0, axis=1),
                         lat, 45.0, 75.0, axis=-1)

    # annual and seasonal reductions
    out = {}
    yrs, out["mass_flux"] = C.annual_means(mflux, yr)
    _, out["w_star"] = C.annual_means(w70 * 1e3, yr)
    out["_annual_years"] = yrs

    sN, out["vortex_NH"] = C.seasonal_means(u60n, yr, mo, [12, 1, 2], C.season_year)
    _, out["polar_cap_T_NH"] = C.seasonal_means(tcapN, yr, mo, [12, 1, 2], C.season_year)
    _, out["heat_flux_100"] = C.seasonal_means(vpt100, yr, mo, [12, 1, 2], C.season_year)
    out["_djf_years"] = sN

    sS, out["vortex_SH"] = C.seasonal_means(u60s, yr, mo, [6, 7, 8])
    _, out["polar_cap_T_SH"] = C.seasonal_means(tcapS, yr, mo, [6, 7, 8])
    out["_jja_years"] = sS

    # daily DJF vortex, for the daily sigma and the percentiles
    djf = C.djf_mask(mo)
    out["_u60n_djf_daily"] = u60n[djf]

    # shape checks - the per-year series a threshold on the seasonal cycle, the
    # daily distribution and the tropical w* profile needs. Derived here because
    # this is the tape stage; 16, 17 and 18 stay JSON-only and take seconds.
    daily = {"mass_flux": mflux, "w_star": w70 * 1e3, "vortex_NH": u60n,
             "vortex_SH": u60s, "polar_cap_T_NH": tcapN, "polar_cap_T_SH": tcapS}
    my, mrows = monthly_by_year(daily, yr, mo)
    out["_monthly_years"] = my
    for k, v in mrows.items():
        out[f"_monthly_{k}"] = v

    py, prof = profile_by_year(wstar, p, lat, yr)
    out["_profile_years"] = py
    for k, v in prof.items():
        out[f"_profile_{k}"] = v

    wy, pct = winter_percentiles(u60n, yr, mo)
    out["_djf_pctl_years"] = wy
    out["_djf_pctl"] = pct

    # SSW
    ssw = C.detect_ssw_cp07(u60n, yr, mo, dy, "NH")
    nwin = len(np.unique(C.season_year(yr, mo)[np.isin(mo, [11, 12, 1, 2, 3])])) - 1
    out["_ssw_count"] = len(ssw)
    out["_ssw_winters"] = nwin
    out["_ssw_seasons"] = [e["season"] for e in ssw]
    return out


def stats_of(x, detrend=True):
    x = np.asarray(x, dtype=float)
    t = np.arange(len(x), dtype=float)
    sl, ic, r, p, se = stats.linregress(t, x)
    det = x - (ic + sl * t)
    sig_raw = float(x.std(ddof=1))
    sig_det = float(det.std(ddof=1))
    use = sig_det if (detrend and p < 0.05) else sig_raw
    return dict(mean=float(x.mean()), n=len(x), sigma_raw=sig_raw,
                sigma_detrended=sig_det, sigma_used=use,
                trend_per_decade=float(sl * 10), trend_p=float(p),
                trend_significant=bool(p < 0.05))


SCALARS = [
    ("mass_flux", "1e9 kg/s", True),
    ("w_star", "mm/s", True),
    ("vortex_NH", "m/s", False),
    ("vortex_SH", "m/s", False),
    ("polar_cap_T_NH", "K", True),
    ("polar_cap_T_SH", "K", True),
]


def main():
    os.makedirs(C.OUTDIR, exist_ok=True)
    D = {}
    for seg in (TRAIN, TEST):
        print(f"loading {seg} ...", flush=True)
        D[seg] = diagnostics_for(seg)

    res = {"train": TRAIN, "test": TEST, "rows": []}

    # persist every series so the figure script does not recompute them
    res["series"] = {}
    for seg in (TRAIN, TEST):
        res["series"][seg] = {k: (np.asarray(v).tolist() if not np.isscalar(v) else v)
                              for k, v in D[seg].items()}

    print(f"\n{'='*100}")
    print(f"TARGETS SET ON {TRAIN}  ->  VERIFIED ON {TEST}")
    print(f"{'='*100}")
    print(f"{'diagnostic':17s} {'train mean':>11} {'test mean':>11} "
          f"{'diff':>9} {'target':>9} {'sigma':>8}  verdict")
    print("-" * 100)

    for name, unit, detr in SCALARS:
        a = stats_of(D[TRAIN][name], detrend=detr)
        b = stats_of(D[TEST][name], detrend=detr)
        tol = 0.5 * a["sigma_used"]
        diff = b["mean"] - a["mean"]
        ok = abs(diff) <= tol
        # sigma ratio, with the 95% sampling window for these two sample sizes
        rr = np.sqrt(1 / (2 * (b["n"] - 1)) + 1 / (2 * (a["n"] - 1)))
        ratio = b["sigma_used"] / a["sigma_used"]
        ratio_ok = abs(ratio - 1) <= 1.96 * rr
        row = dict(diagnostic=name, units=unit, train=a, test=b,
                   target_0p5sigma=float(tol), difference=float(diff),
                   passes_bias=bool(ok), sigma_ratio=float(ratio),
                   sigma_ratio_window=[float(1 - 1.96 * rr), float(1 + 1.96 * rr)],
                   passes_sigma_ratio=bool(ratio_ok))
        res["rows"].append(row)
        print(f"{name:17s} {a['mean']:11.4f} {b['mean']:11.4f} {diff:+9.4f} "
              f"{tol:9.4f} {a['sigma_used']:8.4f}  "
              f"{'PASS' if ok else 'FAIL'}"
              f"   sigma ratio {ratio:.2f} "
              f"[{1-1.96*rr:.2f},{1+1.96*rr:.2f}] "
              f"{'ok' if ratio_ok else 'OUT'}")

    # ---------------------------------------------------------------- SSW
    ka, na = D[TRAIN]["_ssw_count"], D[TRAIN]["_ssw_winters"]
    kb, nb = D[TEST]["_ssw_count"], D[TEST]["_ssw_winters"]
    lam = ka / na
    lo, hi = 0.5 * lam, 1.5 * lam
    ok = lo <= kb / nb <= hi
    # Poisson prediction interval for the TEST count given the TRAIN rate
    pred_lo = stats.poisson.ppf(0.025, lam * nb)
    pred_hi = stats.poisson.ppf(0.975, lam * nb)
    res["ssw"] = dict(train_count=ka, train_winters=na, train_rate=lam,
                      test_count=kb, test_winters=nb, test_rate=kb / nb,
                      window_pm50=[lo, hi], passes=bool(ok),
                      poisson_pred_interval_counts=[int(pred_lo), int(pred_hi)])
    print(f"\n{'SSW NH':17s} {lam:11.3f} {kb/nb:11.3f} {kb/nb-lam:+9.3f} "
          f"{'+/-50%':>9} {'':8}  {'PASS' if ok else 'FAIL'}"
          f"   ({kb} events in {nb} winters; "
          f"predicted {int(pred_lo)}-{int(pred_hi)} from the train rate)")

    # ---------------------------------------------------------------- daily sigma
    da, db = D[TRAIN]["_u60n_djf_daily"], D[TEST]["_u60n_djf_daily"]
    tau = 14.0
    na_e, nb_e = len(da) / tau, len(db) / tau
    rr = np.sqrt(1 / (2 * na_e) + 1 / (2 * nb_e))
    ratio = db.std(ddof=1) / da.std(ddof=1)
    ok = abs(ratio - 1) <= 1.96 * rr
    res["daily_sigma_DJF"] = dict(
        train=float(da.std(ddof=1)), test=float(db.std(ddof=1)),
        ratio=float(ratio), window=[float(1 - 1.96 * rr), float(1 + 1.96 * rr)],
        passes=bool(ok),
        train_p05=float(np.percentile(da, 5)), train_p95=float(np.percentile(da, 95)),
        test_p05=float(np.percentile(db, 5)), test_p95=float(np.percentile(db, 95)))
    print(f"{'daily sigma DJF':17s} {da.std(ddof=1):11.3f} {db.std(ddof=1):11.3f} "
          f"{db.std(ddof=1)-da.std(ddof=1):+9.3f} {'ratio':>9} {'':8}  "
          f"{'PASS' if ok else 'FAIL'}   ratio {ratio:.3f} "
          f"[{1-1.96*rr:.2f},{1+1.96*rr:.2f}]")
    print(f"{'p05 / p95 DJF':17s} {res['daily_sigma_DJF']['train_p05']:5.1f} /"
          f"{res['daily_sigma_DJF']['train_p95']:5.1f} "
          f"{res['daily_sigma_DJF']['test_p05']:5.1f} /"
          f"{res['daily_sigma_DJF']['test_p95']:5.1f}")

    # ---------------------------------------------------------------- mechanism
    print(f"\n{'mechanism slopes':17s} {'train':>18} {'test':>18}   overlap")
    res["mechanism"] = {}
    for tag, xk, yk in (("R1 wave->vortex", "heat_flux_100", "vortex_NH"),
                        ("R2 thermal wind", "polar_cap_T_NH", "vortex_NH")):
        out = {}
        for seg in (TRAIN, TEST):
            x = np.asarray(D[seg][xk]); y = np.asarray(D[seg][yk])
            n = min(len(x), len(y))
            b, a0 = np.polyfit(x[:n], y[:n], 1)
            rng = np.random.default_rng(20260813)
            sl = np.array([np.polyfit(*(lambda k: (x[:n][k], y[:n][k]))(
                rng.integers(0, n, n)), 1)[0] for _ in range(4000)])
            out[seg] = dict(slope=float(b), ci95=[float(np.percentile(sl, 2.5)),
                                                  float(np.percentile(sl, 97.5))],
                            r=float(np.corrcoef(x[:n], y[:n])[0, 1]), n=int(n))
        A, B = out[TRAIN], out[TEST]
        overlap = not (B["ci95"][1] < A["ci95"][0] or B["ci95"][0] > A["ci95"][1])
        inside = A["ci95"][0] <= B["slope"] <= A["ci95"][1]
        out["test_slope_inside_train_CI"] = bool(inside)
        out["CIs_overlap"] = bool(overlap)
        res["mechanism"][tag] = out
        print(f"{tag:17s} {A['slope']:7.3f} [{A['ci95'][0]:6.3f},{A['ci95'][1]:6.3f}] "
              f"{B['slope']:7.3f} [{B['ci95'][0]:6.3f},{B['ci95'][1]:6.3f}]   "
              f"{'PASS' if inside else 'FAIL'}")

    # ---------------------------------------------------------------- summary
    nb_pass = sum(1 for r in res["rows"] if r["passes_bias"])
    print(f"\n{'='*100}")
    print(f"bias targets: {nb_pass}/{len(res['rows'])} pass out of sample")
    for r in res["rows"]:
        if not r["passes_bias"]:
            print(f"  FAIL {r['diagnostic']}: difference {r['difference']:+.4f} "
                  f"{r['units']} vs target {r['target_0p5sigma']:.4f}; "
                  f"train trend {r['train']['trend_per_decade']:+.4f}/decade "
                  f"(p={r['train']['trend_p']:.4f}), "
                  f"test trend {r['test']['trend_per_decade']:+.4f}/decade "
                  f"(p={r['test']['trend_p']:.4f})")

    with open(OUT, "w") as f:
        json.dump(res, f, indent=2, default=float)
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
