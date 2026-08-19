"""
16 - The operational anchor: CESM 1970-2014, the whole record.

The anchors this replaced were deliberately short. Tier 1 sat on 1996-2014 and
tier 2 on 1980-2014, because a threshold has to be verified on CESM output that
did NOT set it (D12), and that costs half the record. Both are in stale/.

This script answers the other question. Once the method is trusted, the best
estimate of CESM's mean and variability is the one that uses every year available,
and a candidate model that is not CESM does not sit inside the anchor at all.
So: one fixed anchor on the full 45 years, for scoring anything that is not this
CESM run.

What that buys, and what it costs:

  BUYS   the largest sample the run affords - 44 annual years, 42 DJF and 45 JJA
         seasons - so sigma is the best estimated and the variance windows the
         tightest of any anchor in the repo.

  COSTS  the out-of-sample verification. Every CESM year is inside this anchor,
         so nothing is left over to test it with. The archived split-anchor
         scripts in stale/ remain the evidence that the METHOD holds up out of
         sample; this anchor cannot re-earn it, and does not claim to. Section 2
         of the protocol states that trade.

  CARRIES the forced BDC trend. Over the pooled record the mass-flux trend is
         significant at p = 2e-09 and its raw sigma is 54% larger than its
         detrended sigma. Sigma is detrended wherever the trend is significant
         (window_stats, the same p < 0.05 rule as 07), so the band measures
         variability rather than trend. The MEAN still sits mid-trend, which is
         why a period-matched comparison stays the better test for a rollout
         that overlaps CESM in time (D9).

It also carries the tier thresholds, which is what 14 did against the split
anchors: the per-year screening band at tier 1, and the mean tolerance, variance
windows, count intervals and resolvability limits at tier 2. Every number in
EVALUATION_PROTOCOL.md sections 1 and 3 is transcribed from this file.

Reads output/07_period_split.json only - both segments, for their series.
Output: output/16_anchors_45yr.json
"""

import json
import os
import numpy as np
from scipy.stats import norm as _n, poisson as _p

import aide_val_common as C

_norm_sf, _pois_ppf = _n.sf, _p.ppf

OUT = os.path.join(C.OUTDIR, "16_anchors_45yr.json")
ANCHOR = (1970, 2014)
K_SCREEN = 3.0                # tier-1 per-year band, in sigma (protocol section 1)
SCREEN, FULL = 5, 35          # the two tier lengths the protocol scores at
TAU_DAYS = 14.0               # DJF decorrelation time, for the daily sigma (D7)

SERIES = [
    ("mass_flux", "_annual_years", "1e9 kg/s"),
    ("w_star", "_annual_years", "mm/s"),
    ("vortex_NH", "_djf_years", "m/s"),
    ("vortex_SH", "_jja_years", "m/s"),
    ("polar_cap_T_NH", "_djf_years", "K"),
    ("polar_cap_T_SH", "_jja_years", "K"),
]

MECHANISM = [("R1 wave->vortex", "heat_flux_100", "vortex_NH", "_djf_years"),
             ("R2 thermal wind", "polar_cap_T_NH", "vortex_NH", "_djf_years")]


def slope_ci(x, y, seed=20260813, n_boot=4000):
    n = min(len(x), len(y))
    x, y = x[:n], y[:n]
    b = float(np.polyfit(x, y, 1)[0])
    rng = np.random.default_rng(seed)
    draws = np.array([np.polyfit(*(lambda k: (x[k], y[k]))(rng.integers(0, n, n)), 1)[0]
                      for _ in range(n_boot)])
    return dict(slope=b, ci95=[float(np.percentile(draws, 2.5)),
                               float(np.percentile(draws, 97.5))],
                r=float(np.corrcoef(x, y)[0, 1]), n=int(n))


def main():
    ps = json.load(open(os.path.join(C.OUTDIR, "07_period_split.json")))
    S = ps["series"]
    segs = [ps["test"], ps["train"]]                   # 1970-1995, then 1996-2014
    lo, hi = ANCHOR

    res = {"anchor_period": f"{lo}-{hi}", "anchor_years": list(ANCHOR),
           "segments_joined": segs, "k_screen": K_SCREEN,
           "sigma_rule": "detrended within the window where the trend has p < 0.05",
           "diagnostics": {}}

    print(f"OPERATIONAL ANCHOR - CESM {lo}-{hi}, the whole record")
    print(f"{'diagnostic':16s} {'n':>3} {'mean':>10} {'sig_raw':>8} {'sig_det':>8} "
          f"{'trend/dec':>10} {'p':>9} {'det':>4} "
          f"{f'+/-{K_SCREEN:.0f} sigma band':>24} {'worst yr':>9}")
    print("-" * 118)

    for key, ykey, unit in SERIES:
        y, v = C.join_segments(S, segs, key, ykey)
        A = C.window_stats(y, v, lo, hi)
        sig, mu = A["sigma_used"], A["mean"]
        band = [mu - K_SCREEN * sig, mu + K_SCREEN * sig]
        m = (y >= lo) & (y <= hi)
        z = (v[m] - mu) / sig
        worst = int(y[m][np.argmax(np.abs(z))])
        A.update(units=unit, screening_band=band,
                 k_screen=K_SCREEN,
                 years_outside_band=[int(t) for t in y[m][np.abs(z) > K_SCREEN]],
                 widest_excursion_sigma=float(np.abs(z).max()),
                 widest_excursion_year=worst,
                 tolerance_0p5sigma=float(0.5 * sig))
        res["diagnostics"][key] = A
        print(f"{key:16s} {A['n']:3d} {mu:10.4f} {A['sigma_raw']:8.4f} "
              f"{A['sigma_detrended']:8.4f} {A['trend_per_decade']:+10.4f} "
              f"{A['trend_p']:9.2e} {'yes' if A['detrended'] else 'no':>4} "
              f"{band[0]:11.4f} to {band[1]:8.4f} "
              f"{np.abs(z).max():+6.2f}sig")

    n_out = sum(len(d["years_outside_band"]) for d in res["diagnostics"].values())
    n_chk = sum(d["n"] for d in res["diagnostics"].values())
    res["band_self_consistency"] = dict(
        checks=n_chk, outside=n_out,
        note="anchor years scored against their own band; a non-zero count is the "
             "forced trend, not scatter")
    print(f"\n  the anchor's own years against its own band: {n_out} of {n_chk} outside "
          f"+/-{K_SCREEN:.0f} sigma")

    # ----------------------------------------------------------------- daily DJF
    daily = np.array(sum((S[s]["_u60n_djf_daily"] for s in segs), []), dtype=float)
    res["daily_DJF_u60N"] = dict(
        sigma=float(daily.std(ddof=1)), n_days=int(len(daily)),
        effective_n=float(len(daily) / TAU_DAYS), tau_days=TAU_DAYS,
        p05=float(np.percentile(daily, 5)), p95=float(np.percentile(daily, 95)))
    d = res["daily_DJF_u60N"]
    print(f"  daily DJF u 60N: sigma {d['sigma']:.3f} m/s over {d['n_days']} days "
          f"({d['effective_n']:.0f} effective samples), p05/p95 "
          f"{d['p05']:.1f}/{d['p95']:.1f}")

    # ----------------------------------------------------------------- SSW rate
    ssw = np.array(sum((S[s]["_ssw_seasons"] for s in segs), []), dtype=float)
    winters = sum(S[s]["_ssw_winters"] for s in segs)
    count = int(((ssw >= lo) & (ssw <= hi)).sum())
    res["ssw_NH"] = dict(count=count, winters=int(winters),
                         rate_per_winter=float(count / winters))
    print(f"  major NH SSW: {count} in {winters} winters, "
          f"{count / winters:.3f}/winter")

    # ----------------------------------------------------------------- mechanism
    res["mechanism"] = {}
    for tag, xk, yk, ykey in MECHANISM:
        yx, vx = C.join_segments(S, segs, xk, ykey)
        yy, vy = C.join_segments(S, segs, yk, ykey)
        m = (yx >= lo) & (yx <= hi)
        res["mechanism"][tag] = slope_ci(vx[m], vy[m])
        M = res["mechanism"][tag]
        print(f"  {tag:17s} slope {M['slope']:+.3f} "
              f"[{M['ci95'][0]:+.3f}, {M['ci95'][1]:+.3f}]  r = {M['r']:+.2f}  "
              f"n = {M['n']}")

    # ------------------------------------------------------------------- tiers
    # The thresholds themselves, at the two lengths the protocol scores at. Tier 1
    # gates individual years; only its mean row depends on the 5-year length.
    res["tiers"] = {"screen_length": SCREEN, "full_length": FULL, "mean": {},
                    "variance": {}, "trend": {}, "mechanism": {}}
    print(f"\nTIER THRESHOLDS from this anchor - tier 1 at n = {SCREEN}, "
          f"tier 2 at n = {FULL}")
    print(f"{'diagnostic':16s} {'sigma':>8} {'t1 mean tol':>12} {'t2 mean tol':>12} "
          f"{'t2 tol %mean':>13} {'branch at t2':>13}")
    for key, ykey, unit in SERIES:
        d = res["diagnostics"][key]
        sig, mu = d["sigma_used"], d["mean"]
        t1, t2 = C.bias_target(sig, SCREEN), C.bias_target(sig, FULL)
        res["tiers"]["mean"][key] = dict(
            units=unit, sigma=sig, anchor_mean=mu,
            tier1_advisory_tolerance=float(t1),
            tier1_advisory_in_sigma=float(t1 / sig),
            tier2_tolerance=float(t2), tier2_in_sigma=float(t2 / sig),
            tier2_pct_of_mean=float(100 * t2 / abs(mu)),
            binding_branch=("0.5 sigma" if 0.5 * sig >= 1.96 * sig / np.sqrt(FULL)
                            else "detection"))
        m = res["tiers"]["mean"][key]
        print(f"{key:16s} {sig:8.4f} {t1:12.4f} {t2:12.4f} "
              f"{m['tier2_pct_of_mean']:12.1f}% {m['binding_branch']:>13}")
        res["tiers"]["variance"][key] = dict(
            interannual_ratio_window=list(C.ratio_window(FULL, d["n"])),
            anchor_n=d["n"])

    dd = res["daily_DJF_u60N"]
    per_year = dd["effective_n"] / res["diagnostics"]["vortex_NH"]["n"]
    res["tiers"]["variance"]["daily_DJF_u60N"] = dict(
        ratio_window=list(C.ratio_window(FULL, res["diagnostics"]["vortex_NH"]["n"],
                                         per_year)),
        effective_samples_per_winter=float(per_year))
    iw = res["tiers"]["variance"]["mass_flux"]["interannual_ratio_window"]
    dw = res["tiers"]["variance"]["daily_DJF_u60N"]["ratio_window"]
    print(f"\n  sigma ratio windows for a {FULL}-yr candidate: "
          f"interannual {iw[0]:.2f}-{iw[1]:.2f}, daily DJF {dw[0]:.2f}-{dw[1]:.2f} "
          f"({per_year:.1f} effective samples per winter)")

    # False alarms: 6 diagnostics x 5 years is 30 simultaneous checks per screen.
    n_chk = len(SERIES) * SCREEN
    res["tiers"]["false_alarm"] = {"n_checks_per_screen": n_chk, "by_band": {}}
    for k in (2.0, 2.5, 3.0, 3.5):
        per = float(2 * _norm_sf(k))
        res["tiers"]["false_alarm"]["by_band"][f"{k:.1f}"] = dict(
            per_check=per, any_flag=float(1 - (1 - per) ** n_chk))
    print(f"  a perfect Gaussian model trips the +/-{K_SCREEN:.0f} sigma gate on "
          f"{100 * res['tiers']['false_alarm']['by_band']['3.0']['any_flag']:.1f}% "
          f"of screening runs ({n_chk} checks each)")

    lam = res["ssw_NH"]["rate_per_winter"]
    res["tiers"]["ssw_count"] = {}
    for tag, n in (("tier1", SCREEN), ("tier2", FULL)):
        res["tiers"]["ssw_count"][tag] = dict(
            winters=n, expected=float(lam * n),
            interval=[int(_pois_ppf(0.025, lam * n)),
                      int(_pois_ppf(0.975, lam * n))])
        c = res["tiers"]["ssw_count"][tag]
        print(f"  major NH SSW over {n:2d} winters: expect {c['expected']:5.1f}, "
              f"accept {c['interval'][0]}-{c['interval'][1]}")

    for key, ykey, unit in SERIES:
        d = res["diagnostics"][key]
        se = d["trend_se_per_decade"] * (d["n"] / FULL) ** 1.5
        res["tiers"]["trend"][key] = dict(
            units=unit + " per decade", trend=d["trend_per_decade"],
            se_at_full=float(se),
            resolvable_at_full=bool(abs(d["trend_per_decade"]) > 1.96 * se))
        t = res["tiers"]["trend"][key]
        print(f"  trend {key:16s} {d['trend_per_decade']:+.4f} +/- {1.96 * se:.4f} "
              f"per decade at {FULL} yr  "
              + ("resolvable" if t["resolvable_at_full"] else "NOT resolvable"))

    for tag, M in res["mechanism"].items():
        half = 0.5 * (M["ci95"][1] - M["ci95"][0])
        h = half * np.sqrt(M["n"] / FULL)
        res["tiers"]["mechanism"][tag] = dict(
            slope=M["slope"], half_width_full=float(h),
            resolvable_at_full=bool(abs(M["slope"]) > h))
        print(f"  {tag:17s} slope {M['slope']:+.3f} +/-{h:.2f} at {FULL} yr  "
              + ("resolvable" if abs(M["slope"]) > h else "NOT resolvable"))

    with open(OUT, "w") as f:
        json.dump(res, f, indent=2)
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
