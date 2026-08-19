"""
16 - The operational anchor: CESM 1970-2014, the whole record.

The protocol's own anchors are deliberately short. Tier 1 sits on 1996-2014 and
tier 2 on 1980-2014, because a threshold has to be verified on CESM output that
did NOT set it (D12), and that costs half the record.

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
         so nothing is left over to test it with. Sections 2 and 3.3 of the
         protocol remain the evidence that the METHOD holds up out of sample;
         this anchor cannot re-earn it, and does not claim to.

  CARRIES the forced BDC trend. Over the pooled record the mass-flux trend is
         significant at p = 2e-09 and its raw sigma is 54% larger than its
         detrended sigma. Sigma is detrended wherever the trend is significant
         (window_stats, the same p < 0.05 rule as 07), so the band measures
         variability rather than trend. The MEAN still sits mid-trend, which is
         why a period-matched comparison stays the better test for a rollout
         that overlaps CESM in time (D9).

Reads output/07_period_split.json only - both segments, for their series.
Output: output/16_anchors_45yr.json
"""

import json
import os
import numpy as np

import aide_val_common as C

OUT = os.path.join(C.OUTDIR, "16_anchors_45yr.json")
ANCHOR = (1970, 2014)
K_SCREEN = 3.0                # tier-1 per-year band, in sigma (protocol section 1)
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

    with open(OUT, "w") as f:
        json.dump(res, f, indent=2)
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
