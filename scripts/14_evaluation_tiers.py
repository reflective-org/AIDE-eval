"""
14 - The two-tier evaluation protocol: 5-year screening, 35-year validation.

Script 08 answers "how does a target scale with rollout length?" for one focus
length, n = 10. This script answers a different question: what are the numbers for
the two lengths the evaluation protocol actually uses?

  TIER 1   n = 5    every new model version. A regression check, not a validation.
                    Individual years are scored against a LOOSE envelope, taken
                    empirically from the range CESM's own years occupy. The point
                    is to catch a broken configuration, not to certify a good one.

  TIER 2   n = 35   once a configuration looks promising. Mean AND variance, plus
                    the metrics that only become resolvable at that length.

Two things are deliberately different between the tiers.

  * Tier 1 scores INDIVIDUAL YEARS; tier 2 scores the ROLLOUT MEAN. A single year
    carries the full interannual sigma, so a per-year band has to be several sigma
    wide or a perfect model fails it about a third of the time (see D5, and the
    "most common misreading" section of the technical report).

  * Tier 1 thresholds are empirical and provisional. They are the observed spread
    of CESM, widened by nothing, and they are expected to TIGHTEN as real failure
    modes are found. Tier 2 thresholds are the derived ones and are not adjustable
    by inspection.

Reads only output/07_period_split.json. Nothing here feeds 05 or the PDFs, so the
existing artefacts are untouched.

Output: output/14_evaluation_tiers.json
"""

import json
import os
import numpy as np
from scipy import stats

import aide_val_common as C

OUT = os.path.join(C.OUTDIR, "14_evaluation_tiers.json")
SCREEN, FULL = 5, 35
SAMPLES_PER_WINTER = 6.4      # daily DJF, decorrelation time 14 d (D7)

# Tier 1 stays on the 1996-2014 anchor that sets every sigma elsewhere in the repo.
# Tier 2 is anchored on a full 35 years, the length of the rollout it scores, and is
# tested against 1970-1994. Two caveats, both stated in the protocol document:
#   * the windows OVERLAP by 15 years (1980-1994), so the test is not independent the
#     way 07's 1996-2014 / 1970-1995 split is;
#   * 1980-2014 spans the 1995/96 restart, so the anchor mixes two separate CESM runs.
T2_ANCHOR = (1980, 2014)
T2_TEST = (1970, 1994)
K_SCREEN = 3.0                # per-year screening band, in sigma. See the printout:
                              # no CESM anchor year of any diagnostic reaches 2.6 sigma,
                              # so this is loose by construction and is the number to
                              # tighten first once real failure modes appear.

SERIES = [
    ("mass_flux", "_annual_years", "1e9 kg/s"),
    ("w_star", "_annual_years", "mm/s"),
    ("vortex_NH", "_djf_years", "m/s"),
    ("vortex_SH", "_jja_years", "m/s"),
    ("polar_cap_T_NH", "_djf_years", "K"),
    ("polar_cap_T_SH", "_jja_years", "K"),
]


def bias_target(sigma, n):
    """D5, EVALUATION_PROTOCOL.md appendix A."""
    return max(0.5 * sigma, 1.96 * sigma / np.sqrt(n))


def ratio_window(n, n_ref, per_year=1.0):
    """95% window on a sigma ratio between an n-year rollout and the n_ref anchor."""
    ne, nr = n * per_year, n_ref * per_year
    rr = np.sqrt(1 / (2 * (ne - 1)) + 1 / (2 * (nr - 1)))
    return float(1 - 1.96 * rr), float(1 + 1.96 * rr)


def pooled(S, segs, key, ykey):
    """The 1970-2014 record for one diagnostic, both segments joined and sorted."""
    y = np.array(sum((S[s][ykey] for s in segs), []), dtype=float)
    v = np.array(sum((S[s][key] for s in segs), []), dtype=float)
    o = np.argsort(y)
    return y[o], v[o]


def window_stats(y, v, lo, hi):
    """Mean and sigma over [lo, hi], detrended when the trend is significant.

    Same rule as 07_period_split.stats_of, so the two are directly comparable.
    """
    m = (y >= lo) & (y <= hi)
    x = v[m]
    t = np.arange(len(x), dtype=float)
    sl, ic, r, p, se = stats.linregress(t, x)
    det = x - (ic + sl * t)
    sig_raw, sig_det = float(x.std(ddof=1)), float(det.std(ddof=1))
    return dict(n=int(m.sum()), mean=float(x.mean()), sigma_raw=sig_raw,
                sigma_detrended=sig_det,
                sigma_used=(sig_det if p < 0.05 else sig_raw),
                trend_p=float(p), detrended=bool(p < 0.05))


def main():
    ps = json.load(open(os.path.join(C.OUTDIR, "07_period_split.json")))
    train, test = ps["train"], ps["test"]
    rows = {r["diagnostic"]: r for r in ps["rows"]}
    S = ps["series"]

    res = {"screen_length": SCREEN, "full_length": FULL, "anchor": train,
           "note": "tier 1 = per-year empirical envelope; tier 2 = mean and variance"}

    # ------------------------------------------------------- tier 1, n = 5
    print(f"TIER 1 - {SCREEN}-year screening, INDIVIDUAL YEARS vs an empirical envelope")
    print(f"{'diagnostic':16s} {'anchor mean':>11} {'sigma':>8} "
          f"{'anchor yrs, sigma':>18} {'pooled yrs, sigma':>18} "
          f"{f'+/-{K_SCREEN:.0f} sigma band':>26}")
    print("-" * 104)

    res["tier1"] = {}
    for key, ykey, unit in SERIES:
        sig = rows[key]["train"]["sigma_used"]
        mu = rows[key]["train"]["mean"]
        pool_v = np.array(S[test][key] + S[train][key], float)
        anchor = np.array(S[train][key], float)
        env = (float(pool_v.min()), float(pool_v.max()))
        env_a = (float(anchor.min()), float(anchor.max()))
        res["tier1"][key] = dict(
            units=unit, anchor_mean=mu, sigma=sig, n_years_pooled=len(pool_v),
            envelope_pooled=list(env), envelope_anchor_only=list(env_a),
            envelope_pooled_in_sigma=[(env[0] - mu) / sig, (env[1] - mu) / sig],
            envelope_anchor_in_sigma=[(env_a[0] - mu) / sig, (env_a[1] - mu) / sig],
            widest_excursion_sigma=float(max(abs(env[0] - mu), abs(env[1] - mu)) / sig),
            widest_anchor_excursion_sigma=float(
                max(abs(env_a[0] - mu), abs(env_a[1] - mu)) / sig),
            screening_band=[float(mu - K_SCREEN * sig), float(mu + K_SCREEN * sig)],
            mean_tolerance=float(bias_target(sig, SCREEN)),
            mean_tolerance_in_sigma=float(bias_target(sig, SCREEN) / sig))
        b = res["tier1"][key]["screening_band"]
        print(f"{key:16s} {mu:11.4f} {sig:8.4f} "
              f"{(env_a[0]-mu)/sig:+8.2f} {(env_a[1]-mu)/sig:+8.2f} "
              f"{(env[0]-mu)/sig:+9.2f} {(env[1]-mu)/sig:+8.2f} "
              f"{b[0]:12.4f} to {b[1]:10.4f}")

    worst_a = max(v["widest_anchor_excursion_sigma"] for v in res["tier1"].values())
    worst_p = max(v["widest_excursion_sigma"] for v in res["tier1"].values())
    res["tier1_summary"] = dict(
        k_screen=K_SCREEN, widest_anchor_excursion_sigma=float(worst_a),
        widest_pooled_excursion_sigma=float(worst_p),
        mean_tolerance_in_sigma=float(bias_target(1.0, SCREEN)))
    print(f"\n  widest excursion of any single year: {worst_a:.2f} sigma over the "
          f"{train} anchor, {worst_p:.2f} sigma over the pooled 1970-2014 record")
    print(f"  the pooled figure is the BDC trend, not natural scatter - it is "
          f"mass_flux and w_star that blow out (D9). Score a rollout against the "
          f"CESM years it covers.")
    print(f"  +/-{K_SCREEN:.0f} sigma therefore clears every anchor year with margin; "
          f"tighten this number, not the derived ones.")
    print(f"  a {SCREEN}-year MEAN is held to "
          f"{bias_target(1.0, SCREEN):.3f} sigma (advisory at this tier)")

    # How often would a PERFECT model trip the per-year gate? 6 diagnostics x 5 years
    # is 30 simultaneous checks, so the multiplicity matters more than intuition says.
    nchk = len(SERIES) * SCREEN
    res["tier1_summary"]["n_checks_per_screen"] = nchk
    res["tier1_summary"]["false_alarm_rate"] = {}
    print(f"\n  false-alarm rate of the per-year gate, {nchk} checks per screening run,")
    print(f"  for a model that is PERFECT and merely Gaussian:")
    for k in (2.0, 2.5, 3.0, 3.5):
        per = 2 * stats.norm.sf(k)
        anyflag = 1 - (1 - per) ** nchk
        res["tier1_summary"]["false_alarm_rate"][f"{k:.1f}"] = dict(
            per_check=float(per), any_flag=float(anyflag))
        print(f"    +/-{k:.1f} sigma   {100*per:6.2f}% per check   "
              f"{100*anyflag:5.1f}% chance of at least one flag")

    lam = ps["ssw"]["train_rate"]
    res["tier1_summary"]["ssw_count"] = dict(
        expected=float(lam * SCREEN),
        interval=[int(stats.poisson.ppf(0.025, lam * SCREEN)),
                  int(stats.poisson.ppf(0.975, lam * SCREEN))])
    s1 = res["tier1_summary"]["ssw_count"]
    print(f"\n  major SSW over {SCREEN} winters: expect {s1['expected']:.1f}, "
          f"95% interval [{s1['interval'][0]}, {s1['interval'][1]}] - wide enough that "
          f"only a dead or a hyperactive vortex trips it, which is what tier 1 is for")

    # ------------------------------------------------------ tier 2, n = 35
    a_lo, a_hi = T2_ANCHOR
    t_lo, t_hi = T2_TEST
    segs = [test, train]
    print(f"\nTIER 2 - {FULL}-year validation, anchored on {a_lo}-{a_hi}, "
          f"tested on {t_lo}-{t_hi}")
    print(f"  NOTE the windows overlap by {t_hi - a_lo + 1} years "
          f"({a_lo}-{t_hi}); the test is not independent.")
    print(f"  NOTE {a_lo}-{a_hi} spans the 1995/96 restart between two separate runs.")
    print(f"\n{'diagnostic':16s} {'anchor mean':>11} {'sigma':>8} {'mean tol':>9} "
          f"{'test mean':>10} {'offset':>8} {'verdict':>8} {'sig ratio':>10} {'':>7}")
    print("-" * 100)

    res["tier2"] = {"anchor_period": f"{a_lo}-{a_hi}", "test_period": f"{t_lo}-{t_hi}",
                    "overlap_years": [a_lo, t_hi], "mean": {}, "test": {}}
    for key, ykey, unit in SERIES:
        y, v = pooled(S, segs, key, ykey)
        A = window_stats(y, v, a_lo, a_hi)
        T = window_stats(y, v, t_lo, t_hi)
        sig, mu = A["sigma_used"], A["mean"]
        t = bias_target(sig, FULL)
        branch = "0.5 sigma" if 0.5 * sig >= 1.96 * sig / np.sqrt(FULL) else "detection"
        res["tier2"]["mean"][key] = dict(
            units=unit, sigma=sig, anchor_mean=mu, anchor_n=A["n"],
            anchor_detrended=A["detrended"], tolerance=float(t),
            tolerance_in_sigma=float(t / sig),
            tolerance_pct_of_mean=float(100 * t / abs(mu)), binding_branch=branch)

        rlo, rhi = ratio_window(T["n"], A["n"])
        ratio = T["sigma_used"] / sig
        res["tier2"]["test"][key] = dict(
            test_mean=T["mean"], test_n=T["n"], test_sigma=T["sigma_used"],
            test_detrended=T["detrended"],
            difference=float(T["mean"] - mu),
            difference_in_sigma=float((T["mean"] - mu) / sig),
            passes_mean=bool(abs(T["mean"] - mu) <= t),
            sigma_ratio=float(ratio), sigma_ratio_window=[rlo, rhi],
            passes_variance=bool(rlo <= ratio <= rhi))
        R = res["tier2"]["test"][key]
        print(f"{key:16s} {mu:11.4f} {sig:8.4f} {t:9.4f} {T['mean']:10.4f} "
              f"{R['difference_in_sigma']:+8.2f}σ {'PASS' if R['passes_mean'] else 'FAIL':>8} "
              f"{ratio:10.2f} {'pass' if R['passes_variance'] else 'FAIL':>7}")

    n_a = res["tier2"]["mean"]["mass_flux"]["anchor_n"]
    n_t = res["tier2"]["test"]["mass_flux"]["test_n"]
    lo, hi = ratio_window(FULL, n_a)
    dlo, dhi = ratio_window(FULL, n_a, SAMPLES_PER_WINTER)
    res["tier2"]["variance"] = {
        "interannual_sigma_ratio_window": [lo, hi],
        "daily_sigma_ratio_window_DJF": [dlo, dhi],
        "reference_n": n_a,
        "daily_samples_per_winter": SAMPLES_PER_WINTER}
    print(f"\n  a {FULL}-yr rollout vs the {n_a}-yr anchor: interannual sigma ratio "
          f"{lo:.2f} - {hi:.2f}  (width {hi-lo:.2f})")
    print(f"  daily DJF sigma ratio:                                       "
          f"{dlo:.2f} - {dhi:.2f}  (width {dhi-dlo:.2f})")

    # SSW rate re-counted on the tier-2 anchor window.
    ssw_years = np.array(sum((S[s]["_ssw_seasons"] for s in segs), []), dtype=float)
    n_ssw_a = int(((ssw_years >= a_lo) & (ssw_years <= a_hi)).sum())
    n_ssw_t = int(((ssw_years >= t_lo) & (ssw_years <= t_hi)).sum())
    lam = n_ssw_a / (a_hi - a_lo + 1)
    exp_t = lam * (t_hi - t_lo + 1)
    res["tier2"]["ssw_count"] = dict(
        rate_per_winter=float(lam), anchor_count=n_ssw_a, expected=float(lam * FULL),
        interval=[int(stats.poisson.ppf(0.025, lam * FULL)),
                  int(stats.poisson.ppf(0.975, lam * FULL))],
        resolvable_pct=float(100 * 1.96 / np.sqrt(lam * FULL)),
        test_count=n_ssw_t, test_expected=float(exp_t),
        test_interval=[int(stats.poisson.ppf(0.025, exp_t)),
                       int(stats.poisson.ppf(0.975, exp_t))])
    s = res["tier2"]["ssw_count"]
    s["test_passes"] = bool(s["test_interval"][0] <= n_ssw_t <= s["test_interval"][1])
    print(f"  major SSW, rate {lam:.2f}/winter over {a_lo}-{a_hi}: a {FULL}-winter "
          f"rollout should give {s['expected']:.1f}, 95% [{s['interval'][0]}, "
          f"{s['interval'][1]}]")
    print(f"  {t_lo}-{t_hi} has {n_ssw_t} against an expected {exp_t:.1f} "
          f"[{s['test_interval'][0]}, {s['test_interval'][1]}] - "
          f"{'PASS' if s['test_passes'] else 'FAIL'}")

    # The forced BDC trend. An OLS slope SE scales as n^-3/2, so rescale the
    # 44-year fit in 02b to the rollout length.
    tr = json.load(open(os.path.join(C.OUTDIR, "02b_trends.json")))
    res["tier2"]["trend"] = {}
    print()
    for tag, d in tr.items():
        if not isinstance(d, dict) or "trend_se_per_decade" not in d:
            continue
        se = d["trend_se_per_decade"] * (d["n"] / FULL) ** 1.5
        res["tier2"]["trend"][tag] = dict(
            units=d["units"] + " per decade", trend=d["trend_per_decade"],
            se_at_full=float(se), t_at_full=float(abs(d["trend_per_decade"]) / se),
            resolvable_at_full=bool(abs(d["trend_per_decade"]) > 1.96 * se))
        r = res["tier2"]["trend"][tag]
        print(f"  trend {tag:24s} {d['trend_per_decade']:+.4f} +/- {1.96*se:.4f} "
              f"per decade at {FULL} yr  "
              + ("resolvable" if r["resolvable_at_full"] else "NOT resolvable"))

    res["tier2"]["mechanism"] = {}
    for tag, m in ps["mechanism"].items():
        A = m[train]
        half = 0.5 * (A["ci95"][1] - A["ci95"][0])
        h5, h35 = half * np.sqrt(A["n"] / SCREEN), half * np.sqrt(A["n"] / FULL)
        res["tier2"]["mechanism"][tag] = dict(
            anchor_slope=A["slope"], anchor_n=A["n"],
            half_width_screen=float(h5), half_width_full=float(h35),
            resolvable_at_full=bool(abs(A["slope"]) > h35))
        print(f"  {tag:17s} slope {A['slope']:+.3f}  "
              f"+/-{h35:.2f} at {FULL} yr  "
              + ("resolvable" if abs(A["slope"]) > h35 else "NOT resolvable"))

    with open(OUT, "w") as f:
        json.dump(res, f, indent=2)
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
