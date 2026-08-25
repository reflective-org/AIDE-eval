"""
17 - Score a climate model against the 45-year operational anchor.

Usage:  17_validate.py [--climate-model NAME] [first_year last_year]

  --climate-model   name of the model being scored. Default is derived from the
                    year window, so a bare run still names itself.
  first_year last_year   default 1996 2014

This is the scoring step, kept separate from the anchor (16) on purpose: the
anchor is fixed, the climate model changes. Everything the climate-model side
needs arrives through climate_model_series() - one function, one place to swap
when the model is a rollout rather than a window of CESM's own record.

Running it on 1996-2014 exercises the pipeline end to end. That window is INSIDE
the anchor, so the result is a self-consistency check, not a validation: the
verdicts show the machinery works and the thresholds are the right size, and any
pass is weaker evidence than the corresponding failure. Real use is a climate
model the anchor has never seen.

Every artefact is stamped `__<climate model>__<production date>` so a result can
never be read without knowing what was scored and when.

Reads output/16_anchors_45yr.json and output/07_period_split.json.
Output: output/17_validation__<stamp>.json
        validation_results/validation_result__<stamp>.md
"""

import argparse
import datetime
import json
import os
import re
import numpy as np
from scipy import stats

import aide_val_common as C

RESULTS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "validation_results")


def stamp_for(climate_model, produced):
    """`<climate model>__<date>`, reduced to characters that are safe in a path."""
    slug = re.sub(r"[^A-Za-z0-9.+-]+", "-", climate_model).strip("-")
    return f"{slug}__{produced}"


def stamped(name, stamp):
    """Insert the stamp before the extension: a.png -> a__<stamp>.png"""
    root, ext = os.path.splitext(name)
    return f"{root}__{stamp}{ext}"

SERIES = [
    ("mass_flux", "_annual_years", "Mass flux, 70 hPa", "10^9 kg/s", 4),
    ("w_star", "_annual_years", "w*, 10S-10N", "mm/s", 4),
    ("vortex_NH", "_djf_years", "Vortex NH, DJF", "m/s", 2),
    ("vortex_SH", "_jja_years", "Vortex SH, JJA", "m/s", 2),
    ("polar_cap_T_NH", "_djf_years", "Polar cap T, NH", "K", 2),
    ("polar_cap_T_SH", "_jja_years", "Polar cap T, SH", "K", 2),
]
MECHANISM = [("R1 wave->vortex", "heat_flux_100", "vortex_NH", "_djf_years"),
             ("R2 thermal wind", "polar_cap_T_NH", "vortex_NH", "_djf_years")]
BLOCK = 5                     # tier-1 screening runs are 5 years long

# Written by 18_validation_figures.py into the same directory as the report.
# Base names; 18 stamps them the same way this script stamps the report.
FIGURES = [
    ("tier1_screening.png",
     "tier 1 — every climate-model year against the ±3σ band, six diagnostics"),
    ("tier2_mean.png",
     "tier 2 mean — offset from the anchor against the ±0.5σ tolerance"),
    ("tier2_variance.png",
     "tier 2 variance — σ ratios against their 95% windows, interannual and daily"),
    ("counts_and_relations.png",
     "SSW count against its Poisson interval; R1 and R2 slopes against the anchor fit"),
    ("shape_seasonal.png",
     "seasonal cycle — the 12-month climatology, and the annual-harmonic "
     "amplitude and phase against their tolerances"),
    ("shape_daily_distribution.png",
     "daily distribution — per-winter percentiles of u at 60°N against the anchor"),
    ("shape_w_star_profile.png",
     "tropical w* profile — absolute (advisory) and normalised (gated)"),
]


def climate_model_series(ps, lo, hi):
    """The climate model's scalar series, one entry per diagnostic key.

    Here the climate model is a window of CESM's own record, so the series come
    from 07's JSON. A model rollout would replace this function: return the same dict
    of {key: (years, values)} plus the daily DJF u at 60N, computed through
    aide_val_common.tem_residual on the model's own grid (protocol section 5,
    rules 1 and 2).
    """
    S = ps["series"]
    segs = [ps["test"], ps["train"]]
    out = {}
    for key, ykey, _, _, _ in SERIES:
        y, v = C.join_segments(S, segs, key, ykey)
        m = (y >= lo) & (y <= hi)
        out[key] = (y[m], v[m])
    for _, xk, yk, ykey in MECHANISM:
        for k in (xk, yk):
            if k not in out:
                y, v = C.join_segments(S, segs, k, ykey)
                m = (y >= lo) & (y <= hi)
                out[k] = (y[m], v[m])

    # Shape checks. All three are year-labelled, so unlike the daily series below
    # they subset cleanly to any window. A rollout returns the same keys: the
    # per-year 12-month climatology of each diagnostic, the per-year tropical w*
    # at each profile level, and the per-winter DJF percentiles of u at 60N.
    for key, _, _, _, _ in SERIES:
        y, v = C.join_segments(S, segs, f"_monthly_{key}", "_monthly_years")
        m = (y >= lo) & (y <= hi)
        out[f"_monthly_{key}"] = (y[m], np.asarray(v, float)[m])
    for pl in C.PROFILE_LEVELS:
        y, v = C.join_segments(S, segs, f"_profile_{pl:g}", "_profile_years")
        m = (y >= lo) & (y <= hi)
        out[f"_profile_{pl:g}"] = (y[m], np.asarray(v, float)[m])
    y, v = C.join_segments(S, segs, "_djf_pctl", "_djf_pctl_years")
    m = (y >= lo) & (y <= hi)
    out["_djf_pctl"] = (y[m], np.asarray(v, float)[m])

    # Daily DJF is stored per segment with no year labels, so it can only be
    # subset when the window lines up with whole segments.
    daily, aligned = [], True
    for s in segs:
        yy = np.asarray(S[s]["_annual_years"], float)
        s_lo, s_hi = yy.min(), yy.max()
        if s_hi < lo or s_lo > hi:
            continue
        if s_lo >= lo and s_hi <= hi:
            daily += S[s]["_u60n_djf_daily"]
        else:
            aligned = False
    ssw = np.array(sum((S[s]["_ssw_seasons"] for s in segs), []), float)
    winters = sum(S[s]["_ssw_winters"] for s in segs
                  if np.asarray(S[s]["_annual_years"], float).min() >= lo
                  and np.asarray(S[s]["_annual_years"], float).max() <= hi)
    return out, (np.array(daily, float) if aligned else None), ssw, winters


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1].strip())
    ap.add_argument("--climate-model", dest="climate_model", default=None,
                    help="name of the model being scored; default from the window")
    ap.add_argument("years", nargs="*", type=int, default=[],
                    help="first_year last_year (default 1996 2014)")
    args = ap.parse_args()
    lo, hi = (args.years[0], args.years[1]) if len(args.years) > 1 else (1996, 2014)

    climate_model = args.climate_model or f"CESM2.1.5-WACCM6 histSST {lo}-{hi}"
    produced = datetime.date.today().isoformat()
    stamp = stamp_for(climate_model, produced)
    out_j = os.path.join(C.OUTDIR, stamped("17_validation.json", stamp))
    out_m = os.path.join(RESULTS, stamped("validation_result.md", stamp))

    A = json.load(open(os.path.join(C.OUTDIR, "16_anchors_45yr.json")))
    ps = json.load(open(os.path.join(C.OUTDIR, "07_period_split.json")))
    cm, daily, ssw_years, winters = climate_model_series(ps, lo, hi)

    res = {"climate_model": climate_model, "produced": produced, "stamp": stamp,
           "climate_model_period": [lo, hi], "anchor_period": A["anchor_period"],
           "estimator": "aide_val_common.tem_residual",
           "inside_anchor": bool(lo >= A["anchor_years"][0] and hi <= A["anchor_years"][1]),
           "tier1": {}, "tier2_mean": {}, "tier2_variance": {}, "counts": {},
           "mechanism": {}}

    print(f"CLIMATE MODEL {climate_model}  {lo}-{hi}  vs  ANCHOR {A['anchor_period']}")
    if res["inside_anchor"]:
        print("  NOTE the climate model lies inside the anchor: self-consistency check, "
              "not an independent validation.")

    # ------------------------------------------------------------------ tier 1
    print(f"\nTIER 1 - individual years vs the +/-{A['k_screen']:.0f} sigma band")
    for key, ykey, label, unit, dp in SERIES:
        a = A["diagnostics"][key]
        y, v = cm[key]
        band = a["screening_band"]
        z = (v - a["mean"]) / a["sigma_used"]
        outside = np.abs(z) > a["k_screen"]
        nb = len(v) // BLOCK
        blocks_failed = sum(1 for b in range(nb)
                            if outside[b * BLOCK:(b + 1) * BLOCK].any())
        res["tier1"][key] = dict(
            units=unit, band=band, n_years=int(len(v)),
            n_outside=int(outside.sum()),
            years_outside=[int(t) for t in y[outside]],
            worst_sigma=float(z[np.argmax(np.abs(z))]),
            worst_year=int(y[np.argmax(np.abs(z))]),
            n_blocks=int(nb), blocks_failed=int(blocks_failed),
            passes=bool(not outside.any()))
        r = res["tier1"][key]
        print(f"  {label:20s} {r['n_outside']:2d}/{r['n_years']:2d} outside   "
              f"worst {r['worst_sigma']:+5.2f}s ({r['worst_year']})   "
              f"blocks {blocks_failed}/{nb}   "
              f"{'PASS' if r['passes'] else 'FAIL'}")

    # ------------------------------------------------------------------ tier 2
    print(f"\nTIER 2 - rollout mean and variance")
    for key, ykey, label, unit, dp in SERIES:
        a = A["diagnostics"][key]
        y, v = cm[key]
        T = C.window_stats(y, v, lo, hi)
        sig, mu, n = a["sigma_used"], a["mean"], T["n"]
        tol = C.bias_target(sig, n)
        off = T["mean"] - mu
        rlo, rhi = C.ratio_window(n, a["n"])
        ratio = T["sigma_used"] / sig
        res["tier2_mean"][key] = dict(
            units=unit, anchor_mean=mu, sigma=sig, n=n, tolerance=float(tol),
            tolerance_in_sigma=float(tol / sig),
            binding_branch=("0.5 sigma" if 0.5 * sig >= 1.96 * sig / np.sqrt(n)
                            else "detection"),
            climate_model_mean=T["mean"], offset=float(off),
            offset_in_sigma=float(off / sig), passes=bool(abs(off) <= tol))
        res["tier2_variance"][key] = dict(
            climate_model_sigma=T["sigma_used"], climate_model_detrended=T["detrended"],
            ratio=float(ratio), window=[rlo, rhi],
            passes=bool(rlo <= ratio <= rhi))
        m, vr = res["tier2_mean"][key], res["tier2_variance"][key]
        print(f"  {label:20s} mean {m['offset_in_sigma']:+5.2f}s / "
              f"{m['tolerance_in_sigma']:.2f}s  {'PASS' if m['passes'] else 'FAIL'}"
              f"   sigma ratio {vr['ratio']:.2f} "
              f"[{vr['window'][0]:.2f},{vr['window'][1]:.2f}] "
              f"{'pass' if vr['passes'] else 'FAIL'}")

    # daily DJF sigma ratio
    d = A["daily_DJF_u60N"]
    if daily is not None and len(daily) > 0:
        ne, nr = len(daily) / d["tau_days"], d["effective_n"]
        rr = np.sqrt(1 / (2 * ne) + 1 / (2 * nr))
        ratio = float(daily.std(ddof=1) / d["sigma"])
        res["tier2_variance"]["daily_DJF_u60N"] = dict(
            climate_model_sigma=float(daily.std(ddof=1)), anchor_sigma=d["sigma"],
            ratio=ratio, window=[float(1 - 1.96 * rr), float(1 + 1.96 * rr)],
            effective_n=float(ne),
            passes=bool(abs(ratio - 1) <= 1.96 * rr),
            climate_model_p05=float(np.percentile(daily, 5)),
            climate_model_p95=float(np.percentile(daily, 95)))
        dv = res["tier2_variance"]["daily_DJF_u60N"]
        print(f"  {'daily DJF u 60N':20s} sigma ratio {ratio:.2f} "
              f"[{dv['window'][0]:.2f},{dv['window'][1]:.2f}] "
              f"{'pass' if dv['passes'] else 'FAIL'}")
    else:
        res["tier2_variance"]["daily_DJF_u60N"] = dict(
            available=False,
            reason="window does not line up with whole CESM segments")
        print(f"  {'daily DJF u 60N':20s} not available - window is not "
              f"segment-aligned")

    # ------------------------------------------------------------- shape checks
    # Reported at tier 1, gated at tier 2. The tolerance is always taken at the
    # climate model's own n, as for the mean above.
    print(f"\nSHAPE CHECKS - seasonal cycle, daily distribution, w* profile")

    res["shape_seasonal"] = {}
    for key, ykey, label, unit, dp in SERIES:
        a = A["seasonal_cycle"][key]
        y, m12 = cm[f"_monthly_{key}"]
        amp = np.array([C.first_harmonic(r)[0] for r in m12])
        pha = np.array([C.first_harmonic(r)[1] for r in m12])
        n = len(amp)

        sa, ma = a["amplitude"]["sigma_used"], a["amplitude"]["mean"]
        tol_a = C.bias_target(sa, n)
        off_a = float(amp.mean() - ma)

        sp, mp = a["phase"]["sigma_used"], a["phase"]["mean"]
        tol_p = C.bias_target(sp, n)
        off_p = float(C.wrap_months(C.circ_mean_months(pha) - mp))

        res["shape_seasonal"][key] = dict(
            units=unit, n=int(n),
            climate_model_climatology=[float(v) for v in m12.mean(0)],
            anchor_climatology=a["monthly_climatology"],
            month_of_max_observed=int(np.argmax(m12.mean(0)) + 1),
            amplitude=dict(anchor=ma, sigma=sa, climate_model=float(amp.mean()),
                           tolerance=float(tol_a), offset=off_a,
                           offset_in_sigma=float(off_a / sa),
                           passes=bool(abs(off_a) <= tol_a)),
            phase=dict(units="months, 0 = mid-January", anchor=mp, sigma=sp,
                       climate_model=C.circ_mean_months(pha),
                       tolerance=float(tol_p), offset=off_p,
                       offset_in_sigma=float(off_p / sp),
                       passes=bool(abs(off_p) <= tol_p)))
        r = res["shape_seasonal"][key]
        print(f"  {label:20s} amp {r['amplitude']['offset_in_sigma']:+5.2f}s "
              f"{'PASS' if r['amplitude']['passes'] else 'FAIL'}   "
              f"phase {r['phase']['offset_in_sigma']:+5.2f}s "
              f"{'PASS' if r['phase']['passes'] else 'FAIL'}")

    res["shape_daily_distribution"] = {}
    ywq, wq = cm["_djf_pctl"]
    for i, q in enumerate(C.DJF_PCTL):
        a = A["daily_distribution"]["percentiles"][f"p{q}"]
        v = wq[:, i]
        n = len(v)
        tol = C.bias_target(a["sigma_used"], n)
        off = float(v.mean() - a["mean"])
        res["shape_daily_distribution"][f"p{q}"] = dict(
            units="m/s", anchor=a["mean"], sigma=a["sigma_used"], n=int(n),
            climate_model=float(v.mean()), tolerance=float(tol), offset=off,
            offset_in_sigma=float(off / a["sigma_used"]),
            passes=bool(abs(off) <= tol))
        r = res["shape_daily_distribution"][f"p{q}"]
        print(f"  {'u 60N DJF p' + str(q):20s} {r['offset_in_sigma']:+5.2f}s / "
              f"{tol / a['sigma_used']:.2f}s  {'PASS' if r['passes'] else 'FAIL'}")

    res["shape_w_star_profile"] = dict(
        normalisation=A["w_star_profile"]["normalisation"],
        absolute_is_advisory=True, levels={})
    mat = np.array([cm[f"_profile_{pl:g}"][1] for pl in C.PROFILE_LEVELS])
    norm = mat / mat.mean(axis=0, keepdims=True)
    for i, pl in enumerate(C.PROFILE_LEVELS):
        a = A["w_star_profile"]["levels"][f"{pl:g}"]
        an, aa = a["normalised_gated"], a["absolute_advisory"]
        n = mat.shape[1]
        tol = C.bias_target(an["sigma_used"], n)
        off = float(norm[i].mean() - an["mean"])
        res["shape_w_star_profile"]["levels"][f"{pl:g}"] = dict(
            normalised=dict(anchor=an["mean"], sigma=an["sigma_used"], n=int(n),
                            climate_model=float(norm[i].mean()),
                            tolerance=float(tol), offset=off,
                            offset_in_sigma=float(off / an["sigma_used"]),
                            passes=bool(abs(off) <= tol)),
            absolute_advisory=dict(units="mm/s", anchor=aa["mean"],
                                   climate_model=float(mat[i].mean()),
                                   offset=float(mat[i].mean() - aa["mean"])))
        r = res["shape_w_star_profile"]["levels"][f"{pl:g}"]["normalised"]
        print(f"  {'w* ' + f'{pl:g}' + ' hPa norm':20s} {r['offset_in_sigma']:+5.2f}s / "
              f"{tol / an['sigma_used']:.2f}s  {'PASS' if r['passes'] else 'FAIL'}")

    # ------------------------------------------------------------------ counts
    s = A["ssw_NH"]
    k = int(((ssw_years >= lo) & (ssw_years <= hi)).sum())
    exp = s["rate_per_winter"] * winters
    ilo, ihi = int(stats.poisson.ppf(0.025, exp)), int(stats.poisson.ppf(0.975, exp))
    res["counts"]["ssw_NH"] = dict(
        climate_model_count=k, climate_model_winters=int(winters),
        anchor_rate=s["rate_per_winter"], expected=float(exp),
        interval=[ilo, ihi], passes=bool(ilo <= k <= ihi))
    print(f"\n  major NH SSW: {k} in {winters} winters, expected {exp:.1f} "
          f"[{ilo}, {ihi}]  {'PASS' if res['counts']['ssw_NH']['passes'] else 'FAIL'}")

    # -------------------------------------------------------------- mechanism
    for tag, xk, yk, ykey in MECHANISM:
        M = A["mechanism"][tag]
        x, yv = cm[xk][1], cm[yk][1]
        n = min(len(x), len(yv))
        b = float(np.polyfit(x[:n], yv[:n], 1)[0])
        res["mechanism"][tag] = dict(
            anchor_slope=M["slope"], anchor_ci95=M["ci95"],
            climate_model_slope=b, n=int(n),
            passes=bool(M["ci95"][0] <= b <= M["ci95"][1]))
        r = res["mechanism"][tag]
        print(f"  {tag:20s} slope {b:+.3f} vs anchor {M['slope']:+.3f} "
              f"[{M['ci95'][0]:+.3f},{M['ci95'][1]:+.3f}]  "
              f"{'PASS' if r['passes'] else 'FAIL'}")

    # ------------------------------------------------------------------ totals
    t1 = [v["passes"] for v in res["tier1"].values()]
    t2m = [v["passes"] for v in res["tier2_mean"].values()]
    t2v = [v["passes"] for v in res["tier2_variance"].values() if "passes" in v]
    shape = ([v["amplitude"]["passes"] for v in res["shape_seasonal"].values()]
             + [v["phase"]["passes"] for v in res["shape_seasonal"].values()]
             + [v["passes"] for v in res["shape_daily_distribution"].values()]
             + [v["normalised"]["passes"]
                for v in res["shape_w_star_profile"]["levels"].values()])
    res["summary"] = dict(
        tier1_pass=int(sum(t1)), tier1_total=len(t1),
        tier2_mean_pass=int(sum(t2m)), tier2_mean_total=len(t2m),
        tier2_variance_pass=int(sum(t2v)), tier2_variance_total=len(t2v),
        ssw_pass=res["counts"]["ssw_NH"]["passes"],
        mechanism_pass=int(sum(v["passes"] for v in res["mechanism"].values())),
        mechanism_total=len(res["mechanism"]),
        shape_pass=int(sum(shape)), shape_total=len(shape))
    S = res["summary"]
    print(f"\n  tier 1 {S['tier1_pass']}/{S['tier1_total']}   "
          f"tier 2 mean {S['tier2_mean_pass']}/{S['tier2_mean_total']}   "
          f"tier 2 variance {S['tier2_variance_pass']}/{S['tier2_variance_total']}   "
          f"mechanism {S['mechanism_pass']}/{S['mechanism_total']}   "
          f"shape {S['shape_pass']}/{S['shape_total']}")

    with open(out_j, "w") as f:
        json.dump(res, f, indent=2)
    write_markdown(res, A, out_m)
    print(f"\nwrote {out_j}\nwrote {out_m}")


def verdict(ok):
    return "PASS" if ok else "**FAIL**"


def write_markdown(res, A, out_m):
    lo, hi = res["climate_model_period"]
    S = res["summary"]
    L = []
    w = L.append

    w("# Validation result")
    w("")
    w("| | |")
    w("|---|---|")
    w(f"| Climate model | {res['climate_model']} |")
    w(f"| Scored period | {lo}-{hi} |")
    w(f"| Produced | {res['produced']} |")
    w(f"| Anchor | CESM {res['anchor_period']}, "
      f"{A['diagnostics']['vortex_SH']['n']} JJA / "
      f"{A['diagnostics']['vortex_NH']['n']} DJF seasons, "
      f"{A['diagnostics']['mass_flux']['n']} annual years |")
    w(f"| Estimator | `{res['estimator']}` |")
    w(f"| Anchor σ | {A['sigma_rule']} |")
    w(f"| Independent of anchor | {'no — the model lies inside the anchor' if res['inside_anchor'] else 'yes'} |")
    w("| Generated by | `scripts/17_validate.py`, from "
      "`output/16_anchors_45yr.json` |")
    w("")
    w(f"| Tier | Result |")
    w("|---|---|")
    w(f"| 1 · screening, individual years | {S['tier1_pass']}/{S['tier1_total']} |")
    w(f"| 2 · mean | {S['tier2_mean_pass']}/{S['tier2_mean_total']} |")
    w(f"| 2 · variance | {S['tier2_variance_pass']}/{S['tier2_variance_total']} |")
    w(f"| 2 · SSW count | {verdict(S['ssw_pass'])} |")
    w(f"| 2 · mechanism slopes | {S['mechanism_pass']}/{S['mechanism_total']} |")
    w(f"| 2 · shape — cycle, distribution, profile | "
      f"{S['shape_pass']}/{S['shape_total']} |")
    w("")
    w(f"## Tier 1 — individual years, ±{A['k_screen']:.0f}σ")
    w("")
    w("| Diagnostic | Unit | Accept | Years | Outside | Worst | 5-yr blocks failed | Verdict |")
    w("|---|---|---|---|---|---|---|---|")
    for key, _, label, unit, dp in SERIES:
        r, a = res["tier1"][key], A["diagnostics"][key]
        w(f"| {label} | {unit} | {a['screening_band'][0]:.{dp}f} – "
          f"{a['screening_band'][1]:.{dp}f} | {r['n_years']} | {r['n_outside']} | "
          f"{r['worst_sigma']:+.2f}σ ({r['worst_year']}) | "
          f"{r['blocks_failed']}/{r['n_blocks']} | {verdict(r['passes'])} |")
    w("")
    w("## Tier 2 — mean")
    w("")
    w("| Diagnostic | Unit | Anchor mean | σ | Tolerance | Climate model | Offset | Verdict |")
    w("|---|---|---|---|---|---|---|---|")
    for key, _, label, unit, dp in SERIES:
        m = res["tier2_mean"][key]
        w(f"| {label} | {unit} | {m['anchor_mean']:.{dp}f} | {m['sigma']:.{dp}f} | "
          f"±{m['tolerance']:.{dp}f} ({m['tolerance_in_sigma']:.2f}σ) | "
          f"{m['climate_model_mean']:.{dp}f} | {m['offset_in_sigma']:+.2f}σ | "
          f"{verdict(m['passes'])} |")
    w("")
    w("## Tier 2 — variance")
    w("")
    w("| Metric | Anchor σ | Climate-model σ | Ratio | 95% window | Verdict |")
    w("|---|---|---|---|---|---|")
    for key, _, label, unit, dp in SERIES:
        v, a = res["tier2_variance"][key], A["diagnostics"][key]
        w(f"| {label} | {a['sigma_used']:.{dp}f} | {v['climate_model_sigma']:.{dp}f} | "
          f"{v['ratio']:.2f} | {v['window'][0]:.2f} – {v['window'][1]:.2f} | "
          f"{verdict(v['passes'])} |")
    d = res["tier2_variance"]["daily_DJF_u60N"]
    if d.get("available", True):
        w(f"| Daily DJF u, 60°N | {d['anchor_sigma']:.2f} | "
          f"{d['climate_model_sigma']:.2f} | {d['ratio']:.2f} | "
          f"{d['window'][0]:.2f} – {d['window'][1]:.2f} | {verdict(d['passes'])} |")
    else:
        w(f"| Daily DJF u, 60°N | — | — | — | — | n/a: {d['reason']} |")
    w("")
    w("## Tier 2 — shape")
    w("")
    w("Three things the mean and variance tests cannot see: the annual march, the "
      "shape of the daily distribution, and the vertical structure of the tropical "
      "upwelling.")
    w("")
    w("### Seasonal cycle — annual harmonic")
    w("")
    w("Amplitude is the harmonic's half range. Phase is its month of maximum, in "
      "months from mid-January, compared circularly. The twelve monthly means are "
      "not scored individually — that would be a twelve-way multiplicity problem.")
    w("")
    w("| Diagnostic | Unit | Amp anchor | Amp offset | Amp | Phase anchor | "
      "Phase offset | Phase |")
    w("|---|---|---|---|---|---|---|---|")
    for key, _, label, unit, dp in SERIES:
        r = res["shape_seasonal"][key]
        a, p = r["amplitude"], r["phase"]
        w(f"| {label} | {unit} | {a['anchor']:.{dp}f} | "
          f"{a['offset_in_sigma']:+.2f}σ | {verdict(a['passes'])} | "
          f"{p['anchor']:.2f} | {p['offset_in_sigma']:+.2f}σ | "
          f"{verdict(p['passes'])} |")
    w("")
    w("### Daily distribution — per-winter percentiles of u at 60°N")
    w("")
    w("Each percentile is taken **within** each DJF winter, so the sample is winters "
      "and the 14-day decorrelation time of the daily series never enters. A squashed "
      "model shows p5 too high and p95 too low at once, which no mean test sees.")
    w("")
    w("| Percentile | Anchor | σ | Tolerance | Climate model | Offset | Verdict |")
    w("|---|---|---|---|---|---|---|")
    for q in C.DJF_PCTL:
        r = res["shape_daily_distribution"][f"p{q}"]
        w(f"| p{q} | {r['anchor']:.2f} | {r['sigma']:.2f} | ±{r['tolerance']:.2f} | "
          f"{r['climate_model']:.2f} | {r['offset_in_sigma']:+.2f}σ | "
          f"{verdict(r['passes'])} |")
    w("")
    w("### Tropical w* profile — 10°S–10°N")
    w("")
    w("The gate is on the profile divided by its own vertical mean, which cancels a "
      "multiplicative estimator bias that is uniform in height. The absolute values "
      "are **advisory**: appendix C records a grid error on w* at 70 hPa larger than "
      "the tier-2 tolerance, so an absolute per-level target would fail a correct "
      "model on grid choice alone. Height-uniformity is assumed, not measured.")
    w("")
    w("| Level | Normalised anchor | Tolerance | Normalised model | Offset | "
      "Verdict | Absolute anchor (advisory) | Absolute model |")
    w("|---|---|---|---|---|---|---|---|")
    for pl in C.PROFILE_LEVELS:
        r = res["shape_w_star_profile"]["levels"][f"{pl:g}"]
        nm, ab = r["normalised"], r["absolute_advisory"]
        w(f"| {pl:g} hPa | {nm['anchor']:.4f} | ±{nm['tolerance']:.4f} | "
          f"{nm['climate_model']:.4f} | {nm['offset_in_sigma']:+.2f}σ | "
          f"{verdict(nm['passes'])} | {ab['anchor']:.4f} mm/s | "
          f"{ab['climate_model']:.4f} mm/s |")
    w("")
    w("## Tier 2 — counts and relations")
    w("")
    w("| Check | Anchor | Climate model | Accept | Verdict |")
    w("|---|---|---|---|---|")
    s = res["counts"]["ssw_NH"]
    w(f"| Major NH SSW, count | {s['anchor_rate']:.2f}/winter | "
      f"{s['climate_model_count']} in {s['climate_model_winters']} winters | "
      f"{s['interval'][0]}–{s['interval'][1]} (expect {s['expected']:.1f}) | "
      f"{verdict(s['passes'])} |")
    for tag in res["mechanism"]:
        m = res["mechanism"][tag]
        w(f"| {tag} | {m['anchor_slope']:+.3f} | {m['climate_model_slope']:+.3f} | "
          f"{m['anchor_ci95'][0]:+.3f} – {m['anchor_ci95'][1]:+.3f} | "
          f"{verdict(m['passes'])} |")
    w("")
    w("## Flags")
    w("")
    w("| Check | Reading |")
    w("|---|---|")
    any_flag = False
    for key, _, label, unit, dp in SERIES:
        if not res["tier1"][key]["passes"]:
            r = res["tier1"][key]
            any_flag = True
            shared = sorted(set(r["years_outside"])
                            & set(A["diagnostics"][key]["years_outside_band"]))
            note = (f"; {len(shared)} of them the anchor also flags against its own "
                    f"record ({', '.join(str(t) for t in shared)})" if shared else "")
            w(f"| Tier 1, {label} | {r['n_outside']} of {r['n_years']} years outside; "
              f"worst {r['worst_sigma']:+.2f}σ in {r['worst_year']}{note} |")
        if not res["tier2_mean"][key]["passes"]:
            m = res["tier2_mean"][key]
            any_flag = True
            w(f"| Tier 2 mean, {label} | offset {m['offset_in_sigma']:+.2f}σ "
              f"against a {m['tolerance_in_sigma']:.2f}σ tolerance |")
        if not res["tier2_variance"][key]["passes"]:
            v = res["tier2_variance"][key]
            any_flag = True
            w(f"| Tier 2 variance, {label} | ratio {v['ratio']:.2f}, outside "
              f"{v['window'][0]:.2f}–{v['window'][1]:.2f} |")
    if not res["counts"]["ssw_NH"]["passes"]:
        s_ = res["counts"]["ssw_NH"]
        any_flag = True
        w(f"| Tier 2 SSW count | {s_['climate_model_count']} in "
          f"{s_['climate_model_winters']} winters, outside "
          f"{s_['interval'][0]}–{s_['interval'][1]} |")
    for tag, m in res["mechanism"].items():
        if not m["passes"]:
            any_flag = True
            w(f"| {tag} | slope {m['climate_model_slope']:+.3f}, outside the anchor CI "
              f"{m['anchor_ci95'][0]:+.3f} – {m['anchor_ci95'][1]:+.3f} |")
    if not any_flag:
        w("| — | no check failed |")
    w("")
    w("| Standing condition | Value |")
    w("|---|---|")
    w(f"| Anchor carries the forced BDC trend | mass flux "
      f"{A['diagnostics']['mass_flux']['trend_per_decade']:+.4f} per decade, "
      f"p = {A['diagnostics']['mass_flux']['trend_p']:.0e} |")
    w(f"| Anchor years outside their own band | "
      f"{A['band_self_consistency']['outside']} of "
      f"{A['band_self_consistency']['checks']} |")
    w(f"| Climate model inside anchor | "
      f"{'yes' if res['inside_anchor'] else 'no'} |")
    w("")
    w("## Evidence")
    w("")
    w("| Figure | Covers |")
    w("|---|---|")
    for fn, cap in FIGURES:
        sfn = stamped(fn, res["stamp"])
        w(f"| [{sfn}]({sfn}) | {cap} |")
    w("")
    w("Figures are written by `scripts/18_validation_figures.py`, from the same two "
      "JSON files as the tables above.")
    w("")

    with open(out_m, "w") as f:
        f.write("\n".join(L))


if __name__ == "__main__":
    main()
