"""
23 - Score the ERA5 monthly series against tier 1, and draw it.

RESTRICTED SCORECARD. Tier 1 gates six diagnostics; a monthly archive can
supply four. Mass flux and w* need the eddy heat flux v'theta', which cannot
be formed once the fields have been time-averaged, and the SSW count and daily
DJF sigma need the daily series. Those are reported here as NOT EVALUABLE
rather than omitted, so the figure cannot be mistaken for a full tier-1 pass.

The four that are scored go through the same band test 17_validate.py applies:
z = (value - anchor mean) / anchor sigma, flagged when |z| > k_screen, with the
5-year block count alongside. Thresholds come from output/16_anchors_45yr.json,
untouched - this scores ERA5 against the CESM anchor as it stands.

Run with the pinned analysis environment, after 22:
  ../../.AIDE-eval_env/bin/python 23_era5_monthly_tier1.py [--climate-model NAME]

Output: output/23_era5_tier1__<stamp>.json
        validation_results/era5_tier1_screening__<stamp>.png
        validation_results/era5_tier1_sigma__<stamp>.png
        validation_results/era5_tier1__<stamp>.md
"""
import os
import re
import sys
import json
import argparse
import datetime
import importlib.util

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import aide_val_common as C
from report_layout import SURFACE, INK, INK2, MUTED, RULE, BLUE, ORANGE, AQUA

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RESULTS = os.path.join(ROOT, "validation_results")
BLOCK = 5

# Reuse 18's figure furniture rather than restate it, so these plots sit in the
# same visual system as the rest of the suite. The module name starts with a
# digit, so it cannot be imported by name.
_spec = importlib.util.spec_from_file_location(
    "figs18", os.path.join(HERE, "18_validation_figures.py"))
figs = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(figs)
import matplotlib.pyplot as plt

SCORED = [
    ("vortex_NH", "_djf_years", "Vortex NH, DJF", "m s⁻¹", 2),
    ("vortex_SH", "_jja_years", "Vortex SH, JJA", "m s⁻¹", 2),
    ("polar_cap_T_NH", "_djf_years", "Polar cap T, NH", "K", 2),
    ("polar_cap_T_SH", "_jja_years", "Polar cap T, SH", "K", 2),
]
NOT_EVALUABLE = [
    ("Mass flux, 70 hPa", "needs v'θ′, which cannot be formed after time averaging"),
    ("w̄*, 10°S–10°N", "needs v'θ′"),
    ("Major NH SSW count", "needs the day the wind reverses"),
    ("Daily DJF σ ratio", "needs the daily series"),
]


def score(A, series, years_of):
    out = {}
    for key, ykey, label, unit, dp in SCORED:
        a = A["diagnostics"][key]
        v = np.asarray(series[key], dtype=float)
        y = np.asarray(years_of[ykey], dtype=int)
        z = (v - a["mean"]) / a["sigma_used"]
        outside = np.abs(z) > a["k_screen"]
        nb = len(v) // BLOCK
        blocks_failed = sum(1 for b in range(nb)
                            if outside[b * BLOCK:(b + 1) * BLOCK].any())
        out[key] = dict(
            units=unit, band=a["screening_band"], anchor_mean=a["mean"],
            anchor_sigma=a["sigma_used"], k_screen=a["k_screen"],
            years=[int(t) for t in y], values=[float(x) for x in v],
            sigma=[float(x) for x in z],
            n_years=int(len(v)), n_outside=int(outside.sum()),
            years_outside=[int(t) for t in y[outside]],
            worst_sigma=float(z[np.argmax(np.abs(z))]),
            worst_year=int(y[np.argmax(np.abs(z))]),
            n_blocks=int(nb), blocks_failed=int(blocks_failed),
            passes=bool(not outside.any()))
    return out


def fig_screening(A, res, ps):
    fig = figs.newfig(
        9.6, 6.0, "Tier 1 (restricted) — ERA5 monthly against the ±3σ band",
        "Grey: CESM anchor seasons. Blue: ERA5. Band: anchor mean ±3σ. Only the "
        "four diagnostics a monthly archive can supply are shown.",
        "a season outside the band is a flag, not a verdict; two of the six "
        "tier-1 diagnostics are absent entirely.")
    for i, (key, ykey, label, unit, dp) in enumerate(SCORED):
        ax = figs.style(fig.add_subplot(2, 2, i + 1))
        a, r = A["diagnostics"][key], res["tier1"][key]
        band, mu = a["screening_band"], a["mean"]
        ax.axhspan(band[0], band[1], color=BLUE, alpha=0.06, lw=0)
        for b in band:
            ax.axhline(b, color=BLUE, lw=0.8, ls=(0, (4, 3)), alpha=0.7)
        ax.axhline(mu, color=MUTED, lw=0.7)
        ay, av = C.join_segments(ps["series"], [ps["test"], ps["train"]], key, ykey)
        ax.plot(ay, av, "o", ms=2.4, color=RULE, mec="none")
        y = np.array(r["years"]); v = np.array(r["values"])
        ax.plot(y, v, "o", ms=4.2, color=BLUE, mec="none")
        bad = np.isin(y, r["years_outside"])
        if bad.any():
            ax.plot(y[bad], v[bad], "o", ms=8, mfc="none", mec=ORANGE, mew=1.4)
            for t, val in zip(y[bad], v[bad]):
                ax.annotate(f"{int(t)}", (t, val), textcoords="offset points",
                            xytext=(7, 0), color=ORANGE, fontsize=6.5, va="center")
        figs.style(ax, f"{'abcd'[i]}   {label}",
                   f"{unit} · {r['n_outside']}/{r['n_years']} outside · "
                   f"worst {r['worst_sigma']:+.2f}σ ({r['worst_year']}) · "
                   f"{'PASS' if r['passes'] else 'FAIL'}")
        ax.set_xlim(1968, 2016)
        span = [min(band[0], av.min(), v.min()), max(band[1], av.max(), v.max())]
        pad = 0.10 * (span[1] - span[0])
        ax.set_ylim(span[0] - pad, span[1] + pad)
    fig.subplots_adjust(left=0.085, right=0.975, top=fig.text_bottom - 0.070,
                        bottom=0.105, hspace=0.55, wspace=0.21)
    figs.footer(fig, res)
    figs.save(fig, "era5_tier1_screening.png", res["stamp"])


def fig_sigma(A, res):
    fig = figs.newfig(
        9.6, 4.4, "Tier 1 (restricted) — every ERA5 season in units of anchor σ",
        "Each marker is one ERA5 season. The gate is ±3σ; the shaded band is ±1σ "
        "for scale.",
        "this is the same data as the previous figure, reduced to the quantity "
        "the gate actually tests.")
    ax = figs.style(fig.add_subplot(111))
    ax.axhspan(-1, 1, color=RULE, alpha=0.45, lw=0)
    ax.axhline(0, color=MUTED, lw=0.8)
    for k in (-3, 3):
        ax.axhline(k, color=ORANGE, lw=0.9, ls=(0, (4, 3)))
    ax.text(len(SCORED) - 0.45, 3.05, "±3σ gate", color=ORANGE, fontsize=7,
            ha="right", va="bottom")
    for i, (key, ykey, label, unit, dp) in enumerate(SCORED):
        r = res["tier1"][key]
        z = np.array(r["sigma"])
        x = i + np.linspace(-0.17, 0.17, len(z))
        ok = np.abs(z) <= r["k_screen"]
        ax.plot(x[ok], z[ok], "o", ms=5, color=AQUA, mec="none")
        if (~ok).any():
            ax.plot(x[~ok], z[~ok], "o", ms=6.5, color=ORANGE, mec="none")
        ax.plot([i - 0.28, i + 0.28], [z.mean()] * 2, "-", color=INK2, lw=1.4)
    ax.set_xticks(range(len(SCORED)))
    ax.set_xticklabels([s[2] for s in SCORED], fontsize=8, color=INK)
    ax.set_ylabel("deviation from anchor mean, in anchor σ", color=INK2, fontsize=8)
    lim = max(3.6, 1.15 * max(abs(np.array(res["tier1"][k]["sigma"])).max()
                              for k, *_ in SCORED))
    ax.set_ylim(-lim, lim)
    ax.set_xlim(-0.55, len(SCORED) - 0.45)
    ax.plot([], [], "-", color=INK2, lw=1.4, label="series mean")
    lg = ax.legend(frameon=False, fontsize=7.5, loc="lower left",
                   bbox_to_anchor=(0.0, 1.01))
    for t in lg.get_texts():
        t.set_color(INK2)
    fig.subplots_adjust(left=0.095, right=0.975, top=fig.text_bottom - 0.085,
                        bottom=0.135)
    figs.footer(fig, res)
    figs.save(fig, "era5_tier1_sigma.png", res["stamp"])


def report(res, path):
    L = []
    w = L.append
    w(f"# Tier 1 (restricted) — {res['climate_model']}")
    w("")
    w("| | |")
    w("|---|---|")
    w(f"| Source | {res['source']} |")
    w(f"| Scored | {res['climate_model_period'][0]}–{res['climate_model_period'][1]} |")
    w(f"| Anchor | CESM {res['anchor_period']} |")
    w(f"| Produced | {res['produced']} |")
    w(f"| Scope | 4 of 6 tier-1 diagnostics |")
    w("")
    w("## Scored")
    w("")
    w("| Diagnostic | Unit | Accept | Seasons | Outside | Worst | 5-yr blocks | Verdict |")
    w("|---|---|---|---|---|---|---|---|")
    for key, ykey, label, unit, dp in SCORED:
        r = res["tier1"][key]
        w(f"| {label} | {unit} | {r['band'][0]:.{dp}f} – {r['band'][1]:.{dp}f} "
          f"| {r['n_years']} | {r['n_outside']} | "
          f"{r['worst_sigma']:+.2f}σ ({r['worst_year']}) | "
          f"{r['blocks_failed']}/{r['n_blocks']} | "
          f"{'PASS' if r['passes'] else '**FAIL**'} |")
    w("")
    w("## Not evaluable from a monthly archive")
    w("")
    w("| Diagnostic | Why |")
    w("|---|---|")
    for label, why in NOT_EVALUABLE:
        w(f"| {label} | {why} |")
    w("")
    w("## Caveats")
    w("")
    w("- ERA5 is a reanalysis of the real atmosphere; the anchor is fixed-SST "
      "CESM with no ENSO. Differences here are not by themselves an ERA5 error.")
    w("- The polar-cap layer uses "
      f"{res['n_levels_in_polar_cap_layer']} ERA5 levels inside 10–50 hPa "
      "against CESM's 8.")
    w("- Seasonal means from monthly data are exact, not approximate: every "
      "reduction involved is linear, and the day weighting is preserved.")
    w("")
    open(path, "w").write("\n".join(L) + "\n")
    print(f"wrote {path}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--climate-model", dest="climate_model", default=None)
    args = ap.parse_args()

    S = json.load(open(os.path.join(C.OUTDIR, "22_era5_monthly_series.json")))
    A = json.load(open(os.path.join(C.OUTDIR, "16_anchors_45yr.json")))
    ps = json.load(open(os.path.join(C.OUTDIR, "07_period_split.json")))

    lo, hi = S["years"]
    name = args.climate_model or f"ERA5 monthly {lo}-{hi}"
    produced = datetime.date.today().isoformat()
    stamp = f"{re.sub(r'[^A-Za-z0-9.-]+', '-', name).strip('-')}__{produced}"

    res = dict(
        climate_model=name, produced=produced, stamp=stamp,
        climate_model_period=[lo, hi], anchor_period=A["anchor_period"],
        estimator="aide_val_common (linear reductions only; no TEM)",
        inside_anchor=bool(lo >= A["anchor_years"][0] and hi <= A["anchor_years"][1]),
        source=S["source"], scope=S["scope"],
        n_levels_in_polar_cap_layer=S["n_levels_in_polar_cap_layer"],
        not_evaluable={k: v for k, v in NOT_EVALUABLE},
        tier1=score(A, S["series"], S["series"]))

    npass = sum(res["tier1"][k]["passes"] for k, *_ in SCORED)
    print(f"\nTIER 1 (restricted)  {name}  vs anchor {A['anchor_period']}")
    for key, ykey, label, unit, dp in SCORED:
        r = res["tier1"][key]
        print(f"  {label:18s} {r['n_outside']}/{r['n_years']} outside  "
              f"worst {r['worst_sigma']:+5.2f}s ({r['worst_year']})  "
              f"{'PASS' if r['passes'] else 'FAIL'}")
    print(f"  -> {npass}/{len(SCORED)} scored diagnostics pass; "
          f"{len(NOT_EVALUABLE)} not evaluable from monthly data")

    os.makedirs(RESULTS, exist_ok=True)
    out_j = os.path.join(C.OUTDIR, f"23_era5_tier1__{stamp}.json")
    json.dump(res, open(out_j, "w"), indent=2)
    print(f"\nwrote {out_j}")
    fig_screening(A, res, ps)
    fig_sigma(A, res)
    report(res, os.path.join(RESULTS, f"era5_tier1__{stamp}.md"))


if __name__ == "__main__":
    main()
