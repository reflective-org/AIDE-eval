"""
18 - The figures behind validation_results/validation_result.md.

One figure per test family, so every verdict in the report has something to look
at rather than a number to trust:

  tier1_screening.png       every climate-model year against the +/-3 sigma band
  tier2_mean.png            offset from the anchor against the +/-0.5 sigma tolerance
  tier2_variance.png        sigma ratios against their 95% windows
  counts_and_relations.png  the major NH SSW count against its Poisson interval

Reads output/16_anchors_45yr__<stamp>.json, written by 17, plus the series in
output/07_period_split.json. Run after 17.

  18_validation_figures.py [--climate-model NAME]

With no argument it takes the most recently written 17_validation__*.json and
says which one. The stamp comes out of that file, so the figure names and the
report that links them cannot disagree.

Output: PNGs in validation_results/, stamped `__<climate model>__<date>`
"""

import argparse
import glob
import json
import os
import re
import textwrap
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import aide_val_common as C
from report_layout import SURFACE, INK, INK2, MUTED, RULE, BLUE, ORANGE, AQUA

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "validation_results")

SERIES = [
    ("mass_flux", "_annual_years", "Mass flux, 70 hPa", "10⁹ kg s⁻¹"),
    ("w_star", "_annual_years", "w̄*, 10°S–10°N", "mm s⁻¹"),
    ("vortex_NH", "_djf_years", "Vortex NH, DJF", "m s⁻¹"),
    ("vortex_SH", "_jja_years", "Vortex SH, JJA", "m s⁻¹"),
    ("polar_cap_T_NH", "_djf_years", "Polar cap T, NH", "K"),
    ("polar_cap_T_SH", "_jja_years", "Polar cap T, SH", "K"),
]
FAIL, PASS = ORANGE, AQUA


def style(ax, title=None, sub=None):
    ax.set_facecolor(SURFACE)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(RULE)
    ax.tick_params(colors=MUTED, labelsize=7, length=3)
    if title:
        ax.set_title(title, loc="left", color=INK, fontsize=8.5,
                     pad=16 if sub else 5)
    if sub:
        ax.text(0, 1.035, sub, transform=ax.transAxes, color=MUTED, fontsize=7)
    return ax


def newfig(w, h, title, lead, interp=None):
    """Title, a descriptive lead, and optionally an interpretation line.

    The lead says what is plotted; the interpretation line says what it means,
    kept visually separate the way the protocol separates the two.
    `fig.text_bottom` is where the block ends, so the axes can sit under it.
    """
    fig = plt.figure(figsize=(w, h), facecolor=SURFACE)
    fig.text(0.045, 0.965, title, color=INK, fontsize=12.5, va="top")
    wrap = int(w * 13.5)
    lead_t = textwrap.fill(lead, wrap)
    fig.text(0.045, 0.925, lead_t, color=INK2, fontsize=8, va="top", linespacing=1.5)

    def block_h(text, size):
        return (text.count("\n") + 1) * size * 1.5 / (h * 72)

    y = 0.925 - block_h(lead_t, 8)
    if interp:
        y -= 0.010
        interp_t = textwrap.fill(f"Interpretation — {interp}", wrap)
        fig.text(0.045, y, interp_t, color=MUTED, fontsize=7.5, style="italic",
                 va="top", linespacing=1.5)
        y -= block_h(interp_t, 7.5)
    fig.text_bottom = y
    return fig


def footer(fig, res):
    lo, hi = res["climate_model_period"]
    # two lines: a long climate-model name would run off the page on one
    fig.text(0.045, 0.033,
             f"{res['climate_model']}  ·  scored {lo}–{hi}  ·  "
             f"produced {res['produced']}",
             color=MUTED, fontsize=6.5)
    fig.text(0.045, 0.014,
             f"anchor CESM {res['anchor_period']}  ·  "
             f"estimator {res['estimator']}  ·  "
             + ("model lies INSIDE the anchor: self-consistency check"
                if res["inside_anchor"] else "model independent of the anchor"),
             color=MUTED, fontsize=6.5)


def save(fig, name, stamp):
    root, ext = os.path.splitext(name)
    p = os.path.join(RESULTS, f"{root}__{stamp}{ext}")
    fig.savefig(p, dpi=200, facecolor=SURFACE)
    plt.close(fig)
    print(f"wrote {p}")


# ---------------------------------------------------------------- tier 1
def fig_tier1(A, res, S, segs):
    fig = newfig(11.0, 7.2, "Tier 1 — individual years against the ±3σ band",
                 "Grey: anchor years outside the model window. Blue: climate-model "
                 "years. Band: anchor mean ±3σ.",
                 "a year outside the band is a flag, not a verdict.")
    lo, hi = res["climate_model_period"]
    for i, (key, ykey, label, unit) in enumerate(SERIES):
        ax = style(fig.add_subplot(3, 2, i + 1))
        a, r = A["diagnostics"][key], res["tier1"][key]
        y, v = C.join_segments(S, segs, key, ykey)
        band, mu = a["screening_band"], a["mean"]
        ax.axhspan(band[0], band[1], color=BLUE, alpha=0.06, lw=0)
        for b in band:
            ax.axhline(b, color=BLUE, lw=0.8, ls=(0, (4, 3)), alpha=0.7)
        ax.axhline(mu, color=MUTED, lw=0.7)
        out = (y >= lo) & (y <= hi)
        ax.plot(y[~out], v[~out], "o", ms=2.6, color=RULE, mec="none")
        ax.plot(y[out], v[out], "o", ms=3.4, color=BLUE, mec="none")
        bad = np.isin(y, r["years_outside"])
        if bad.any():
            ax.plot(y[bad], v[bad], "o", ms=7, mfc="none", mec=FAIL, mew=1.4)
            for t, val in zip(y[bad], v[bad]):
                ax.annotate(f"{int(t)}", (t, val), textcoords="offset points",
                            xytext=(7, 0), color=FAIL, fontsize=6.5, va="center")
        v_ok = "PASS" if r["passes"] else "FAIL"
        style(ax, f"{'abcdef'[i]}   {label}",
              f"{unit} · {r['n_outside']}/{r['n_years']} outside · "
              f"worst {r['worst_sigma']:+.2f}σ ({r['worst_year']}) · {v_ok}")
        ax.set_xlim(1968, 2016)
        span = [min(band[0], v.min()), max(band[1], v.max())]
        pad = 0.10 * (span[1] - span[0])
        ax.set_ylim(span[0] - pad, span[1] + pad)
    fig.subplots_adjust(left=0.075, right=0.975, top=fig.text_bottom - 0.062, bottom=0.075,
                        hspace=0.62, wspace=0.19)
    footer(fig, res)
    save(fig, "tier1_screening.png", res["stamp"])


# ---------------------------------------------------------------- tier 2 mean
def fig_mean(A, res):
    fig = newfig(9.2, 5.4, "Tier 2 — the rollout mean",
                 "Offset of the climate-model mean from the anchor, in units of the "
                 "anchor σ. Shaded: the ±0.5σ tolerance. All six diagnostics are "
                 "scaled by σ, so they share one axis.")
    ax = style(fig.add_subplot(111))
    keys = [k for k, *_ in SERIES][::-1]
    labels = [l for _, _, l, _ in SERIES][::-1]
    ax.axvspan(-0.5, 0.5, color=AQUA, alpha=0.10, lw=0)
    for x in (-0.5, 0.5):
        ax.axvline(x, color=AQUA, lw=0.9, ls=(0, (4, 3)))
    ax.axvline(0, color=MUTED, lw=0.8)
    for i, k in enumerate(keys):
        _offset_row(ax, i, res["tier2_mean"][k])
    ax.set_yticks(range(len(keys)))
    ax.set_yticklabels(labels, fontsize=8, color=INK)
    ax.set_xlabel("offset from the anchor mean, in anchor σ", color=INK2, fontsize=8)
    ax.set_ylim(-0.7, len(keys) - 0.3)
    lim = max(1.3, 1.15 * max(abs(res["tier2_mean"][k]["offset_in_sigma"])
                             for k in keys))
    ax.set_xlim(-lim, lim)
    ax.text(0.5, 1.0, " tolerance ±0.5σ", transform=ax.get_xaxis_transform(),
            color=AQUA, fontsize=7, va="top")
    fig.subplots_adjust(left=0.20, right=0.965, top=fig.text_bottom - 0.030, bottom=0.135)
    footer(fig, res)
    save(fig, "tier2_mean.png", res["stamp"])


# ------------------------------------------------------------ tier 2 variance
def fig_variance(A, res):
    # samples per winter comes from 16, not a literal: it follows the anchor length
    eff = A["tiers"]["variance"]["daily_DJF_u60N"]["effective_samples_per_winter"]
    tau = A["daily_DJF_u60N"]["tau_days"]
    fig = newfig(9.2, 5.4, "Tier 2 — the variance",
                 "Ratio of climate-model σ to anchor σ, against the 95% window that "
                 "sampling alone allows.",
                 f"the daily DJF ratio is the sharper test: {eff:.1f} independent "
                 f"samples per winter, at a {tau:.0f}-day decorrelation time, "
                 "against one for the annual series.")
    ax = style(fig.add_subplot(111))
    rows = [(l, res["tier2_variance"][k]) for k, _, l, _ in SERIES]
    d = res["tier2_variance"]["daily_DJF_u60N"]
    if d.get("available", True):
        rows.append(("Daily DJF u, 60°N", d))
    rows = rows[::-1]
    ax.axvline(1.0, color=MUTED, lw=0.8)
    for i, (label, v) in enumerate(rows):
        c = PASS if v["passes"] else FAIL
        w = v["window"]
        ax.plot(w, [i, i], color=RULE, lw=5, solid_capstyle="butt", zorder=1)
        for b in w:
            ax.plot([b, b], [i - 0.22, i + 0.22], color=MUTED, lw=0.9, zorder=2)
        ax.plot(v["ratio"], i, "o", ms=8, color=c, mec="none", zorder=3)
        ax.annotate(f"{v['ratio']:.2f}", (v["ratio"], i),
                    textcoords="offset points", xytext=(0, 11), ha="center",
                    color=c, fontsize=7.5)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([l for l, _ in rows], fontsize=8, color=INK)
    ax.set_xlabel("climate-model σ ÷ anchor σ", color=INK2, fontsize=8)
    ax.set_ylim(-0.7, len(rows) - 0.3)
    ax.text(0.02, 0.02, "grey bar: 95% window from sampling error alone",
            transform=ax.transAxes, color=MUTED, fontsize=7)
    fig.subplots_adjust(left=0.20, right=0.965, top=fig.text_bottom - 0.030, bottom=0.135)
    footer(fig, res)
    save(fig, "tier2_variance.png", res["stamp"])


# ------------------------------------------------------------------- the count
def fig_counts(A, res):
    s = res["counts"]["ssw_NH"]
    fig = newfig(6.6, 4.6, "Tier 2 — the major NH SSW count",
                 "Major sudden stratospheric warmings against the Poisson interval "
                 "implied by the anchor rate. Blue bars are the 95% acceptance "
                 "interval; grey bars fall outside it.",
                 "a count, not a rate: a rate band needs about 23 winters before it "
                 "is resolvable at all, while a count interval is well defined at any "
                 "length (D6).")
    ax = style(fig.add_subplot(111), None,
               f"{s['climate_model_count']} in {s['climate_model_winters']} winters · "
               f"expect {s['expected']:.1f} · anchor rate "
               f"{s['anchor_rate']:.3f}/winter · "
               f"{'PASS' if s['passes'] else 'FAIL'}")
    ilo, ihi = s["interval"]
    xs = np.arange(max(0, ilo - 4), ihi + 5)
    from scipy.stats import poisson
    ax.bar(xs, poisson.pmf(xs, s["expected"]), width=0.85,
           color=[BLUE if ilo <= x <= ihi else RULE for x in xs], lw=0)
    c = PASS if s["passes"] else FAIL
    ax.axvline(s["climate_model_count"], color=c, lw=1.6)
    ax.annotate(f"climate model {s['climate_model_count']}",
                (s["climate_model_count"], 1.0),
                xycoords=("data", "axes fraction"), textcoords="offset points",
                xytext=(4, -8), color=c, fontsize=7.5)
    ax.annotate(f"accept {ilo}–{ihi} events", (0.025, 0.95),
                xycoords="axes fraction", color=BLUE, fontsize=7.5, ha="left")
    ax.set_xlabel("events", color=INK2, fontsize=8)
    ax.set_ylabel("Poisson probability", color=INK2, fontsize=8)
    fig.subplots_adjust(left=0.135, right=0.965, top=fig.text_bottom - 0.075,
                        bottom=0.135)
    footer(fig, res)
    save(fig, "counts_and_relations.png", res["stamp"])


# ------------------------------------------------------------------ shape: cycle
MONTH_INITIALS = "JFMAMJJASOND"


def _offset_row(ax, i, d, fmt="{:+.2f}σ"):
    """One offset marker. Advisory checks are drawn grey and hollow: their
    tolerance is below the n = 15.4 crossover, so neither colour would be honest."""
    off = d["offset_in_sigma"]
    adv = d.get("advisory")
    c = MUTED if adv else (PASS if d["passes"] else FAIL)
    ax.plot([0, off], [i, i], color=c, lw=1.2, alpha=0.5)
    ax.plot(off, i, "o", ms=8, color=c, mec=c,
            mfc="none" if adv else c, mew=1.4 if adv else 0)
    ax.annotate(fmt.format(off) + ("  advisory" if adv else ""), (off, i),
                textcoords="offset points", xytext=(0, 11), ha="center",
                color=c, fontsize=7.5)


def _tolerance_axis(ax, offsets, xlabel, n):
    ax.axvspan(-0.5, 0.5, color=AQUA, alpha=0.10, lw=0)
    for x in (-0.5, 0.5):
        ax.axvline(x, color=AQUA, lw=0.9, ls=(0, (4, 3)))
    ax.axvline(0, color=MUTED, lw=0.8)
    ax.set_xlabel(xlabel, color=INK2, fontsize=8)
    ax.set_ylim(-0.7, n - 0.3)
    lim = max(1.3, 1.15 * max(abs(o) for o in offsets))
    ax.set_xlim(-lim, lim)


def fig_seasonal(A, res):
    fig = newfig(11.0, 6.6, "Tier 2 shape — the seasonal cycle",
                 "The 12-month climatology of each diagnostic, anchor in grey and "
                 "climate model in blue. Monthly means averaged over each period.",
                 "a model can reproduce the annual mean and still have the annual "
                 "march wrong, which is what this shows and no test of the mean can. "
                 "The amplitude and phase of the annual harmonic are what is actually "
                 "scored; they are in the report tables, not drawn here.")
    gs = fig.add_gridspec(2, 3, hspace=0.52, wspace=0.28)
    x = np.arange(1, 13)
    for i, (key, _, label, unit) in enumerate(SERIES):
        ax = style(fig.add_subplot(gs[i // 3, i % 3]), f"{'abcdef'[i]}   {label}")
        r = res["shape_seasonal"][key]
        ax.plot(x, r["anchor_climatology"], "-o", color=MUTED, lw=1.3, ms=3.0,
                label="anchor")
        ax.plot(x, r["climate_model_climatology"], "-o", color=BLUE, lw=1.6, ms=3.0,
                label="climate model")
        ax.set_xticks(x)
        ax.set_xticklabels(list(MONTH_INITIALS), fontsize=6.5)
        ax.set_ylabel(unit, color=INK2, fontsize=7.5)
        if i == 0:
            lg = ax.legend(fontsize=6.6, frameon=False, loc="best")
            for t in lg.get_texts():
                t.set_color(INK2)
    fig.subplots_adjust(left=0.065, right=0.975, top=fig.text_bottom - 0.055,
                        bottom=0.085)
    footer(fig, res)
    save(fig, "shape_seasonal.png", res["stamp"])


# ------------------------------------------------------- shape: daily distribution
def fig_daily_distribution(A, res, S, segs):
    a = A["daily_DJF_u60N"]
    m = res["tier2_variance"]["daily_DJF_u60N"]
    lo, hi = res["climate_model_period"]

    anchor = np.array(sum((S[s]["_u60n_djf_daily"] for s in segs), []), float)
    # Same rule 17 applies: the daily series carries no year labels, so it can only
    # be subset when the window covers whole segments. 17 has already decided that,
    # so follow its verdict rather than re-deriving it.
    model = None
    if m.get("available", True):
        model = np.array(sum((S[s]["_u60n_djf_daily"] for s in segs
                              if np.asarray(S[s]["_annual_years"], float).min() >= lo
                              and np.asarray(S[s]["_annual_years"], float).max() <= hi),
                             []), float)

    fig = newfig(8.4, 5.0, "Tier 2 shape — the daily distribution",
                 f"Every DJF day of u at 60°N, 10 hPa: {len(anchor)} anchor days in "
                 f"grey, the climate model outlined in blue. Densities, so the two "
                 f"are comparable at different lengths.",
                 "the width is the point. A model that has smoothed away its own "
                 "weather is narrow here while its mean and its seasonal cycle stay "
                 "right. What is scored is the five per-winter percentiles, in the "
                 "report table; this is what they summarise.")
    ax = style(fig.add_subplot(111), None,
               f"σ ratio {m['ratio']:.2f} · window {m['window'][0]:.2f}–"
               f"{m['window'][1]:.2f} · "
               f"{'PASS' if m.get('passes') else 'n/a'}"
               if m.get("available", True) else
               "climate-model days unavailable: window is not segment-aligned")

    span = (anchor.min(), anchor.max()) if model is None else (
        min(anchor.min(), model.min()), max(anchor.max(), model.max()))
    bins = np.linspace(span[0], span[1], 46)
    ax.hist(anchor, bins=bins, density=True, color=RULE, lw=0,
            label=f"anchor CESM {res['anchor_period']}")
    if model is not None:
        ax.hist(model, bins=bins, density=True, histtype="step", color=BLUE,
                lw=1.5, label=f"climate model {lo}–{hi}")
    for v, c, ls in ((a["p05"], MUTED, (0, (3, 3))), (a["p95"], MUTED, (0, (3, 3)))):
        ax.axvline(v, color=c, lw=0.9, ls=ls)
    ax.annotate(f"anchor p5 {a['p05']:.1f}", (a["p05"], 1.0),
                xycoords=("data", "axes fraction"), textcoords="offset points",
                xytext=(3, -9), color=MUTED, fontsize=7)
    ax.annotate(f"p95 {a['p95']:.1f}", (a["p95"], 1.0),
                xycoords=("data", "axes fraction"), textcoords="offset points",
                xytext=(3, -9), color=MUTED, fontsize=7)
    ax.set_xlabel("u 60°N, 10 hPa, daily DJF  (m s⁻¹)", color=INK2, fontsize=8)
    ax.set_ylabel("density", color=INK2, fontsize=8)
    lg = ax.legend(frameon=False, fontsize=7, loc="upper left")
    for t in lg.get_texts():
        t.set_color(INK2)
    fig.subplots_adjust(left=0.105, right=0.965, top=fig.text_bottom - 0.075,
                        bottom=0.135)
    footer(fig, res)
    save(fig, "shape_daily_distribution.png", res["stamp"])


# ------------------------------------------------------------- shape: w* profile
def fig_w_star_profile(A, res):
    fig = newfig(11.0, 5.8, "Tier 2 shape — the tropical w* profile",
                 "10°S–10°N w* at six pressure levels. Left: absolute, which is "
                 "ADVISORY only. Centre: each level divided by the profile's own "
                 "vertical mean, which is what the gate is set on, with the "
                 "tolerance shaded. Right: the offsets.",
                 "the absolute profile cannot carry a gate - the grid error on w* at "
                 "70 hPa is larger than the tier-2 tolerance, so an absolute target "
                 "would fail a correct model on grid choice alone. Normalising "
                 "cancels a multiplicative bias that is uniform in height, and that "
                 "uniformity is assumed, not measured. The profile is non-monotonic, "
                 "so a single 70 hPa value constrains none of this structure.")
    gs = fig.add_gridspec(1, 3, wspace=0.36)
    keys = [k for k in A["w_star_profile"]["levels"]]
    levs = [float(k) for k in keys]
    L = res["shape_w_star_profile"]["levels"]

    ax = style(fig.add_subplot(gs[0, 0]), "absolute — advisory")
    ax.plot([L[k]["absolute_advisory"]["anchor"] for k in keys], levs, "-o",
            color=MUTED, lw=1.4, ms=4, label="anchor")
    ax.plot([L[k]["absolute_advisory"]["climate_model"] for k in keys], levs, "-s",
            color=BLUE, lw=1.4, ms=4, label="climate model")
    ax.set_yscale("log"); ax.invert_yaxis()
    ax.set_yticks(levs); ax.set_yticklabels([f"{v:g}" for v in levs], fontsize=7)
    ax.set_ylabel("hPa", color=INK2, fontsize=8)
    ax.set_xlabel("w*  (mm s⁻¹)", color=INK2, fontsize=8)
    ax.legend(fontsize=6.5, frameon=False, loc="upper left")

    ax2 = style(fig.add_subplot(gs[0, 1]), "normalised — gated")
    an = np.array([L[k]["normalised"]["anchor"] for k in keys])
    tl = np.array([L[k]["normalised"]["tolerance"] for k in keys])
    ax2.fill_betweenx(levs, an - tl, an + tl, color=AQUA, alpha=0.13, lw=0)
    ax2.plot(an, levs, "-o", color=MUTED, lw=1.4, ms=4)
    ax2.plot([L[k]["normalised"]["climate_model"] for k in keys], levs, "-s",
             color=BLUE, lw=1.4, ms=4)
    ax2.set_yscale("log"); ax2.invert_yaxis()
    ax2.set_yticks(levs); ax2.set_yticklabels([f"{v:g}" for v in levs], fontsize=7)
    ax2.set_xlabel("w* ÷ profile mean", color=INK2, fontsize=8)

    ax3 = style(fig.add_subplot(gs[0, 2]), "offsets")
    rk = keys[::-1]
    offs = [L[k]["normalised"]["offset_in_sigma"] for k in rk]
    for i, k in enumerate(rk):
        _offset_row(ax3, i, L[k]["normalised"])
    _tolerance_axis(ax3, offs, "offset from the anchor, in anchor σ", len(rk))
    ax3.set_yticks(range(len(rk)))
    ax3.set_yticklabels([f"{k} hPa" for k in rk], fontsize=7.5, color=INK)

    fig.subplots_adjust(left=0.070, right=0.975, top=fig.text_bottom - 0.035,
                        bottom=0.125)
    footer(fig, res)
    save(fig, "shape_w_star_profile.png", res["stamp"])


def latest_validation(climate_model=None):
    """The newest 17_validation__*.json, optionally restricted to one model."""
    hits = sorted(glob.glob(os.path.join(C.OUTDIR, "17_validation__*.json")),
                  key=os.path.getmtime, reverse=True)
    if climate_model:
        slug = re.sub(r"[^A-Za-z0-9.+-]+", "-", climate_model).strip("-")
        hits = [h for h in hits
                if os.path.basename(h).startswith(f"17_validation__{slug}__")]
    if not hits:
        raise SystemExit("no 17_validation__*.json in output/ - run 17_validate.py first"
                         + (f" for {climate_model}" if climate_model else ""))
    return hits[0]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1].strip())
    ap.add_argument("--climate-model", dest="climate_model", default=None,
                    help="draw the newest result for this model; default the newest of any")
    args = ap.parse_args()

    src = latest_validation(args.climate_model)
    print(f"using {os.path.basename(src)}")
    A = json.load(open(os.path.join(C.OUTDIR, "16_anchors_45yr.json")))
    res = json.load(open(src))
    ps = json.load(open(os.path.join(C.OUTDIR, "07_period_split.json")))
    S, segs = ps["series"], [ps["test"], ps["train"]]
    os.makedirs(RESULTS, exist_ok=True)
    fig_tier1(A, res, S, segs)
    fig_mean(A, res)
    fig_variance(A, res)
    fig_counts(A, res)
    fig_seasonal(A, res)
    fig_daily_distribution(A, res, S, segs)
    fig_w_star_profile(A, res)


if __name__ == "__main__":
    main()
