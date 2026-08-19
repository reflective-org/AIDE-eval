"""
18 - The figures behind validation_results/validation_result.md.

One figure per test family, so every verdict in the report has something to look
at rather than a number to trust:

  tier1_screening.png       every candidate year against the +/-3 sigma band
  tier2_mean.png            offset from the anchor against the +/-0.5 sigma tolerance
  tier2_variance.png        sigma ratios against their 95% windows
  counts_and_relations.png  the SSW count, and the two mechanism slopes
  trends.png                anchor and candidate trends, and what this n resolves

Reads output/16_anchors_45yr.json, output/17_validation.json and the series in
output/07_period_split.json. Run after 17.

Output: PNGs in validation_results/
"""

import json
import os
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


def newfig(w, h, title, lead):
    fig = plt.figure(figsize=(w, h), facecolor=SURFACE)
    fig.text(0.045, 0.965, title, color=INK, fontsize=12.5, va="top")
    fig.text(0.045, 0.925, textwrap.fill(lead, int(w * 13.5)), color=INK2,
             fontsize=8, va="top", linespacing=1.5)
    return fig


def footer(fig, res):
    lo, hi = res["candidate_period"]
    fig.text(0.045, 0.022,
             f"candidate {lo}–{hi}  ·  anchor CESM {res['anchor_period']}  ·  "
             f"estimator {res['estimator']}  ·  "
             + ("candidate lies INSIDE the anchor: self-consistency check"
                if res["inside_anchor"] else "candidate independent of the anchor"),
             color=MUTED, fontsize=6.5)


def save(fig, name):
    p = os.path.join(RESULTS, name)
    fig.savefig(p, dpi=200, facecolor=SURFACE)
    plt.close(fig)
    print(f"wrote {p}")


# ---------------------------------------------------------------- tier 1
def fig_tier1(A, res, S, segs):
    fig = newfig(11.0, 7.2, "Tier 1 — individual years against the ±3σ band",
                 "Grey: anchor years outside the candidate. Blue: candidate years. "
                 "Band: anchor mean ±3σ. A year outside it is a flag, not a verdict.")
    lo, hi = res["candidate_period"]
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
    fig.subplots_adjust(left=0.075, right=0.975, top=0.855, bottom=0.075,
                        hspace=0.62, wspace=0.19)
    footer(fig, res)
    save(fig, "tier1_screening.png")


# ---------------------------------------------------------------- tier 2 mean
def fig_mean(A, res):
    fig = newfig(9.2, 5.4, "Tier 2 — the rollout mean",
                 "Offset of the candidate mean from the anchor, in units of the "
                 "anchor σ. Shaded: the ±0.5σ tolerance. Everything is scaled by σ, "
                 "so the six diagnostics share one axis.")
    ax = style(fig.add_subplot(111))
    keys = [k for k, *_ in SERIES][::-1]
    labels = [l for _, _, l, _ in SERIES][::-1]
    ax.axvspan(-0.5, 0.5, color=AQUA, alpha=0.10, lw=0)
    for x in (-0.5, 0.5):
        ax.axvline(x, color=AQUA, lw=0.9, ls=(0, (4, 3)))
    ax.axvline(0, color=MUTED, lw=0.8)
    for i, k in enumerate(keys):
        m = res["tier2_mean"][k]
        c = PASS if m["passes"] else FAIL
        ax.plot([0, m["offset_in_sigma"]], [i, i], color=c, lw=1.2, alpha=0.5)
        ax.plot(m["offset_in_sigma"], i, "o", ms=8, color=c, mec="none")
        ax.annotate(f"{m['offset_in_sigma']:+.2f}σ",
                    (m["offset_in_sigma"], i), textcoords="offset points",
                    xytext=(0, 11), ha="center", color=c, fontsize=7.5)
    ax.set_yticks(range(len(keys)))
    ax.set_yticklabels(labels, fontsize=8, color=INK)
    ax.set_xlabel("offset from the anchor mean, in anchor σ", color=INK2, fontsize=8)
    ax.set_ylim(-0.7, len(keys) - 0.3)
    lim = max(1.3, 1.15 * max(abs(res["tier2_mean"][k]["offset_in_sigma"])
                             for k in keys))
    ax.set_xlim(-lim, lim)
    ax.text(0.5, 1.0, " tolerance ±0.5σ", transform=ax.get_xaxis_transform(),
            color=AQUA, fontsize=7, va="top")
    fig.subplots_adjust(left=0.20, right=0.965, top=0.80, bottom=0.135)
    footer(fig, res)
    save(fig, "tier2_mean.png")


# ------------------------------------------------------------ tier 2 variance
def fig_variance(A, res):
    fig = newfig(9.2, 5.4, "Tier 2 — the variance",
                 "Ratio of candidate σ to anchor σ, against the 95% window that "
                 "sampling alone allows. The daily DJF ratio is the sharper test: "
                 "≈7 effective samples per winter instead of one.")
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
    ax.set_xlabel("candidate σ ÷ anchor σ", color=INK2, fontsize=8)
    ax.set_ylim(-0.7, len(rows) - 0.3)
    ax.text(0.02, 0.02, "grey bar: 95% window from sampling error alone",
            transform=ax.transAxes, color=MUTED, fontsize=7)
    fig.subplots_adjust(left=0.20, right=0.965, top=0.80, bottom=0.135)
    footer(fig, res)
    save(fig, "tier2_variance.png")


# ------------------------------------------------------- counts and relations
def fig_counts(A, res, S, segs):
    fig = newfig(11.0, 4.3, "Tier 2 — the count, and the mechanism relations",
                 "Left: major NH sudden stratospheric warmings against the Poisson "
                 "interval implied by the anchor rate. Right: the two mechanism "
                 "slopes, candidate fit against the anchor's bootstrap CI.")
    s = res["counts"]["ssw_NH"]
    ax = style(fig.add_subplot(1, 3, 1), "a   Major NH SSW",
               f"{s['candidate_count']} in {s['candidate_winters']} winters · "
               f"expect {s['expected']:.1f} · "
               f"{'PASS' if s['passes'] else 'FAIL'}")
    ilo, ihi = s["interval"]
    xs = np.arange(max(0, ilo - 4), ihi + 5)
    from scipy.stats import poisson
    ax.bar(xs, poisson.pmf(xs, s["expected"]), width=0.85,
           color=[BLUE if ilo <= x <= ihi else RULE for x in xs], lw=0)
    c = PASS if s["passes"] else FAIL
    ax.axvline(s["candidate_count"], color=c, lw=1.6)
    ax.annotate(f"candidate {s['candidate_count']}", (s["candidate_count"], 1.0),
                xycoords=("data", "axes fraction"), textcoords="offset points",
                xytext=(4, -8), color=c, fontsize=7.5)
    ax.set_xlabel("events", color=INK2, fontsize=8)
    ax.set_ylabel("Poisson probability", color=INK2, fontsize=8)

    lo, hi = res["candidate_period"]
    panels = [("R1 wave->vortex", "heat_flux_100", "vortex_NH",
               "v'θ' 100 hPa, 45–75°N  [K m s⁻¹]", "u 60°N DJF  [m s⁻¹]"),
              ("R2 thermal wind", "polar_cap_T_NH", "vortex_NH",
               "polar cap T, NH  [K]", "u 60°N DJF  [m s⁻¹]")]
    for j, (tag, xk, yk, xlab, ylab) in enumerate(panels):
        m = res["mechanism"][tag]
        ax = style(fig.add_subplot(1, 3, j + 2))
        yx, vx = C.join_segments(S, segs, xk, "_djf_years")
        yy, vy = C.join_segments(S, segs, yk, "_djf_years")
        n = min(len(vx), len(vy))
        inc = (yx[:n] >= lo) & (yx[:n] <= hi)
        ax.plot(vx[:n][~inc], vy[:n][~inc], "o", ms=3, color=RULE, mec="none")
        ax.plot(vx[:n][inc], vy[:n][inc], "o", ms=3.6, color=BLUE, mec="none")
        xx = np.linspace(vx[:n].min(), vx[:n].max(), 2)
        a_int = vy[:n].mean() - m["anchor_slope"] * vx[:n].mean()
        ax.plot(xx, m["anchor_slope"] * xx + a_int, color=MUTED, lw=1.1,
                label=f"anchor {m['anchor_slope']:+.2f}")
        c_int = (vy[:n][inc].mean() - m["candidate_slope"] * vx[:n][inc].mean())
        col = PASS if m["passes"] else FAIL
        ax.plot(xx, m["candidate_slope"] * xx + c_int, color=col, lw=1.4,
                label=f"candidate {m['candidate_slope']:+.2f}")
        lg = ax.legend(frameon=False, fontsize=7, loc="best")
        for t in lg.get_texts():
            t.set_color(INK2)
        ax.set_xlabel(xlab, color=INK2, fontsize=7.5)
        ax.set_ylabel(ylab, color=INK2, fontsize=7.5)
        style(ax, f"{'bc'[j]}   {tag}",
              f"anchor CI {m['anchor_ci95'][0]:+.2f} to "
              f"{m['anchor_ci95'][1]:+.2f} · "
              f"{'PASS' if m['passes'] else 'FAIL'}")
    fig.subplots_adjust(left=0.06, right=0.985, top=0.75, bottom=0.19, wspace=0.30)
    footer(fig, res)
    save(fig, "counts_and_relations.png")


# ---------------------------------------------------------------------- trends
def fig_trends(A, res):
    n = res["tier2_mean"]["mass_flux"]["n"]
    fig = newfig(9.2, 5.4, f"Trends, and what n = {n} resolves",
                 "Anchor and candidate trend per decade, each scaled by the 1.96σ "
                 "the candidate length can resolve. Inside the shaded band a trend "
                 "cannot be told from zero, so neither a match nor a mismatch means "
                 "anything.")
    ax = style(fig.add_subplot(111))
    keys = [k for k, *_ in SERIES][::-1]
    labels = [l for _, _, l, _ in SERIES][::-1]
    ax.axvspan(-1, 1, color=RULE, alpha=0.5, lw=0)
    for x in (-1, 1):
        ax.axvline(x, color=MUTED, lw=0.7, ls=(0, (4, 3)))
    ax.axvline(0, color=MUTED, lw=0.8)
    for i, k in enumerate(keys):
        t = res["trend"][k]
        h = 1.96 * t["se_at_candidate_n"]
        ax.axhline(i, color=RULE, lw=0.5, zorder=0)
        ax.plot(t["anchor"] / h, i, "s", ms=6.5, color=INK2, mec="none", zorder=3)
        ax.plot(t["candidate"] / h, i, "o", ms=7.5,
                mfc=BLUE if t["resolvable"] else "none",
                mec=BLUE, mew=1.3, zorder=3)
        ax.annotate("resolvable" if t["resolvable"] else "not resolvable",
                    (1.01, i), xycoords=("axes fraction", "data"),
                    va="center", fontsize=6.5,
                    color=INK2 if t["resolvable"] else MUTED)
    ax.plot([], [], "s", ms=6.5, color=INK2, label="anchor trend")
    ax.plot([], [], "o", ms=7.5, mfc=BLUE, mec=BLUE, label="candidate trend")
    lg = ax.legend(frameon=False, fontsize=7.5, loc="lower left",
                   bbox_to_anchor=(0.0, 1.01), ncol=2, handletextpad=0.4,
                   columnspacing=1.6)
    for t_ in lg.get_texts():
        t_.set_color(INK2)
    ax.set_yticks(range(len(keys)))
    ax.set_yticklabels(labels, fontsize=8, color=INK)
    ax.set_xlabel(f"trend per decade, in units of the 1.96σ resolvable at n = {n}",
                  color=INK2, fontsize=8)
    ax.set_ylim(-0.7, len(keys) - 0.3)
    zs = [res["trend"][k]["anchor"] / (1.96 * res["trend"][k]["se_at_candidate_n"])
          for k in keys] + \
         [res["trend"][k]["candidate"] / (1.96 * res["trend"][k]["se_at_candidate_n"])
          for k in keys]
    lim = 1.25 * max(1.2, max(abs(z) for z in zs))
    ax.set_xlim(-lim, lim)
    ax.text(-0.5 * lim, -0.52, "shaded: cannot be told from zero at this n",
            color=MUTED, fontsize=7, ha="center")
    fig.subplots_adjust(left=0.20, right=0.845, top=0.775, bottom=0.135)
    footer(fig, res)
    save(fig, "trends.png")


def main():
    A = json.load(open(os.path.join(C.OUTDIR, "16_anchors_45yr.json")))
    res = json.load(open(os.path.join(C.OUTDIR, "17_validation.json")))
    ps = json.load(open(os.path.join(C.OUTDIR, "07_period_split.json")))
    S, segs = ps["series"], [ps["test"], ps["train"]]
    os.makedirs(RESULTS, exist_ok=True)
    fig_tier1(A, res, S, segs)
    fig_mean(A, res)
    fig_variance(A, res)
    fig_counts(A, res, S, segs)
    fig_trends(A, res)


if __name__ == "__main__":
    main()
