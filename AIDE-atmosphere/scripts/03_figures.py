"""
03 - Figures of the scorecard.

Reads the transport_scorecard__<set>__<scoring>__<date>.json files of one production
date (the newest, or the one named with --date) and writes the PNGs beside them,
stamped the same way, so the figures cannot drift from the reports:

  per set and scoring
    transport_score_table      group scores and passes per run
    transport_heatmap          every metric x run, % off the reference; absolute
                               metrics as a fraction of their bound, scaled so the
                               bound sits at 25 %; the colour steps at +-25 %;
                               metrics outside the scoring faded
  per set (from the extensive scorecard)
    transport_wstar_profile    tropical (15S-15N) w* profile of each run against WACCM6
    transport_age_profile      tropical (10S-10N) entry-age profile of each run against
                               CLaMS (entry offset removed) and WACCM6 REF-D1

  03_figures.py [--date YYYY-MM-DD]
"""

import argparse
import glob
import json
import os
import re
import sys
import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import strat_eval_common as SE

RESULTS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "validation_results", "transport")

# Colours: one series hue (categorical slot 1) against neutral-ink references; a
# blue - gray - red diverging ramp, lightness monotonic on each arm; status tints
# only together with a symbol.
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#dcdbd7"
MODEL = "#2a78d6"
REF1, REF2 = INK, "#8a8984"
DIVERGING = ["#104281", "#256abf", "#6da7ec", "#b7d3f6", "#f0efec",
             "#f6c2c1", "#ec8584", "#cf4746", "#8f1f1f"]
BINS = [-100, -50, -37.5, -25, -12.5, 12.5, 25, 37.5, 50, 100]
STATUS = {"pass": ("#0ca30c", "✓"), "near": ("#fab219", "~"), "fail": ("#d03b3b", "✗")}
GROUP_ROWS = [("A", "A  age of air"), ("B", "B  residual circulation"), ("C", "C  polar vortex"),
              ("combined", "combined A-C"), ("D", "D  variability")]

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "text.color": INK, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "axes.edgecolor": INK2, "axes.linewidth": 0.6, "font.size": 9, "axes.titlesize": 10,
    "lines.linewidth": 1.5, "grid.color": GRID, "grid.linewidth": 0.5,
})


def stamped(name, stamp):
    root, ext = os.path.splitext(name)
    return f"{root}__{stamp}{ext}"


def window_text(res):
    return f"{res['window'][0]}-{res['window'][1]}"


def footer(fig, res):
    """Source line at the bottom, broken at its separators if the figure is narrow."""
    parts = [f"transport_scorecard__{res['stamp']}.json", "w* interim (AIDE-eval issue #5)",
             "scripts/03_figures.py"]
    width = int(fig.get_figwidth() * 72 / (0.55 * 7))
    lines = [parts[0]]
    for q in parts[1:]:
        if len(lines[-1]) + 3 + len(q) <= width:
            lines[-1] += f" · {q}"
        else:
            lines.append(q)
    fig.text(0.01, 0.005, "\n".join(lines), fontsize=7, color=INK2, ha="left", va="bottom")


def wrap(text, width_in, size=10):
    return textwrap.wrap(text, int(width_in * 72 / (0.55 * size)), break_on_hyphens=False)


def title(fig, text, width_in, y_in=0.12, size=10):
    """Left-aligned title, wrapped to the figure width; returns how far down it
    reaches, in inches from the top."""
    lines = wrap(text, width_in, size)
    h = fig.get_figheight()
    fig.text(0.01, 1 - y_in / h, "\n".join(lines), ha="left", va="top", fontsize=size)
    return y_in + len(lines) * size * 1.35 / 72


def save(fig, stamp, name):
    out = os.path.join(RESULTS, stamped(name, stamp))
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"wrote {out}")


def x_scaled(mid, value):
    """% off the reference; an absolute metric as a fraction of its bound, bound = 25 %."""
    if mid in SE.ABS_BOUND:
        return 100.0 * SE.TOL * value / SE.ABS_BOUND[mid]
    return value


# --------------------------------------------------------------------- score table
def score_table(res):
    runs, S = res["runs"], {s["run"]: s for s in res["scores"]}
    rows = GROUP_ROWS + [("passed", "passed A-C"), ("D_passed", "passed D")]
    W, H = max(7.6, len(runs) + 3.1), 4.3
    fig = plt.figure(figsize=(W, H))
    top = title(fig, f"Transport scorecard {res['name']}, {res['scoring']} scoring: group score = "
                "mean err of the group's scored metrics (err = |model/ref - 1| / 0.25; lower is "
                "better)", W)
    ax = fig.add_axes([0.01, 0.12, 0.98, 1 - 0.12 - (top + 0.05) / H])
    ax.set_xlim(-2.6, W - 2.6 - 0.1)
    ax.set_ylim(len(rows) + 0.2, -1.0)
    ax.axis("off")
    for j, n in enumerate(runs):
        ax.text(j + 0.5, -0.45, n, ha="center", va="center", fontweight="bold")
    for i, (k, lab) in enumerate(rows):
        bold = k == "combined"
        ax.text(-2.5, i + 0.5, lab, ha="left", va="center", fontweight="bold" if bold else "normal")
        for j, n in enumerate(runs):
            s = S[n]
            if k in ("passed", "D_passed"):
                tot = s["total"] if k == "passed" else s["D_total"]
                ax.text(j + 0.5, i + 0.5, f"{s[k]} / {tot}", ha="center", va="center", color=INK2)
                continue
            col, sym = STATUS[SE.status(s[k])]
            ax.add_patch(plt.Rectangle((j + 0.03, i + 0.06), 0.94, 0.88, color=col, alpha=0.22, lw=0))
            ax.text(j + 0.5, i + 0.5, f"{s[k]:.2f} {sym}", ha="center", va="center",
                    fontweight="bold" if bold else "normal")
    fig.text(0.01, 0.06, "✓ pass (err ≤ 1)   ~ near (≤ 1.5)   ✗ fail (> 1.5)   ·   combined = "
             f"mean of A, B, C; D separate   ·   window {window_text(res)}", fontsize=8, color=INK2)
    footer(fig, res)
    save(fig, res["stamp"], "transport_score_table.png")


# --------------------------------------------------------------------- heatmap
def heatmap(res):
    runs, rows = res["runs"], res["rows"]
    X = np.full((len(rows), len(runs)), np.nan)
    labels, texts = [], []
    for i, r in enumerate(rows):
        lab = f"{r['id']}  {r['label']}"
        if r["id"] in SE.ABS_BOUND:
            lab += f"  [abs: ±{SE.ABS_BOUND[r['id']]:g} {r['unit']} = 25 %]"
        labels.append(lab if r["scored"] else f"◦ {lab}  (not scored)")
        row_t = []
        for j, n in enumerate(runs):
            v = res["rows"][i]["runs"][n]["value"]
            X[i, j] = x_scaled(r["id"], v)
            nd = (2 if r["unit"] in ("yr", "K") else 1) if r["id"] in SE.ABS_BOUND else 0
            row_t.append("n/a" if not np.isfinite(v) else f"{v:+.{nd}f}")
        texts.append(row_t)
    cmap = ListedColormap(DIVERGING)
    norm = BoundaryNorm(BINS, cmap.N)
    LBL, CW = 5.3, (0.8 if len(runs) > 3 else 1.3)
    W, H = LBL + CW * len(runs) + 0.3, 12.8
    fig = plt.figure(figsize=(W, H))
    top = title(fig, f"Every metric, every run: model against reference, {window_text(res)}; "
                f"{res['name']}, {res['scoring']} scoring", W)
    ax = fig.add_axes([LBL / W, 1.25 / H, CW * len(runs) / W, 1 - (1.25 + top + 0.35) / H])
    im = ax.imshow(np.clip(X, -99.9, 99.9), cmap=cmap, norm=norm, aspect="auto")
    # 2 px surface gaps between cells, wider between groups
    for j in range(len(runs) + 1):
        ax.axvline(j - 0.5, color=SURFACE, lw=2, clip_on=False)
    groups = [r["group"] for r in rows]
    for i in range(len(rows) + 1):
        new_group = 0 < i < len(rows) and groups[i] != groups[i - 1]
        ax.axhline(i - 0.5, color=SURFACE, lw=7 if new_group else 2, clip_on=False)
    for i, r in enumerate(rows):
        if not r["scored"]:                       # faded: surface veil over the row
            ax.add_patch(plt.Rectangle((-0.5, i - 0.5), len(runs), 1.0, color=SURFACE,
                                       alpha=0.72, lw=0, zorder=2))
        for j in range(len(runs)):
            x = X[i, j]
            dark = r["scored"] and np.isfinite(x) and abs(x) > 37.5
            ax.text(j, i, texts[i][j], ha="center", va="center", fontsize=7.5, zorder=3,
                    color=SURFACE if dark else (INK if r["scored"] else INK2),
                    fontweight="bold" if r["scored"] and np.isfinite(x) and abs(x) > 25 else "normal")
    ax.set_xticks(range(len(runs)), runs)
    ax.xaxis.tick_top()
    ax.set_yticks(range(len(rows)), labels, fontsize=7.5)
    for t, r in zip(ax.get_yticklabels(), rows):
        t.set_color(INK if r["scored"] else INK2)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    cax = fig.add_axes([0.06, 0.62 / H, 0.88, 0.14 / H])
    cb = fig.colorbar(im, cax=cax, orientation="horizontal", extend="both",
                      ticks=[-50, -37.5, -25, -12.5, 12.5, 25, 37.5, 50])
    cb.ax.set_xticklabels(["−50", "−37.5", "−25", "−12.5", "12.5", "25", "37.5", "50"])
    cb.ax.tick_params(labelsize=8, length=0)
    cb.outline.set_visible(False)
    for v in (-25, 25):
        cb.ax.axvline(v, color=INK, lw=2)
    cb.ax.set_title("% off the reference; black lines and bold text: the ±25 % tolerance. Absolute "
                    "metrics [abs]: text in their own unit, colour scaled so their bound = 25 %",
                    fontsize=7.5, color=INK2, loc="left", pad=4) if W > 9 else \
        cb.ax.set_title("% off the reference; black lines and bold text: ±25 % tolerance.\n"
                        "[abs]: text in own unit, colour scaled so the bound = 25 %",
                        fontsize=7.5, color=INK2, loc="left", pad=4)
    footer(fig, res)
    save(fig, res["stamp"], "transport_heatmap.png")


# --------------------------------------------------------------------- profiles
def profile_figure(runs, title_text, n_legend, ncol):
    """Small multiples, one run per panel, under a wrapped title and a legend."""
    nc = min(4, len(runs))
    nr = int(np.ceil(len(runs) / nc))
    W = max(3.0 * nc, 6.0)
    t_in = 0.12 + len(wrap(title_text, W)) * 10 * 1.35 / 72
    header = t_in + 0.12 + int(np.ceil(n_legend / ncol)) * 0.2 + 0.3
    H = 3.7 * nr + header + 0.85
    fig, axes = plt.subplots(nr, nc, figsize=(W, H), sharex=True, sharey=True, squeeze=False)
    title(fig, title_text, W)
    fig.subplots_adjust(left=0.6 / W, right=1 - 0.1 / W, top=1 - header / H, bottom=0.85 / H,
                        hspace=0.25, wspace=0.08)
    return fig, axes.ravel(), nc, (t_in + 0.08) / H


def finish_panels(fig, axes, runs, legend_top, ncol):
    for ax in axes[len(runs):]:
        ax.axis("off")
    h, lab = axes[0].get_legend_handles_labels()
    fig.legend(h, lab, loc="upper left", bbox_to_anchor=(0.04, 1 - legend_top), ncol=ncol,
               frameon=False, fontsize=8)


def wstar_profile(res):
    runs, B = res["runs"], res["references"]["B"]
    pw, ww = np.array(B["profile_waccm_p"]), np.array(B["profile_waccm_wstar_tropics"])
    targets = [100.0, 70.0, 50.0, 30.0, 10.0]
    rows = {r["id"]: r for r in res["rows"]}
    ncol = 4 if len(runs) >= 4 else 1
    text = (f"Tropical residual vertical velocity w*, {window_text(res)} mean, against WACCM6 "
            f"{B['years'][0]}-{B['years'][1]} (same estimator; interpolated linearly in log p to "
            "100/70/50/30/10 hPa)")
    fig, axes, nc, legend_top = profile_figure(runs, text, 4, ncol)
    for ax, n in zip(axes, runs):
        P = res["profiles"][n]
        for lv in targets:
            ax.axhline(lv, color=GRID, lw=0.8, zorder=0)
        ax.plot(ww, pw, color=REF1, label=f"WACCM6 histSST {B['years'][0]}-{B['years'][1]}")
        ax.plot(P["wstar_tropics"], P["wstar_p"], color=MODEL, label=f"run {window_text(res)}")
        m = [rows[f"B1-{lv:g}"]["runs"][n]["model"] for lv in targets]
        r = [rows[f"B1-{lv:g}"]["runs"][n]["ref"] for lv in targets]
        ax.plot(r, targets, "o", ms=6, mfc=SURFACE, mec=REF1, mew=1.5, ls="none",
                label="WACCM6 at the target pressure")
        ax.plot(m, targets, "o", ms=6, color=MODEL, mec=SURFACE, mew=1.0, ls="none",
                label="run at the target pressure (scored)")
        ax.set_title(n, loc="left")
        ax.grid(axis="x")
    ax0 = axes[0]
    ax0.set_yscale("log")
    ax0.set_ylim(200, 5)
    ax0.set_xlim(0, 1.0)
    ax0.set_xticks([0, 0.2, 0.4, 0.6, 0.8])
    ax0.set_yticks([200, 100, 70, 50, 30, 20, 10, 5], ["200", "100", "70", "50", "30", "20", "10", "5"])
    for ax in axes[::nc]:
        ax.set_ylabel("pressure [hPa]")
    for ax in axes[len(runs) - nc:len(runs)]:
        ax.set_xlabel("w* 15S-15N [mm/s]")
    footer(fig, res)
    finish_panels(fig, axes, runs, legend_top, ncol)
    save(fig, res["set_stamp"], "transport_wstar_profile.png")


def age_profile(res):
    runs, A = res["runs"], res["references"]["A"]
    pc, ac = np.array(A["profile_clams_p"]), np.array(A["profile_clams_entry_tropics"])
    pw, aw = np.array(A["profile_waccm_p"]), np.array(A["profile_waccm_entry_tropics"])
    y0, y1 = A["years"][0], A["years"][-1]
    ncol = 3 if len(runs) >= 4 else 1
    text = (f"Tropical entry age of air, {window_text(res)} mean (grey lines: 70 and 10 hPa for "
            "A4, 55 hPa for A1 and A3)")
    fig, axes, nc, legend_top = profile_figure(runs, text, 3, ncol)
    for ax, n in zip(axes, runs):
        P = res["profiles"][n]
        for lv in (70.0, 55.0, 10.0):
            ax.axhline(lv, color=GRID, lw=0.8, zorder=0)
        ax.plot(ac, pc, color=REF1, label=f"CLaMS v3.1 (ERA5) {y0}-{y1}, minus its 150 hPa tropical age")
        ax.plot(aw, pw, color=REF2, ls=(0, (5, 2)), label=f"WACCM6 REF-D1 entry age {y0}-{y1}")
        ax.plot(P["age_tropics"], P["age_p"], color=MODEL, label=f"run, AOA150 {window_text(res)}")
        ax.set_title(n, loc="left")
        ax.grid(axis="x")
    ax0 = axes[0]
    ax0.set_yscale("log")
    ax0.set_ylim(150, 3)
    ax0.set_xlim(0, 5)
    ax0.set_xticks([0, 1, 2, 3, 4])
    ax0.set_yticks([150, 100, 70, 55, 30, 10, 5, 3], ["150", "100", "70", "55", "30", "10", "5", "3"])
    for ax in axes[::nc]:
        ax.set_ylabel("pressure [hPa]")
    for ax in axes[len(runs) - nc:len(runs)]:
        ax.set_xlabel("entry age 10S-10N [yr]")
    footer(fig, res)
    finish_panels(fig, axes, runs, legend_top, ncol)
    save(fig, res["set_stamp"], "transport_age_profile.png")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1].strip())
    ap.add_argument("--date", help="production date of the scorecards to draw (default: the newest)")
    a = ap.parse_args()
    found = sorted(glob.glob(os.path.join(RESULTS, "transport_scorecard__*__*__*.json")))
    dates = sorted({re.search(r"__(\d{4}-\d{2}-\d{2})\.json$", f).group(1) for f in found})
    if not dates:
        raise SystemExit(f"no transport_scorecard__<set>__<scoring>__<date>.json in {RESULTS}; "
                         "run 02_scorecard.py first")
    date = a.date or dates[-1]
    drawn = set()
    for f in [f for f in found if f.endswith(f"__{date}.json")]:
        res = json.load(open(f))
        score_table(res)
        heatmap(res)
        if res["set_stamp"] in drawn:
            continue
        drawn.add(res["set_stamp"])
        wstar_profile(res)
        age_profile(res)


if __name__ == "__main__":
    main()
