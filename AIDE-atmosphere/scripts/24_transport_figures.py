"""
24 - Figures of the transport scorecard.

Reads the newest validation_results/transport/transport_scorecard__<stamp>.json
(or the one named with --stamp) and writes four PNGs beside it, stamped the same
way, so the figures cannot drift from the report:

  transport_score_table      group scores and passes per run
  transport_heatmap          every metric x run, % off the reference; absolute
                             metrics as a fraction of their bound, scaled so the
                             bound sits at 25 %; the colour steps at +-25 %
  transport_wstar_profile    tropical (15S-15N) w* profile of each run against WACCM6
  transport_age_profile      tropical (10S-10N) entry-age profile of each run against
                             CLaMS (entry offset removed) and WACCM6 REF-D1

  24_transport_figures.py [--stamp STAMP]
"""

import argparse
import glob
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
import numpy as np
import json

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import transport_common as T

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

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "text.color": INK, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "axes.edgecolor": INK2, "axes.linewidth": 0.6, "font.size": 9, "axes.titlesize": 10,
    "lines.linewidth": 1.5, "grid.color": GRID, "grid.linewidth": 0.5,
})


def stamped(name, stamp):
    root, ext = os.path.splitext(name)
    return f"{root}__{stamp}{ext}"


def footer(fig, res, src):
    fig.text(0.01, 0.005, f"transport scorecard {res['stamp']} · w* interim (AIDE-eval "
             f"issue #5) · scripts/24_transport_figures.py from {os.path.basename(src)[:19]}…json",
             fontsize=7, color=INK2, ha="left", va="bottom")


def save(fig, res, name):
    out = os.path.join(RESULTS, stamped(name, res["stamp"]))
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"wrote {out}")


# --------------------------------------------------------------------- score table
def score_table(res, src):
    runs, S = res["runs"], {s["run"]: s for s in res["scores"]}
    rows = [("A", "A  age of air"), ("B", "B  residual circulation"), ("C", "C  polar vortex"),
            ("combined", "combined A-C"), ("D", "D  variability"),
            ("passed", "passed A-C"), ("D_passed", "passed D")]
    fig, ax = plt.subplots(figsize=(10.5, 3.9))
    ax.set_xlim(-2.6, len(runs))
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
            st = T.status(s[k])
            col, sym = STATUS[st]
            ax.add_patch(plt.Rectangle((j + 0.03, i + 0.06), 0.94, 0.88, color=col, alpha=0.22, lw=0))
            ax.text(j + 0.5, i + 0.5, f"{s[k]:.2f} {sym}", ha="center", va="center",
                    fontweight="bold" if bold else "normal")
    ax.set_title("Transport scorecard: group score = mean err of the group's metrics "
                 "(err = |model/ref - 1| / 0.25; lower is better)", loc="left")
    fig.text(0.01, 0.06, "✓ pass (err ≤ 1)   ~ near (≤ 1.5)   ✗ fail (> 1.5)   ·   combined = "
             "mean of A, B, C; D separate", fontsize=8, color=INK2)
    footer(fig, res, src)
    fig.subplots_adjust(left=0.01, right=0.99, top=0.9, bottom=0.12)
    save(fig, res, "transport_score_table.png")


# --------------------------------------------------------------------- heatmap
def heatmap(res, src):
    runs, rows = res["runs"], res["rows"]
    X = np.full((len(rows), len(runs)), np.nan)
    labels, texts = [], []
    for i, r in enumerate(rows):
        absolute = r["id"] in T.ABS_BOUND
        bound = T.ABS_BOUND.get(r["id"])
        lab = f"{r['id']}  {r['label']}"
        if absolute:
            lab += f"  [abs: ±{bound:g} {r['unit']} = 25 %]"
        labels.append(lab)
        row_t = []
        for j, n in enumerate(runs):
            c = r["runs"][n]
            v = c["value"]
            X[i, j] = (100.0 * T.TOL * v / bound) if absolute else v
            nd = (2 if r["unit"] in ("yr", "K") else 1) if absolute else 0
            row_t.append("n/a" if not np.isfinite(v) else f"{v:+.{nd}f}")
        texts.append(row_t)
    cmap = ListedColormap(DIVERGING)
    cmap.set_over(DIVERGING[-1])
    cmap.set_under(DIVERGING[0])
    norm = BoundaryNorm(BINS, cmap.N)
    fig, ax = plt.subplots(figsize=(11.5, 12.5))
    im = ax.imshow(np.clip(X, -99.9, 99.9), cmap=cmap, norm=norm, aspect="auto")
    # 2 px surface gaps between cells, wider between groups
    for j in range(len(runs) + 1):
        ax.axvline(j - 0.5, color=SURFACE, lw=2, clip_on=False)
    groups = [r["group"] for r in rows]
    for i in range(len(rows) + 1):
        new_group = 0 < i < len(rows) and groups[i] != groups[i - 1]
        ax.axhline(i - 0.5, color=SURFACE, lw=7 if new_group else 2, clip_on=False)
    for i in range(len(rows)):
        for j in range(len(runs)):
            x = X[i, j]
            dark = np.isfinite(x) and abs(x) > 37.5
            ax.text(j, i, texts[i][j], ha="center", va="center", fontsize=7.5,
                    color=SURFACE if dark else INK,
                    fontweight="bold" if np.isfinite(x) and abs(x) > 25 else "normal")
    ax.set_xticks(range(len(runs)), runs)
    ax.xaxis.tick_top()
    ax.set_yticks(range(len(rows)), labels, fontsize=7.5)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    cb = fig.colorbar(im, ax=ax, orientation="horizontal", fraction=0.025, pad=0.02,
                      extend="both", ticks=[-50, -37.5, -25, -12.5, 12.5, 25, 37.5, 50])
    cb.ax.set_xticklabels(["−50", "−37.5", "−25", "−12.5", "12.5", "25", "37.5", "50"])
    cb.ax.tick_params(labelsize=8, length=0)
    cb.outline.set_visible(False)
    for v in (-25, 25):
        cb.ax.axvline(v, color=INK, lw=2)
    cb.set_label("% off the reference; black lines and bold text: the ±25 % tolerance\n"
                 "absolute metrics [abs]: text in their own unit, colour scaled so their bound = 25 %",
                 fontsize=8, color=INK2)
    ax.set_title("Every metric, every run: model against reference, 2006-2010", loc="left", pad=28)
    footer(fig, res, src)
    fig.subplots_adjust(left=0.42, right=0.975, top=0.94, bottom=0.08)
    save(fig, res, "transport_heatmap.png")


# --------------------------------------------------------------------- profiles
def panels(runs):
    nc = 4
    nr = int(np.ceil(len(runs) / nc))
    fig, axes = plt.subplots(nr, nc, figsize=(12, 4.2 * nr), sharex=True, sharey=True)
    return fig, np.atleast_1d(axes).ravel()


def wstar_profile(res, src):
    runs, B = res["runs"], res["references"]["B"]
    pw, ww = np.array(B["profile_waccm_p"]), np.array(B["profile_waccm_wstar_tropics"])
    targets = [100.0, 70.0, 50.0, 30.0, 10.0]
    rows = {r["id"]: r for r in res["rows"]}
    fig, axes = panels(runs)
    for ax, n in zip(axes, runs):
        P = res["profiles"][n]
        for lv in targets:
            ax.axhline(lv, color=GRID, lw=0.8, zorder=0)
        ax.plot(ww, pw, color=REF1, label="WACCM6 histSST 1996-2014")
        ax.plot(P["wstar_tropics"], P["wstar_p"], color=MODEL, label="run")
        m = [rows[f"B1-{lv:g}"]["runs"][n]["model"] for lv in targets]
        r = [rows[f"B1-{lv:g}"]["runs"][n]["ref"] for lv in targets]
        ax.plot(r, targets, "o", ms=6, mfc=SURFACE, mec=REF1, mew=1.5, ls="none",
                label="WACCM6 at the target pressure")
        ax.plot(m, targets, "o", ms=6, color=MODEL, mec=SURFACE, mew=1.0, ls="none",
                label="run at the target pressure (scored)")
        ax.set_title(n, loc="left")
        ax.grid(axis="x")
    for ax in axes[len(runs):]:
        ax.axis("off")
    ax0 = axes[0]
    ax0.set_yscale("log")
    ax0.set_ylim(200, 5)
    ax0.set_xlim(0, 1.0)
    ax0.set_xticks([0, 0.2, 0.4, 0.6, 0.8])
    ax0.set_yticks([200, 100, 70, 50, 30, 20, 10, 5], ["200", "100", "70", "50", "30", "20", "10", "5"])
    for ax in axes[::4]:
        ax.set_ylabel("pressure [hPa]")
    for ax in axes[-4:]:
        ax.set_xlabel("w* 15S-15N [mm/s]")
    h, lab = ax0.get_legend_handles_labels()
    fig.legend(h, lab, loc="upper left", bbox_to_anchor=(0.05, 0.955), ncol=4, frameon=False, fontsize=8)
    fig.suptitle("Tropical residual vertical velocity w*, 2006-2010 mean, against WACCM6 "
                 "(same estimator; interpolated linearly in log p to 100/70/50/30/10 hPa)",
                 x=0.01, y=0.99, ha="left", fontsize=10)
    footer(fig, res, src)
    fig.subplots_adjust(left=0.06, right=0.99, top=0.86, bottom=0.1, hspace=0.25, wspace=0.08)
    save(fig, res, "transport_wstar_profile.png")


def age_profile(res, src):
    runs, A = res["runs"], res["references"]["A"]
    pc, ac = np.array(A["profile_clams_p"]), np.array(A["profile_clams_entry_tropics"])
    pw, aw = np.array(A["profile_waccm_p"]), np.array(A["profile_waccm_entry_tropics"])
    y0, y1 = A["years"][0], A["years"][-1]
    fig, axes = panels(runs)
    for ax, n in zip(axes, runs):
        P = res["profiles"][n]
        for lv in (70.0, 55.0, 10.0):
            ax.axhline(lv, color=GRID, lw=0.8, zorder=0)
        ax.plot(ac, pc, color=REF1, label=f"CLaMS v3.1 (ERA5) {y0}-{y1}, minus its 150 hPa tropical age")
        ax.plot(aw, pw, color=REF2, ls=(0, (5, 2)), label=f"WACCM6 REF-D1 entry age {y0}-{y1}")
        ax.plot(P["age_tropics"], P["age_p"], color=MODEL, label="run, AOA150")
        ax.set_title(n, loc="left")
        ax.grid(axis="x")
    for ax in axes[len(runs):]:
        ax.axis("off")
    ax0 = axes[0]
    ax0.set_yscale("log")
    ax0.set_ylim(150, 3)
    ax0.set_xlim(0, 5)
    ax0.set_xticks([0, 1, 2, 3, 4])
    ax0.set_yticks([150, 100, 70, 55, 30, 10, 5, 3], ["150", "100", "70", "55", "30", "10", "5", "3"])
    for ax in axes[::4]:
        ax.set_ylabel("pressure [hPa]")
    for ax in axes[-4:]:
        ax.set_xlabel("entry age 10S-10N [yr]")
    h, lab = ax0.get_legend_handles_labels()
    fig.legend(h, lab, loc="upper left", bbox_to_anchor=(0.05, 0.955), ncol=3, frameon=False, fontsize=8)
    fig.suptitle("Tropical entry age of air, 2006-2010 mean (grey lines: 70 and 10 hPa for A4, "
                 "55 hPa for A1 and A3)", x=0.01, y=0.99, ha="left", fontsize=10)
    footer(fig, res, src)
    fig.subplots_adjust(left=0.06, right=0.99, top=0.86, bottom=0.1, hspace=0.25, wspace=0.08)
    save(fig, res, "transport_age_profile.png")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1].strip())
    ap.add_argument("--stamp", help="stamp of the scorecard to draw (default: the newest)")
    a = ap.parse_args()
    pat = f"transport_scorecard__{a.stamp}.json" if a.stamp else "transport_scorecard__*.json"
    found = sorted(glob.glob(os.path.join(RESULTS, pat)), key=os.path.getmtime)
    if not found:
        raise SystemExit(f"no {pat} in {RESULTS}; run 23_transport_scorecard.py first")
    src = found[-1]
    res = json.load(open(src))
    score_table(res, src)
    heatmap(res, src)
    wstar_profile(res, src)
    age_profile(res, src)


if __name__ == "__main__":
    main()
