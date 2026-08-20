"""
15 - What the protocol looks like when the thing being scored is CESM itself.

Page 1-2  TIER 1. The gate (docs/EVALUATION_PROTOCOL.md) puts every individual year
          of a 5-year rollout inside +/-3 sigma of the 1996-2014 anchor mean. It is
          run here against CESM 1970-1995, treated as five consecutive 5-year
          rollouts.

Page 3    TIER 2. CESM 1970-1994 scored as a rollout against the 35-year 1980-2014
          anchor - mean and variance, the two things tier 2 tests.

This is the out-of-sample discipline script 07 applies to the bias targets: a
threshold CESM's own output cannot meet is not a threshold, it is a trap.

Note on the tier-2 windows: 1980-2014 and 1970-1994 OVERLAP by 15 years, so that
test is not independent the way 07's split is, and 1980-2014 spans the 1995/96
restart between two separate runs. Both are stated in the protocol document.

Outputs: output/15_screen_out_of_sample.json
         stale/AIDE_WACCM_screening_1970-1995.pdf   (archived; see stale/README.md)
"""

import json
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import aide_val_common as C
from report_layout import (Page, figure_page, PAGE_W, PAGE_H, SURFACE, INK,
                           INK2, MUTED, RULE, BLUE, ORANGE)

OUT_J = os.path.join(C.OUTDIR, "15_screen_out_of_sample.json")
# Written beside this script, in stale/, not into docs/ - this figure is archived.
OUT_P = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                     "AIDE_WACCM_screening_1970-1995.pdf")
TRAIN, TEST = "1996-2014", "1970-1995"
K = 3.0                                   # the agreed screening band, in sigma
BLOCKS = [(1970, 1974), (1975, 1979), (1980, 1984), (1985, 1989), (1990, 1994)]
BAND = "#e9e8e3"
BAND1 = "#d6d5cf"
RED = "#c8322a"
GREEN = "#1baf7a"

PANELS = [
    ("mass_flux", "_annual_years", "Mass flux, 70 hPa", "10⁹ kg s⁻¹"),
    ("w_star", "_annual_years", "w̄*, 10°S–10°N", "mm s⁻¹"),
    ("vortex_NH", "_djf_years", "Vortex NH, DJF", "m s⁻¹"),
    ("vortex_SH", "_jja_years", "Vortex SH, JJA", "m s⁻¹"),
    ("polar_cap_T_NH", "_djf_years", "Polar cap T, NH", "K"),
    ("polar_cap_T_SH", "_jja_years", "Polar cap T, SH", "K"),
]

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "axes.edgecolor": RULE, "axes.linewidth": 0.8, "font.size": 8.5,
    "axes.titlesize": 9, "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "axes.labelcolor": INK2,
    "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
})


def load():
    ps = json.load(open(os.path.join(C.OUTDIR, "07_period_split.json")))
    tiers = json.load(open(os.path.join(C.OUTDIR, "14_evaluation_tiers.json")))
    return ps, tiers


def main():
    ps, tiers = load()
    S = ps["series"]

    res = {"band_sigma": K, "anchor": TRAIN, "screened": TEST,
           "blocks": [list(b) for b in BLOCKS], "diagnostics": {}}

    # ------------------------------------------------------------- scoring
    D = {}
    for key, ykey, label, unit in PANELS:
        mu = tiers["tier1"][key]["anchor_mean"]
        sig = tiers["tier1"][key]["sigma"]
        yt = np.array(S[TEST][ykey], float)
        vt = np.array(S[TEST][key], float)
        ya = np.array(S[TRAIN][ykey], float)
        va = np.array(S[TRAIN][key], float)

        z_raw = (vt - mu) / sig
        D[key] = dict(mu=mu, sig=sig, yt=yt, vt=vt, ya=ya, va=va,
                      z_raw=z_raw, label=label, unit=unit)

        blocks = []
        for lo, hi in BLOCKS:
            m = (yt >= lo) & (yt <= hi)
            blocks.append(dict(
                years=[int(a) for a in yt[m]],
                n_years=int(m.sum()),
                fail_raw=[int(a) for a in yt[m & (np.abs(z_raw) > K)]],
                worst_z_raw=float(z_raw[m][np.argmax(np.abs(z_raw[m]))]) if m.any() else None))
        res["diagnostics"][key] = dict(
            units=unit, anchor_mean=mu, sigma=sig, test_mean=float(vt.mean()),
            band=[mu - K * sig, mu + K * sig],
            offset_in_sigma=float((vt.mean() - mu) / sig),
            n_years=int(len(vt)),
            n_fail_raw=int((np.abs(z_raw) > K).sum()),
            fail_years_raw=[int(a) for a in yt[np.abs(z_raw) > K]],
            worst_z_raw=float(z_raw[np.argmax(np.abs(z_raw))]),
            worst_year_raw=int(yt[np.argmax(np.abs(z_raw))]),
            blocks=blocks)

    tot_raw = sum(v["n_fail_raw"] for v in res["diagnostics"].values())
    nblk_raw = sum(1 for v in res["diagnostics"].values()
                   for b in v["blocks"] if b["fail_raw"])
    res["summary"] = dict(year_checks_failed_raw=tot_raw,
                          block_verdicts_failed_raw=nblk_raw,
                          block_verdicts_total=len(BLOCKS) * len(PANELS))

    print(f"TIER 1 - screening CESM {TEST} against a +/-{K:.0f} sigma band "
          f"anchored on {TRAIN}")
    for key, v in res["diagnostics"].items():
        print(f"  {key:16s} offset {v['offset_in_sigma']:+5.2f} sigma   "
              f"{v['n_fail_raw']:2d}/{v['n_years']:2d} years outside   "
              f"worst {v['worst_z_raw']:+5.2f} ({v['worst_year_raw']})")
    print(f"  -> {tot_raw} year-checks fail; "
          f"{nblk_raw}/{len(BLOCKS)*len(PANELS)} block verdicts fail")

    # ------------------------------------------------------------- figures
    pdf = PdfPages(OUT_P)
    foot = (f"AIDE-WACCM validation · the evaluation protocol run against CESM's own "
            f"output · tier 1 band ±{K:.0f}σ on {TRAIN}, tier 2 anchor "
            f"{tiers['tier2']['anchor_period']}")

    # page 1 - the six panels
    fig, axes = plt.subplots(2, 3, figsize=(PAGE_W, PAGE_H))
    fig.subplots_adjust(left=0.065, right=0.955, top=0.735, bottom=0.105,
                        hspace=0.42, wspace=0.26)
    fig.text(0.055, 0.945, "Tier-1 screening, applied to CESM 1970–1995",
             fontsize=14, weight="bold", color=INK, ha="left", va="top")
    fig.add_artist(plt.Line2D([0.055, 0.955], [0.902, 0.902],
                              transform=fig.transFigure, color=RULE, lw=1.0))
    for i, line in enumerate([
            "Every individual year of 1970–1995 (orange) scored against the ±3σ screening band set on 1996–2014 (blue). Dark grey = ±1σ, where",
            "CESM's own years mostly live; light grey = the ±3σ gate. Red rings mark years the gate rejects. Dotted verticals are the 5-year",
            "screening blocks. Only the mass flux fails, and it fails in the 1970s — that is the forced BDC trend (D9), not a broken configuration."]):
        fig.text(0.055, 0.872 - i * 0.026, line, fontsize=8.4, color=MUTED,
                 ha="left", va="top")

    for ax, (key, ykey, label, unit) in zip(axes.ravel(), PANELS):
        d = D[key]
        mu, sig = d["mu"], d["sig"]
        ax.axhspan(mu - K * sig, mu + K * sig, color=BAND, lw=0, zorder=0)
        ax.axhspan(mu - sig, mu + sig, color=BAND1, lw=0, zorder=1)
        ax.axhline(mu, color=INK2, lw=1.0, zorder=2)
        for lo, hi in BLOCKS[1:]:
            ax.axvline(lo - 0.5, color=RULE, lw=0.7, ls=":", zorder=1)
        ax.plot(d["ya"], d["va"], "o", color=BLUE, ms=2.8, mec=SURFACE, mew=0.5,
                alpha=0.85, zorder=4, label="CESM 1996–2014 (anchor)")
        ax.plot(d["yt"], d["vt"], "-o", color=ORANGE, lw=0.8, ms=3.0,
                mec=SURFACE, mew=0.5, alpha=0.9, zorder=5,
                label="CESM 1970–1995 (screened)")
        bad = np.abs(d["z_raw"]) > K
        if bad.any():
            ax.plot(d["yt"][bad], d["vt"][bad], "o", ms=8.5, mfc="none",
                    mec=RED, mew=1.6, zorder=6, label=f"outside ±{K:.0f}σ")
        ax.set_title(f"{label}", loc="left", color=INK, pad=6)
        ax.set_ylabel(unit)
        ax.set_xlim(1968, 2016)
        n = res["diagnostics"][key]["n_fail_raw"]
        ax.text(0.985, 0.045, f"{n} of {res['diagnostics'][key]['n_years']} years outside",
                transform=ax.transAxes, ha="right", va="bottom", fontsize=7.2,
                color=(RED if n else GREEN))
    axes[0, 0].legend(loc="upper left", fontsize=6.4, framealpha=0.9)
    figure_page(pdf, fig, 1, foot)

    # page 2 - block verdict grid
    fig = plt.figure(figsize=(PAGE_W, PAGE_H))
    fig.text(0.055, 0.945, "The same thing as screening-run verdicts",
             fontsize=14, weight="bold", color=INK, ha="left", va="top")
    fig.add_artist(plt.Line2D([0.055, 0.955], [0.902, 0.902],
                              transform=fig.transFigure, color=RULE, lw=1.0))
    fig.text(0.055, 0.872,
             "1970–1995 split into five consecutive 5-year blocks, each treated as one "
             "screening run against the ±3σ band set on 1996–2014.",
             fontsize=8.4, color=MUTED, ha="left", va="top")

    ax = fig.add_axes([0.235, 0.310, 0.545, 0.455])
    ax.set_xlim(0, len(BLOCKS)); ax.set_ylim(0, len(PANELS))
    ax.invert_yaxis(); ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    for j, (lo, hi) in enumerate(BLOCKS):
        ax.text(j + 0.5, -0.13, f"{lo}–{str(hi)[2:]}", ha="center",
                va="bottom", fontsize=8.0, color=MUTED, clip_on=False)
    for i, (key, _yk, label, _u) in enumerate(PANELS):
        ax.text(-0.10, i + 0.5, label, ha="right", va="center",
                fontsize=8.2, color=INK2, clip_on=False)
        for j in range(len(BLOCKS)):
            bad = res["diagnostics"][key]["blocks"][j]["fail_raw"]
            ax.add_patch(plt.Rectangle(
                (j + 0.06, i + 0.10), 0.88, 0.80,
                facecolor=("#f7dedb" if bad else "#e4f4ec"),
                edgecolor=(RED if bad else GREEN), lw=0.9))
            ax.text(j + 0.5, i + 0.50,
                    ("FAIL\n" + " ".join(str(a)[2:] for a in bad)) if bad else "pass",
                    ha="center", va="center", fontsize=7.4,
                    color=(RED if bad else GREEN),
                    weight=("bold" if bad else "normal"))

    U = res["summary"]
    for i, line in enumerate([
            f"{U['block_verdicts_failed_raw']} of {U['block_verdicts_total']} block verdicts "
            f"fail — both mass flux, both before 1980, {U['year_checks_failed_raw']} failing "
            f"year-checks in total. That is the forced BDC trend (D9), not a broken",
            "configuration: the 1970–1995 mass flux sits 1.77σ below the anchor. The other "
            "five diagnostics never flag in 26 years of out-of-sample CESM, which is the",
            "behaviour a screening gate needs — it does not fire on a healthy model."]):
        fig.text(0.055, 0.225 - i * 0.028, line, fontsize=8.8, color=INK2,
                 ha="left", va="top")
    figure_page(pdf, fig, 2, foot)

    # page 3 - tier 2: 1970-1994 scored against the 1980-2014 anchor
    T2 = tiers["tier2"]
    fig = plt.figure(figsize=(PAGE_W, PAGE_H))
    fig.text(0.055, 0.945,
             f"Tier 2 — CESM {T2['test_period']} scored against the "
             f"{T2['anchor_period']} anchor",
             fontsize=14, weight="bold", color=INK, ha="left", va="top")
    fig.add_artist(plt.Line2D([0.055, 0.955], [0.902, 0.902],
                              transform=fig.transFigure, color=RULE, lw=1.0))
    for i, line in enumerate([
            f"The two things tier 2 tests. Left: the rollout MEAN against the ±0.5σ tolerance. "
            f"Right: the rollout σ against the 35-vs-35-year acceptance window.",
            f"σ is detrended within each window where the trend is significant, as in "
            f"07_period_split. The windows overlap by 15 years (1980–1994), so this is a",
            f"consistency check, not an independent out-of-sample test — and {T2['anchor_period']} "
            f"spans the 1995/96 restart between two separate CESM runs."]):
        fig.text(0.055, 0.872 - i * 0.026, line, fontsize=8.4, color=MUTED,
                 ha="left", va="top")

    labels = [p[2] for p in PANELS]
    keys = [p[0] for p in PANELS]

    axm = fig.add_axes([0.235, 0.335, 0.30, 0.42])
    axm.axvspan(-0.5, 0.5, color="#e4f4ec", lw=0, zorder=0)
    axm.axvline(0, color=INK2, lw=1.0, zorder=2)
    for i, k in enumerate(keys):
        t = T2["test"][k]
        ok = t["passes_mean"]
        axm.plot([0, t["difference_in_sigma"]], [i, i], "-",
                 color=(GREEN if ok else RED), lw=1.2, zorder=3)
        axm.plot(t["difference_in_sigma"], i, "o", ms=7,
                 color=(GREEN if ok else RED), mec=SURFACE, mew=0.8, zorder=4)
        axm.text(t["difference_in_sigma"] - 0.09, i + 0.36,
                 f"{t['difference_in_sigma']:+.2f}σ", fontsize=7.4, ha="right",
                 color=(GREEN if ok else RED))
    axm.set_yticks(range(len(keys))); axm.set_yticklabels(labels, fontsize=8.2)
    axm.invert_yaxis()
    axm.set_xlim(-1.75, 0.95)
    axm.set_xlabel("rollout mean − anchor mean, in σ", fontsize=8.2)
    axm.set_title("a   MEAN — tolerance ±0.5σ", loc="left", color=INK, pad=8)

    lo, hi = T2["variance"]["interannual_sigma_ratio_window"]
    axv = fig.add_axes([0.635, 0.335, 0.30, 0.42])
    axv.axvspan(lo, hi, color="#e4f4ec", lw=0, zorder=0)
    axv.axvline(1.0, color=INK2, lw=1.0, zorder=2)
    for i, k in enumerate(keys):
        t = T2["test"][k]
        ok = t["passes_variance"]
        axv.plot([1.0, t["sigma_ratio"]], [i, i], "-",
                 color=(GREEN if ok else RED), lw=1.2, zorder=3)
        axv.plot(t["sigma_ratio"], i, "o", ms=7, color=(GREEN if ok else RED),
                 mec=SURFACE, mew=0.8, zorder=4)
        axv.text(t["sigma_ratio"], i + 0.36, f"{t['sigma_ratio']:.2f}",
                 fontsize=7.4, ha="center", color=(GREEN if ok else RED))
    axv.set_yticks(range(len(keys))); axv.set_yticklabels([])
    axv.invert_yaxis()
    axv.set_xlim(0.60, 1.42)
    axv.set_xlabel("rollout σ ÷ anchor σ", fontsize=8.2)
    axv.set_title(f"b   VARIANCE — window {lo:.2f}–{hi:.2f}", loc="left",
                  color=INK, pad=8)

    n_mean_fail = sum(1 for k in keys if not T2["test"][k]["passes_mean"])
    n_var_fail = sum(1 for k in keys if not T2["test"][k]["passes_variance"])
    ssw = T2["ssw_count"]
    for i, line in enumerate([
            f"{len(keys) - n_mean_fail} of {len(keys)} diagnostics pass the mean test and "
            f"{len(keys) - n_var_fail} of {len(keys)} pass the variance test. The two failures "
            f"are the upwelling pair, and they fail low —",
            f"mass flux at −1.43σ and w̄* at −0.61σ against a ±0.50σ tolerance. That is the "
            f"forced BDC trend again: 1970–1994 predates most of it.",
            f"Variance is the cleaner result. Every diagnostic sits inside the window, so "
            f"CESM's year-to-year spread is stable across the record even where its mean is not.",
            f"Major SSW: {ssw['test_count']} in {T2['test_period']} against an expected "
            f"{ssw['test_expected']:.1f} [{ssw['test_interval'][0]}, {ssw['test_interval'][1]}] "
            f"at the {ssw['rate_per_winter']:.2f}/winter anchor rate — "
            f"{'pass' if ssw['test_passes'] else 'FAIL'}."]):
        fig.text(0.055, 0.245 - i * 0.028, line, fontsize=8.8, color=INK2,
                 ha="left", va="top")
    figure_page(pdf, fig, 3, foot)

    pdf.close()
    with open(OUT_J, "w") as f:
        json.dump(res, f, indent=2)
    print(f"\nwrote {OUT_P}\nwrote {OUT_J}")


if __name__ == "__main__":
    main()
