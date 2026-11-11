"""
23 - Score the transport diagnostics and write the transport scorecard.

The second stage of the transport scorecard (docs/TRANSPORT_SCORECARD.md). Reads
the per-source output of 22 and scores every metric as jcm-strat
relative_scores.py does:

  rel = model / ref - 1,  err = |rel| / 0.25        pass err <= 1, near <= 1.5, fail above
  A2, C1 latitudes, C2 and D4 final-warming dates keep their absolute bound
  A4 is scored against the midpoint 2.75 yr of its 2.3-3.2 yr range

Group score = mean err of the group's metrics; combined = mean of the A, B and C
group scores; D is reported on its own.

  23_transport_scorecard.py [NAME ...] [--name SET]     default: the first set of
                                                       transport_common.SETS

Every artefact is stamped `__<set>__<production date>`, as 17 stamps its report.

Output: validation_results/transport/transport_scorecard__<stamp>.json
        validation_results/transport/transport_scorecard__<stamp>.md
"""

import argparse
import datetime
import json
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import aide_val_common as C
import transport_common as T

RESULTS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "validation_results", "transport")
SYMBOL = {"pass": "✅", "near": "⚠️", "fail": "❌"}
BASIS = {"A2": "< 0.4 yr", "C2-N": "± 2 K", "C2-S": "± 2 K", "D4-fwN": "± 10 days",
         "D4-fwS": "± 10 days"}
GROUPS = (("A", "A age of air"), ("B", "B residual circulation"), ("C", "C polar vortex"),
          ("combined", "**combined A-C**"), ("D", "D variability"))
FIGURES = [
    ("transport_score_table.png", "group scores and passes per run"),
    ("transport_heatmap.png", "every metric, % off its reference (absolute metrics as a "
     "fraction of their bound), ±25 % marked"),
    ("transport_wstar_profile.png", "tropical w* profiles against WACCM6"),
    ("transport_age_profile.png", "tropical entry-age profiles against CLaMS and WACCM6"),
]


def stamp_for(name, produced):
    """`<set>__<date>`, reduced to characters that are safe in a path (as 17)."""
    return f"{re.sub(r'[^A-Za-z0-9.+-]+', '-', name).strip('-')}__{produced}"


def ref_text(ref):
    if isinstance(ref, list):
        return f"{ref[0]:g}-{ref[1]:g} (mid {T.ref_value(ref):.2f})"
    return f"{ref:.3g}"


def score_all(P, names):
    """Rows (one per metric, every run scored) and per-run group scores."""
    rows, per = [], {n: {g: [] for g in "ABCD"} for n in names}
    for m0 in P[names[0]]["metrics"]:
        mid, absolute = m0["id"], m0["id"] in T.ABS_BOUND
        refs = {n: next(x for x in P[n]["metrics"] if x["id"] == mid)["ref"] for n in names}
        same = len({json.dumps(r) for r in refs.values()}) == 1
        r = dict(id=mid, group=m0["group"], label=m0["label"],
                 unit=T.ABS_UNIT.get(mid, "deg") if absolute else "%",
                 basis=(BASIS.get(mid, "± 4 deg") if absolute else f"± {T.TOL:.0%} of reference"),
                 reference=ref_text(m0["ref"]) + ("" if same else " (first run; per model grid)"),
                 runs={})
        for n in names:
            m = next(x for x in P[n]["metrics"] if x["id"] == mid)
            val, err = T.score(mid, m["model"], m["ref"])
            r["runs"][n] = dict(model=m["model"], ref=m["ref"], value=val, err=err,
                                status=T.status(err))
            per[n][m0["group"]].append(err)
        rows.append(r)
    scores = []
    for n in names:
        g = {k: float(np.mean(v)) for k, v in per[n].items()}
        scores.append(dict(run=n, A=g["A"], B=g["B"], C=g["C"],
                           combined=float(np.mean([g["A"], g["B"], g["C"]])),
                           passed=int(sum(e <= 1 for k in "ABC" for e in per[n][k])),
                           total=int(sum(len(per[n][k]) for k in "ABC")),
                           D=g["D"], D_passed=int(sum(e <= 1 for e in per[n]["D"])),
                           D_total=len(per[n]["D"])))
    return rows, scores


def markdown(res, P, R):
    names = res["runs"]
    first = P[names[0]]
    L = ["# Transport scorecard", "",
         "| | |", "|---|---|",
         f"| Scored | {', '.join(names)} |",
         f"| Read by | reader `{'`, `'.join(sorted({P[n]['reader'] for n in names}))}`; "
         f"{'; '.join(sorted({os.path.dirname(P[n]['path']) for n in names}))} |",
         f"| Window | {first['days'][0]} – {first['days'][1]} ({first['n_days']} days) |",
         f"| Produced | {res['produced']} |",
         f"| Scoring | rel = model / ref − 1, err = \\|rel\\| / {T.TOL:g}; pass err ≤ 1, near ≤ 1.5, "
         "fail above; A2, C1 latitudes, C2 and D4 final warming keep absolute bounds |",
         "| Combined | mean of the A, B and C group scores (each the mean err of its metrics); "
         "D reported separately |",
         "| w* estimator | jcm-strat Ψ* route (`transport_common.tem_streamfunction`, "
         "`residual_w`), linear in log p to the target pressure — **interim**, AIDE-eval issue #5 |",
         "| Generated by | `scripts/23_transport_scorecard.py`, from `output/22_transport__<run>.json` |",
         "", "## Group scores", "",
         "Mean err per group (lower is better; ≤ 1 = on average inside the tolerance).", "",
         "| | " + " | ".join(names) + " |", "|---|" + "---|" * len(names)]
    for k, lab in GROUPS:
        L.append(f"| {lab} | " + " | ".join(
            (f"**{s[k]:.2f}**" if k == "combined" else f"{s[k]:.2f}") for s in res["scores"]) + " |")
    L.append("| passed A-C | " + " | ".join(f"{s['passed']} / {s['total']}" for s in res["scores"]) + " |")
    L.append("| passed D | " + " | ".join(f"{s['D_passed']} / {s['D_total']}" for s in res["scores"]) + " |")
    L += ["", "## Metrics", "",
          "% off the reference; absolute metrics as the difference in their own unit "
          "(jet latitudes as \\|model\\| − \\|reference\\|). ✅ pass, ⚠️ near, ❌ fail.", "",
          "| metric | reference | basis | " + " | ".join(names) + " |",
          "|---|---|---|" + "---|" * len(names)]
    for r in res["rows"]:
        u = "%" if r["unit"] == "%" else ""
        nd = 1 if r["unit"] == "%" else 2
        L.append(f"| **{r['id']}** {r['label']} | {r['reference']} | {r['basis']} | " + " | ".join(
            f"{SYMBOL[c['status']]} {c['value']:+.{nd}f}{u}" if np.isfinite(c["value"])
            else f"{SYMBOL[c['status']]} n/a" for c in (r["runs"][n] for n in names)) + " |")
    A, B, D = R["A"], R["B"], R["D"]
    L += ["", "## References", "",
          f"- **A1** WACCM6 REF-D1 entry age {A['years'][0]}–{A['years'][-1]}, 10S–10N at the level "
          "nearest 55 hPa",
          f"- **A2, A3, A4 profile** CLaMS v3.1 (ERA5) mean age {A['years'][0]}–{A['years'][-1]}; "
          f"A2 compares against CLaMS minus its own 10S–10N age at 150 hPa "
          f"({A['clams_offset_150']:.3f} yr). A4 range 2.3–3.2 yr (CLaMS 3.0, WACCM6 2.5)",
          f"- **B1, B2** CESM2.1.5-WACCM6 histSST h6 {B['years'][0]}–{B['years'][1]} "
          f"({B['n_days']} days, {B['n_levels']} levels from {B['top_hPa']:.2g} hPa), same "
          "estimator as the model",
          f"- **C1, C2** ERA5 monthly zonal means of the same months, on each model's grid",
          f"- **D1–D6** ERA5 CDS daily zonal means {D['years'][0]}–{D['years'][-1]}; **D7** CLaMS "
          "v3.1 monthly age at 55 hPa, same years",
          "", "## Figures", ""]
    L += [f"- `{T_stamped(f, res['stamp'])}` — {d}" for f, d in FIGURES]
    L += ["", "Method, contract and decisions: `docs/TRANSPORT_SCORECARD.md`. Agreement with "
          "jcm-strat's own numbers: `port_check.md` in this folder."]
    return "\n".join(L) + "\n"


def T_stamped(name, stamp):
    root, ext = os.path.splitext(name)
    return f"{root}__{stamp}{ext}"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1].strip())
    first = next(iter(T.SETS))
    ap.add_argument("names", nargs="*", help=f"sources to score (default: the set {first})")
    ap.add_argument("--name", default=first, help="name of the scored set, used in the stamp")
    a = ap.parse_args()
    names = a.names or list(T.SETS[first])
    window = T.window_of(names)
    P = {n: json.load(open(os.path.join(C.OUTDIR, f"22_transport__{n}.json"))) for n in names}
    R = json.load(open(os.path.join(C.OUTDIR, f"22_transport__references__{window[0]}-{window[1]}.json")))
    produced = datetime.date.today().isoformat()
    stamp = stamp_for(a.name, produced)
    rows, scores = score_all(P, names)
    res = dict(name=a.name, produced=produced, stamp=stamp, runs=names, window=list(window),
               tolerance=T.TOL, a4_range=list(T.A4_RANGE), abs_bound=T.ABS_BOUND,
               wstar_estimator="jcm-strat Psi* route, log-p interpolated; interim (AIDE-eval issue #5)",
               rows=rows, scores=scores,
               profiles={n: P[n]["profiles"] for n in names},
               references=R)
    os.makedirs(RESULTS, exist_ok=True)
    out_j = os.path.join(RESULTS, T_stamped("transport_scorecard.json", stamp))
    out_m = os.path.join(RESULTS, T_stamped("transport_scorecard.md", stamp))
    json.dump(res, open(out_j, "w"), indent=1)
    open(out_m, "w").write(markdown(res, P, R))
    for s in scores:
        print(f"{s['run']:8s} A {s['A']:.2f}  B {s['B']:.2f}  C {s['C']:.2f}  combined "
              f"{s['combined']:.2f}  passed {s['passed']}/{s['total']}  D {s['D']:.2f} "
              f"({s['D_passed']}/{s['D_total']})")
    print(f"wrote {out_j}\nwrote {out_m}")


if __name__ == "__main__":
    main()
