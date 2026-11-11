"""
23 - Score the transport diagnostics and write the transport scorecard.

The second stage of the transport scorecard (docs/TRANSPORT_SCORECARD.md). Reads
the per-source output of 22 and scores every metric as jcm-strat
relative_scores.py does:

  rel = model / ref - 1,  err = |rel| / 0.25        pass err <= 1, near <= 1.5, fail above
  A2, A2-WACCM, C1 latitudes, C2 and D4 final-warming dates keep their absolute bound
  A4 is scored against the midpoint 2.75 yr of its 2.3-3.2 yr range

Group score = mean err of the group's scored metrics; combined = mean of the A, B
and C group scores; D is reported on its own.

Every result set (transport_common.SETS) is scored twice (transport_common.SCORINGS):
  extensive  the 34 metrics of the phase-29 scorecard
  reduced    jcm-strat phase 31's scored set: A2-WACCM for A2, 14 metrics shown but
             not scored
A set listed in transport_common.COMPARE also gets a run-by-run comparison with
one run of another set (AL7h 2011-2019 against AL7h 2006-2010), both scorings.

  23_transport_scorecard.py [SET ...] [--scoring extensive|reduced]
                                         default: every set, both scorings

Every artefact is stamped `__<set>__<scoring>__<production date>`, as 17 stamps its
report.

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
BASIS = {"A2": "< 0.4 yr", "A2-WACCM": "< 0.4 yr", "C2-N": "± 2 K", "C2-S": "± 2 K",
         "D4-fwN": "± 10 days", "D4-fwS": "± 10 days"}
GROUPS = (("A", "A age of air"), ("B", "B residual circulation"), ("C", "C polar vortex"),
          ("combined", "**combined A-C**"), ("D", "D variability"))
SCORING_TEXT = {
    "extensive": "**extensive**: the 34 metrics of the phase-29 scorecard, every one scored",
    "reduced": "**reduced**: jcm-strat phase 31's scored set — A2-WACCM replaces A2; A2 and "
               f"{', '.join(T.PHASE31_UNSCORED)} are listed but not scored. Phase 31 scored with "
               "absolute bounds; here the relative scoring applies to the reduced set",
}
FIGURES = [
    ("transport_score_table.png", "scoring", "group scores and passes per run"),
    ("transport_heatmap.png", "scoring", "every metric, % off its reference (absolute metrics "
     "as a fraction of their bound), ±25 % marked; unscored metrics faded"),
    ("transport_wstar_profile.png", "set", "tropical w* profiles against WACCM6"),
    ("transport_age_profile.png", "set", "tropical entry-age profiles against CLaMS and WACCM6"),
]
COMPARE_FIGURE = ("transport_window_comparison.png", "set",
                  "AL7h metric by metric, in sample against out of sample, both scorings")


def stamp_for(*parts):
    """`<set>[__<scoring>]__<date>`, each part reduced to characters that are safe in a
    path (as 17)."""
    return "__".join(re.sub(r"[^A-Za-z0-9.+-]+", "-", p).strip("-") for p in parts)


def stamped(name, stamp):
    root, ext = os.path.splitext(name)
    return f"{root}__{stamp}{ext}"


def ref_text(ref):
    if isinstance(ref, list):
        return f"{ref[0]:g}-{ref[1]:g} (mid {T.ref_value(ref):.2f})"
    return f"{ref:.3g}"


def basis(mid):
    if mid in T.ABS_BOUND:
        return BASIS.get(mid, "± 4 deg")
    return f"± {T.TOL:.0%} of reference"


def score_all(P, names, scoring):
    """Rows (one per shown metric, every run scored, `scored` False for a metric
    outside the scoring) and per-run group scores over the scored metrics."""
    sc = T.SCORINGS[scoring]
    rows, per = [], {n: {g: [] for g in "ABCD"} for n in names}
    for m0 in P[names[0]]["metrics"]:
        mid = m0["id"]
        if mid in sc["hidden"]:
            continue
        absolute, scored = mid in T.ABS_BOUND, mid not in sc["unscored"]
        refs = {n: next(x for x in P[n]["metrics"] if x["id"] == mid)["ref"] for n in names}
        same = len({json.dumps(r) for r in refs.values()}) == 1
        r = dict(id=mid, group=m0["group"], label=m0["label"], scored=scored,
                 unit=T.ABS_UNIT.get(mid, "deg") if absolute else "%", basis=basis(mid),
                 reference=ref_text(m0["ref"]) + ("" if same else " (first run; per model grid)"),
                 runs={})
        for n in names:
            m = next(x for x in P[n]["metrics"] if x["id"] == mid)
            val, err = T.score(mid, m["model"], m["ref"])
            r["runs"][n] = dict(model=m["model"], ref=m["ref"], value=val, err=err,
                                status=T.status(err))
            if scored:
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


def comparison(set_name, P):
    """One run of this set against one run of another set, every metric of either
    scoring, and the group scores of both runs under both scorings."""
    cfg = T.COMPARE[set_name]
    other_set, other = cfg["against"]
    pair = {"in": other, "out": cfg["run"]}
    Q = {k: P[n] for k, n in pair.items()}
    metrics = []
    for m0 in Q["out"]["metrics"]:
        mid = m0["id"]
        e = dict(id=mid, group=m0["group"], label=m0["label"],
                 unit=T.ABS_UNIT.get(mid, "deg") if mid in T.ABS_BOUND else "%",
                 scored={s: mid not in d["hidden"] and mid not in d["unscored"]
                         for s, d in T.SCORINGS.items()})
        for k in pair:
            m = next(x for x in Q[k]["metrics"] if x["id"] == mid)
            val, err = T.score(mid, m["model"], m["ref"])
            e[k] = dict(model=m["model"], ref=m["ref"], value=val, err=err, status=T.status(err))
        metrics.append(e)
    groups = {s: {k: score_all({pair[k]: Q[k]}, [pair[k]], s)[1][0] for k in pair}
              for s in T.SCORINGS}
    return dict(run=cfg["run"], against=dict(set=other_set, run=other), labels=list(cfg["labels"]),
                windows={k: Q[k]["window"] for k in pair}, metrics=metrics, groups=groups)


def fmt_value(c, unit):
    if not np.isfinite(c["value"]):
        return f"{SYMBOL[c['status']]} n/a"
    return (f"{SYMBOL[c['status']]} {c['value']:+.1f}%" if unit == "%"
            else f"{SYMBOL[c['status']]} {c['value']:+.2f}")


def metric_table(rows, names, scored):
    L = ["| metric | reference | basis | " + " | ".join(names) + " |",
         "|---|---|---|" + "---|" * len(names)]
    for r in rows:
        if r["scored"] != scored:
            continue
        L.append(f"| **{r['id']}** {r['label']} | {r['reference']} | {r['basis']} | " + " | ".join(
            fmt_value(r["runs"][n], r["unit"]) for n in names) + " |")
    return L


def comparison_md(cmp):
    lab_in, lab_out = cmp["labels"]
    w_in, w_out = (f"{w[0]}-{w[1]}" for w in (cmp["windows"]["in"], cmp["windows"]["out"]))
    run_in, run_out = cmp["against"]["run"], cmp["run"]
    L = ["", f"## {run_in} {w_in} ({lab_in}) against {run_out} {w_out} ({lab_out})", "",
         f"`{run_in}` is scored in the set `{cmp['against']['set']}`; references follow each "
         "window (B: WACCM6 1996–2014 for both). Group scores under both scorings:", "",
         f"| | extensive {w_in} | extensive {w_out} | reduced {w_in} | reduced {w_out} |",
         "|---|---|---|---|---|"]
    G = cmp["groups"]
    for k, lab in GROUPS:
        cells = [G[s][w][k] for s in ("extensive", "reduced") for w in ("in", "out")]
        L.append(f"| {lab} | " + " | ".join(
            (f"**{c:.2f}**" if k == "combined" else f"{c:.2f}") for c in cells) + " |")
    for k, tot, lab in (("passed", "total", "passed A-C"), ("D_passed", "D_total", "passed D")):
        L.append(f"| {lab} | " + " | ".join(f"{G[s][w][k]} / {G[s][w][tot]}"
                                            for s in ("extensive", "reduced") for w in ("in", "out")) + " |")
    L += ["", "Metric by metric: % off the reference (absolute metrics in their own unit). "
          "Scored: ✓ in the scoring's group scores, ◦ listed but not scored, – not part of it. "
          "Change: out of sample − in sample, in % points (or the unit).", "",
          f"| metric | extensive | reduced | reference {w_in} | {w_in} | reference {w_out} | {w_out} | change |",
          "|---|---|---|---|---|---|---|---|"]
    for e in cmp["metrics"]:
        mark = {s: ("✓" if e["scored"][s] else ("–" if e["id"] in T.SCORINGS[s]["hidden"] else "◦"))
                for s in T.SCORINGS}
        d = e["out"]["value"] - e["in"]["value"]
        nd = 1 if e["unit"] == "%" else 2
        L.append(f"| **{e['id']}** {e['label']} | {mark['extensive']} | {mark['reduced']} | "
                 f"{ref_text(e['in']['ref'])} | {fmt_value(e['in'], e['unit'])} | "
                 f"{ref_text(e['out']['ref'])} | {fmt_value(e['out'], e['unit'])} | "
                 f"{d:+.{nd}f} |")
    return L


def markdown(res, P, R):
    names = res["runs"]
    first = P[names[0]]
    A, B, D = R["A"], R["B"], R["D"]
    L = [f"# Transport scorecard: {res['name']}, {res['scoring']} scoring", "",
         "| | |", "|---|---|",
         f"| Scored | {', '.join(names)} |",
         f"| Read by | reader `{'`, `'.join(sorted({P[n]['reader'] for n in names}))}`; "
         f"{'; '.join(sorted({os.path.dirname(P[n]['path']) for n in names}))} |",
         f"| Window | {first['days'][0]} – {first['days'][1]} ({first['n_days']} days) |",
         f"| Reference years | A, D, D7: {A['years'][0]}–{A['years'][-1]} (the window); "
         f"B: WACCM6 h6 {B['years'][0]}–{B['years'][1]}, the only h6 segment |",
         f"| Produced | {res['produced']} |",
         f"| Scoring | {SCORING_TEXT[res['scoring']]} |",
         f"| Score | rel = model / ref − 1, err = \\|rel\\| / {T.TOL:g}; pass err ≤ 1, near ≤ 1.5, "
         "fail above; A2, A2-WACCM, C1 latitudes, C2 and D4 final warming keep absolute bounds |",
         "| Combined | mean of the A, B and C group scores (each the mean err of its scored "
         "metrics); D reported separately |",
         "| w* estimator | jcm-strat Ψ* route (`transport_common.tem_streamfunction`, "
         "`residual_w`), linear in log p to the target pressure — **interim**, AIDE-eval issue #5 |",
         "| Generated by | `scripts/23_transport_scorecard.py`, from `output/22_transport__<run>.json` |",
         "", "## Group scores", "",
         "Mean err over the group's scored metrics (lower is better; ≤ 1 = on average inside "
         "the tolerance).", "",
         "| | " + " | ".join(names) + " |", "|---|" + "---|" * len(names)]
    for k, lab in GROUPS:
        L.append(f"| {lab} | " + " | ".join(
            (f"**{s[k]:.2f}**" if k == "combined" else f"{s[k]:.2f}") for s in res["scores"]) + " |")
    L.append("| passed A-C | " + " | ".join(f"{s['passed']} / {s['total']}" for s in res["scores"]) + " |")
    L.append("| passed D | " + " | ".join(f"{s['D_passed']} / {s['D_total']}" for s in res["scores"]) + " |")
    L += ["", "## Metrics", "",
          "% off the reference; absolute metrics as the difference in their own unit "
          "(jet latitudes as \\|model\\| − \\|reference\\|). ✅ pass, ⚠️ near, ❌ fail.", ""]
    L += metric_table(res["rows"], names, True)
    if any(not r["scored"] for r in res["rows"]):
        L += ["", "### Listed, not scored", "",
              "Outside this scoring's group scores; shown for context with the status they would have.", ""]
        L += metric_table(res["rows"], names, False)
    L += ["", "## References", "",
          f"- **A1** WACCM6 REF-D1 entry age {A['years'][0]}–{A['years'][-1]}, 10S–10N at the level "
          "nearest 55 hPa",
          f"- **A2-WACCM** WACCM6 REF-D1 2-D entry age {A['years'][0]}–{A['years'][-1]}, on its own "
          "grid, 100–5 hPa, \\|lat\\| ≤ 80°, no offset",
          f"- **A2, A3, A4 profile** CLaMS v3.1 (ERA5) mean age {A['years'][0]}–{A['years'][-1]}; "
          f"A2 compares against CLaMS minus its own 10S–10N age at 150 hPa "
          f"({A['clams_offset_150']:.3f} yr). A4 range 2.3–3.2 yr (CLaMS 3.0, WACCM6 2.5)",
          f"- **B1, B2** CESM2.1.5-WACCM6 histSST h6 {B['years'][0]}–{B['years'][1]} "
          f"({B['n_days']} days, {B['n_levels']} levels from {B['top_hPa']:.2g} hPa), same "
          "estimator as the model; the only h6 segment, so not period-matched",
          "- **C1, C2** ERA5 monthly zonal means of the same months, on each model's grid",
          f"- **D1–D6** ERA5 CDS daily zonal means {D['years'][0]}–{D['years'][-1]}; **D7** CLaMS "
          "v3.1 monthly age at 55 hPa, same years"]
    if res.get("comparison"):
        L += comparison_md(res["comparison"])
    L += ["", "## Figures", ""]
    figs = FIGURES + ([COMPARE_FIGURE] if res.get("comparison") else [])
    L += [f"- `{stamped(f, res['stamp'] if kind == 'scoring' else res['set_stamp'])}` — {d}"
          for f, kind, d in figs]
    L += ["", "Method, contract and decisions: `docs/TRANSPORT_SCORECARD.md`. Agreement with "
          "jcm-strat's own numbers: `port_check.md` in this folder."]
    return "\n".join(L) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1].strip())
    ap.add_argument("sets", nargs="*", help="result sets to score (default: every set)")
    ap.add_argument("--scoring", choices=list(T.SCORINGS), help="one scoring (default: both)")
    a = ap.parse_args()
    sets = a.sets or list(T.SETS)
    unknown = [s for s in sets if s not in T.SETS]
    if unknown:
        raise SystemExit(f"unknown set(s) {unknown}; registered: {list(T.SETS)}")
    scorings = [a.scoring] if a.scoring else list(T.SCORINGS)
    produced = datetime.date.today().isoformat()
    os.makedirs(RESULTS, exist_ok=True)
    for set_name in sets:
        names = list(T.SETS[set_name])
        window = T.window_of(names)
        load = [*names] + ([T.COMPARE[set_name]["against"][1]] if set_name in T.COMPARE else [])
        P = {n: json.load(open(os.path.join(C.OUTDIR, f"22_transport__{n}.json"))) for n in load}
        R = json.load(open(os.path.join(C.OUTDIR,
                                        f"22_transport__references__{window[0]}-{window[1]}.json")))
        cmp = comparison(set_name, P) if set_name in T.COMPARE else None
        for scoring in scorings:
            rows, scores = score_all(P, names, scoring)
            res = dict(name=set_name, scoring=scoring, produced=produced,
                       stamp=stamp_for(set_name, scoring, produced),
                       set_stamp=stamp_for(set_name, produced), runs=names, window=list(window),
                       tolerance=T.TOL, a4_range=list(T.A4_RANGE), abs_bound=T.ABS_BOUND,
                       scoring_definition={k: list(v) for k, v in T.SCORINGS[scoring].items()},
                       wstar_estimator="jcm-strat Psi* route, log-p interpolated; interim "
                                       "(AIDE-eval issue #5)",
                       rows=rows, scores=scores, comparison=cmp,
                       profiles={n: P[n]["profiles"] for n in names},
                       references=R)
            out_j = os.path.join(RESULTS, stamped("transport_scorecard.json", res["stamp"]))
            out_m = os.path.join(RESULTS, stamped("transport_scorecard.md", res["stamp"]))
            json.dump(res, open(out_j, "w"), indent=1)
            open(out_m, "w").write(markdown(res, P, R))
            print(f"{set_name}, {scoring}")
            for s in scores:
                print(f"  {s['run']:9s} A {s['A']:.2f}  B {s['B']:.2f}  C {s['C']:.2f}  combined "
                      f"{s['combined']:.2f}  passed {s['passed']}/{s['total']}  D {s['D']:.2f} "
                      f"({s['D_passed']}/{s['D_total']})")
            print(f"wrote {out_j}\nwrote {out_m}")


if __name__ == "__main__":
    main()
