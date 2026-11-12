"""
02 - Score the diagnostics and write the scorecard reports.

The second stage of the transport scorecard (docs/TRANSPORT_SCORECARD.md). Reads
the per-source output of 01 and scores every metric as jcm-strat
relative_scores.py does:

  rel = model / ref - 1,  err = |rel| / 0.25        pass err <= 1, near <= 1.5, fail above
  A2, A2-WACCM, C1 latitudes, C2 and D4 final-warming dates keep their absolute bound
  A4 is scored against the midpoint 2.75 yr of its 2.3-3.2 yr range

Group score = mean err of the group's scored metrics; combined = mean of the A, B
and C group scores; D is reported on its own.

Every result set (strat_eval_common.SETS) is scored twice (strat_eval_common.SCORINGS):
  extensive  the 34 metrics of the phase-29 scorecard
  reduced    jcm-strat phase 31's scored set: A2-WACCM for A2, 14 metrics shown but
             not scored
A set listed in strat_eval_common.COMPARE also gets a run-by-run comparison with
one run of another set (AL7h 2011-2019 against AL7h 2006-2010), both scorings.

  02_scorecard.py [SET ...] [--scoring extensive|reduced]
                                         default: every set, both scorings

Every artefact is stamped `__<set>__<scoring>__<production date>`, as 17 stamps its
report.

Output: validation_results/transport/transport_scorecard__<stamp>.json
        validation_results/transport/transport_scorecard__<stamp>.md
        validation_results/transport/SCORING_EXTENSIVE.md, SCORING_REDUCED.md
            one-page summary per scoring over every set (written when every set is scored)
"""

import argparse
import datetime
import json
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import strat_eval_common as SE

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
               f"{', '.join(SE.PHASE31_UNSCORED)} are listed but not scored. Phase 31 scored with "
               "absolute bounds; here the relative scoring applies to the reduced set",
}
FIGURES = [
    ("transport_score_table.png", "scoring", "group scores and passes per run"),
    ("transport_heatmap.png", "scoring", "every metric, % off its reference (absolute metrics "
     "as a fraction of their bound), ±25 % marked; unscored metrics faded"),
    ("transport_wstar_profile.png", "set", "tropical w* profiles against WACCM6"),
    ("transport_age_profile.png", "set", "tropical entry-age profiles against CLaMS and WACCM6"),
]
SUMMARY = {"extensive": "SCORING_EXTENSIVE.md", "reduced": "SCORING_REDUCED.md"}
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
        return f"{ref[0]:g}-{ref[1]:g} (mid {SE.ref_value(ref):.2f})"
    return f"{ref:.3g}"


def basis(mid):
    if mid in SE.ABS_BOUND:
        return BASIS.get(mid, "± 4 deg")
    return f"± {SE.TOL:.0%} of reference"


def score_all(P, names, scoring):
    """Rows (one per shown metric, every run scored, `scored` False for a metric
    outside the scoring) and per-run group scores over the scored metrics."""
    sc = SE.SCORINGS[scoring]
    rows, per = [], {n: {g: [] for g in "ABCD"} for n in names}
    for m0 in P[names[0]]["metrics"]:
        mid = m0["id"]
        if mid in sc["hidden"]:
            continue
        absolute, scored = mid in SE.ABS_BOUND, mid not in sc["unscored"]
        refs = {n: next(x for x in P[n]["metrics"] if x["id"] == mid)["ref"] for n in names}
        same = len({json.dumps(r) for r in refs.values()}) == 1
        r = dict(id=mid, group=m0["group"], label=m0["label"], scored=scored,
                 unit=SE.ABS_UNIT.get(mid, "deg") if absolute else "%", basis=basis(mid),
                 reference=ref_text(m0["ref"]) + ("" if same else " (first run; per model grid)"),
                 runs={})
        for n in names:
            m = next(x for x in P[n]["metrics"] if x["id"] == mid)
            val, err = SE.score(mid, m["model"], m["ref"])
            r["runs"][n] = dict(model=m["model"], ref=m["ref"], value=val, err=err,
                                status=SE.status(err))
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
    cfg = SE.COMPARE[set_name]
    other_set, other = cfg["against"]
    pair = {"in": other, "out": cfg["run"]}
    Q = {k: P[n] for k, n in pair.items()}
    metrics = []
    for m0 in Q["out"]["metrics"]:
        mid = m0["id"]
        e = dict(id=mid, group=m0["group"], label=m0["label"],
                 unit=SE.ABS_UNIT.get(mid, "deg") if mid in SE.ABS_BOUND else "%",
                 scored={s: mid not in d["hidden"] and mid not in d["unscored"]
                         for s, d in SE.SCORINGS.items()})
        for k in pair:
            m = next(x for x in Q[k]["metrics"] if x["id"] == mid)
            val, err = SE.score(mid, m["model"], m["ref"])
            e[k] = dict(model=m["model"], ref=m["ref"], value=val, err=err, status=SE.status(err))
        metrics.append(e)
    groups = {s: {k: score_all({pair[k]: Q[k]}, [pair[k]], s)[1][0] for k in pair}
              for s in SE.SCORINGS}
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
        mark = {s: ("✓" if e["scored"][s] else ("–" if e["id"] in SE.SCORINGS[s]["hidden"] else "◦"))
                for s in SE.SCORINGS}
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
         f"| Score | rel = model / ref − 1, err = \\|rel\\| / {SE.TOL:g}; pass err ≤ 1, near ≤ 1.5, "
         "fail above; A2, A2-WACCM, C1 latitudes, C2 and D4 final warming keep absolute bounds |",
         "| Combined | mean of the A, B and C group scores (each the mean err of its scored "
         "metrics); D reported separately |",
         "| w* estimator | jcm-strat Ψ* route (`strat_eval_common.tem_streamfunction`, "
         "`residual_w`), linear in log p to the target pressure — **interim**, AIDE-eval issue #5 |",
         "| Generated by | `scripts/02_scorecard.py`, from `output/01_fields__<run>.json` |",
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


def link(name):
    return f"[`{name}`]({name})"


def summary_md(scoring, results):
    """One page per scoring over every result set: what is scored, the group scores,
    the metrics that pass or fail in every run, the AL7h window comparison, links."""
    first = next(iter(results.values()))
    rows0 = first["rows"]
    L = [f"# Transport scorecard — {scoring} scoring", "",
         f"Summary over every result set, produced {first['produced']}. Full tables: the reports "
         "linked below. Method, contract and decisions: "
         "[docs/TRANSPORT_SCORECARD.md](../../docs/TRANSPORT_SCORECARD.md).", "",
         "## What is scored", "",
         f"- {SCORING_TEXT[scoring]}."]
    for g, lab in GROUPS:
        if g == "combined":
            continue
        ids = [r["id"] for r in rows0 if r["group"] == g and r["scored"]]
        L.append(f"- {lab}: {len(ids)} scored — {', '.join(ids)}")
    off = [r["id"] for r in rows0 if not r["scored"]]
    if off:
        L.append(f"- Listed, not scored: {', '.join(off)}")
    L += [f"- Score: rel = model / ref − 1, err = \\|rel\\| / {SE.TOL:g}; pass err ≤ 1, near ≤ 1.5, "
          "fail above. A2, A2-WACCM, C1 latitudes, C2 and D4 final warming keep absolute bounds; "
          "A4 against the midpoint 2.75 yr.",
          "- Group score = mean err of its scored metrics; combined = mean of A, B, C; D separate.",
          "- References follow each set's window (A, C, D); B is WACCM6 h6 1996–2014 for every set.",
          "- w*: jcm-strat Ψ* route, log-p interpolated — interim, AIDE-eval issue #5. "
          "Group C input from daily means — AIDE-eval issue #7."]
    for set_name, res in results.items():
        names, w = res["runs"], res["window"]
        L += ["", f"## {set_name} ({w[0]}–{w[1]})", "",
              "| | " + " | ".join(names) + " |", "|---|" + "---|" * len(names)]
        for k, lab in GROUPS:
            L.append(f"| {lab} | " + " | ".join(
                (f"**{s[k]:.2f}**" if k == "combined" else f"{s[k]:.2f}") for s in res["scores"]) + " |")
        L.append("| passed A-C | " + " | ".join(f"{s['passed']} / {s['total']}" for s in res["scores"]) + " |")
        L.append("| passed D | " + " | ".join(f"{s['D_passed']} / {s['D_total']}" for s in res["scores"]) + " |")
        scored = [r for r in res["rows"] if r["scored"]]
        L.append("")
        if len(names) == 1:
            for st, word in (("fail", "Fail"), ("near", "Near")):
                ids = [r["id"] for r in scored if r["runs"][names[0]]["status"] == st]
                L.append(f"- {word}: {', '.join(ids) if ids else 'none'}")
        else:
            for st, word in (("fail", "Fail in every run"), ("pass", "Pass in every run")):
                ids = [r["id"] for r in scored if all(c["status"] == st for c in r["runs"].values())]
                L.append(f"- {word}: {', '.join(ids) if ids else 'none'}")
        L.append(f"- Report: {link(stamped('transport_scorecard.md', res['stamp']))}; figures "
                 f"{link(stamped('transport_score_table.png', res['stamp']))}, "
                 f"{link(stamped('transport_heatmap.png', res['stamp']))}")
        cmp = res.get("comparison")
        if cmp:
            G = cmp["groups"][scoring]
            (lab_in, lab_out), wi, wo = cmp["labels"], cmp["windows"]["in"], cmp["windows"]["out"]
            L += ["", f"### {cmp['run']} against {cmp['against']['run']}: {lab_in} "
                  f"{wi[0]}–{wi[1]} vs {lab_out} {wo[0]}–{wo[1]}", "",
                  f"| | {wi[0]}–{wi[1]} | {wo[0]}–{wo[1]} |", "|---|---|---|"]
            for k, lab in GROUPS:
                f = (lambda x: f"**{x:.2f}**") if k == "combined" else (lambda x: f"{x:.2f}")
                L.append(f"| {lab} | {f(G['in'][k])} | {f(G['out'][k])} |")
            L.append(f"| passed A-C | {G['in']['passed']} / {G['in']['total']} | "
                     f"{G['out']['passed']} / {G['out']['total']} |")
            changed = [e["id"] for e in cmp["metrics"] if e["scored"][scoring]
                       and e["in"]["status"] != e["out"]["status"]]
            L += ["", f"- Scored metrics changing status: {', '.join(changed) if changed else 'none'}",
                  f"- Metric by metric: {link(stamped('transport_scorecard.md', res['stamp']))}; figure "
                  f"{link(stamped(COMPARE_FIGURE[0], res['set_stamp']))}"]
    return "\n".join(L) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1].strip())
    ap.add_argument("sets", nargs="*", help="result sets to score (default: every set)")
    ap.add_argument("--scoring", choices=list(SE.SCORINGS), help="one scoring (default: both)")
    a = ap.parse_args()
    sets = a.sets or list(SE.SETS)
    unknown = [s for s in sets if s not in SE.SETS]
    if unknown:
        raise SystemExit(f"unknown set(s) {unknown}; registered: {list(SE.SETS)}")
    scorings = [a.scoring] if a.scoring else list(SE.SCORINGS)
    produced = datetime.date.today().isoformat()
    os.makedirs(RESULTS, exist_ok=True)
    done = {s: {} for s in scorings}
    for set_name in sets:
        names = list(SE.SETS[set_name])
        window = SE.window_of(names)
        load = [*names] + ([SE.COMPARE[set_name]["against"][1]] if set_name in SE.COMPARE else [])
        P = {n: json.load(open(os.path.join(SE.OUTDIR, f"01_fields__{n}.json"))) for n in load}
        R = json.load(open(os.path.join(SE.OUTDIR,
                                        f"01_fields__references__{window[0]}-{window[1]}.json")))
        cmp = comparison(set_name, P) if set_name in SE.COMPARE else None
        for scoring in scorings:
            rows, scores = score_all(P, names, scoring)
            res = dict(name=set_name, scoring=scoring, produced=produced,
                       stamp=stamp_for(set_name, scoring, produced),
                       set_stamp=stamp_for(set_name, produced), runs=names, window=list(window),
                       tolerance=SE.TOL, a4_range=list(SE.A4_RANGE), abs_bound=SE.ABS_BOUND,
                       scoring_definition={k: list(v) for k, v in SE.SCORINGS[scoring].items()},
                       wstar_estimator="jcm-strat Psi* route, log-p interpolated; interim "
                                       "(AIDE-eval issue #5)",
                       rows=rows, scores=scores, comparison=cmp,
                       profiles={n: P[n]["profiles"] for n in names},
                       references=R)
            out_j = os.path.join(RESULTS, stamped("transport_scorecard.json", res["stamp"]))
            out_m = os.path.join(RESULTS, stamped("transport_scorecard.md", res["stamp"]))
            json.dump(res, open(out_j, "w"), indent=1)
            open(out_m, "w").write(markdown(res, P, R))
            done[scoring][set_name] = res
            print(f"{set_name}, {scoring}")
            for s in scores:
                print(f"  {s['run']:9s} A {s['A']:.2f}  B {s['B']:.2f}  C {s['C']:.2f}  combined "
                      f"{s['combined']:.2f}  passed {s['passed']}/{s['total']}  D {s['D']:.2f} "
                      f"({s['D_passed']}/{s['D_total']})")
            print(f"wrote {out_j}\nwrote {out_m}")
    if set(sets) != set(SE.SETS):
        print("summaries not written: they cover every set (run 02 without set names)")
        return
    for scoring in scorings:
        out_s = os.path.join(RESULTS, SUMMARY[scoring])
        open(out_s, "w").write(summary_md(scoring, done[scoring]))
        print(f"wrote {out_s}")


if __name__ == "__main__":
    main()
