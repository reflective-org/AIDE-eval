"""
02 - Score the diagnostics and write the scorecard reports.

The second stage of the scorecard (docs/EVALUATION_PROTOCOL.md). Reads
the per-source output of 01 and scores every metric as jcm-strat
relative_scores.py does:

  rel = model / ref - 1,  err = |rel| / 0.25        pass err <= 1, near <= 1.5, fail above
  A2, A2-WACCM, C1 latitudes, C2 and D4 final-warming dates keep their absolute bound
  A4 is scored against the midpoint 2.75 yr of its 2.3-3.2 yr range

Group score = mean err of the group's scored metrics; combined = mean of the A, B
and C group scores; D is reported on its own.

Every result set (sources/<set>.json) is scored twice (strat_eval_common.SCORINGS):
  standard   the headline score and ranking: jcm-strat phase 31's scored set (the
             former "reduced" scoring), A2-WACCM for A2, 15 metrics shown but not scored
  extensive  the full listing: the 34 metrics of the phase-29 scorecard
A set whose window is not the standard 2006-2019 is flagged "non-standard window
(N of 14 standard years)" in its reports and in both summaries.

  02_scorecard.py [SET ...] [--scoring standard|extensive]
                                         default: every set, both scorings

After scoring, the two summaries are rewritten from the scorecard JSON of every set
registered in sources/ (so scoring one set updates them).

Output: validation_results/<set>/scorecard_<scoring>.json, scorecard_<scoring>.md
        validation_results/SCORING_STANDARD.md, SCORING_EXTENSIVE.md
            one page per scoring comparing every set
"""

import argparse
import datetime
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import strat_eval_common as SE

RESULTS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "validation_results")
SYMBOL = {"pass": "✅", "near": "⚠️", "fail": "❌"}
BASIS = {"A2": "< 0.4 yr", "A2-WACCM": "< 0.4 yr", "C2-N": "± 2 K", "C2-S": "± 2 K",
         "D4-fwN": "± 10 days", "D4-fwS": "± 10 days"}
GROUPS = (("A", "A age of air"), ("B", "B residual circulation"), ("C", "C polar vortex"),
          ("combined", "**combined A-C**"), ("D", "D variability"))
SCORING_TEXT = {
    "standard": "**standard** (the headline score and ranking): jcm-strat phase 31's scored set, "
                "the former \"reduced\" scoring — A2-WACCM replaces A2; A2 and "
                f"{', '.join(SE.PHASE31_UNSCORED)} are listed but not scored. Phase 31 scored with "
                "absolute bounds; here the relative scoring applies to that set",
    "extensive": "**extensive** (the full listing): the 34 metrics of the phase-29 scorecard, "
                 "every one scored",
}
# Figures of a set (03_figures.py); {scoring} filled per report
FIGURES = [
    ("score_table_{scoring}.png", "group scores and passes per run"),
    ("heatmap_{scoring}.png", "every metric, % off its reference (absolute metrics as a "
     "fraction of their bound), ±25 % marked; unscored metrics faded"),
    ("wstar_profile.png", "tropical w* profiles against WACCM6"),
    ("age_profile.png", "tropical entry-age profiles against CLaMS and WACCM6"),
]
SUMMARY = {"standard": "SCORING_STANDARD.md", "extensive": "SCORING_EXTENSIVE.md"}


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


def markdown(res, P, R):
    names = res["runs"]
    first = P[names[0]]
    A, B, D = R["A"], R["B"], R["D"]
    flag = f" — **{res['window_flag']}**" if res["window_flag"] else ""
    L = [f"# Transport scorecard: {res['name']}, {res['scoring']} scoring", ""]
    if res["window_flag"]:
        L += [f"> **{res['window_flag'][0].upper() + res['window_flag'][1:]}**: scored on "
              f"{res['window'][0]}–{res['window'][1]}, not on the standard window "
              f"{SE.STANDARD_WINDOW[0]}–{SE.STANDARD_WINDOW[1]}.", ""]
    L += ["| | |", "|---|---|",
         f"| Scored | {', '.join(names)} |",
         f"| Read by | reader `{'`, `'.join(sorted({P[n]['reader'] for n in names}))}`; "
         f"{'; '.join(sorted({os.path.dirname(f) for n in names for f in P[n]['paths']}))} |",
         f"| Window | {first['days'][0]} – {first['days'][1]} ({first['n_days']} days){flag} |",
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
         "| Generated by | `scripts/02_scorecard.py`, from `output/<set>/<run>.json` |",
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
    L += ["", "## Figures", ""]
    L += [f"- {link(f.format(scoring=res['scoring']))} — {d}" for f, d in FIGURES]
    L += ["", "Method, contract and decisions: "
          "[docs/EVALUATION_PROTOCOL.md](../../docs/EVALUATION_PROTOCOL.md). Agreement with "
          "jcm-strat's own numbers: [port_check.md](../port_check.md)."]
    return "\n".join(L) + "\n"


def link(name):
    return f"[`{name}`]({name})"


def summary_md(scoring, results, missing):
    """One page per scoring over every result set: what is scored, the group scores,
    the metrics that pass or fail in every run, links. `missing`: registered sets with
    no scorecard yet."""
    first = next(iter(results.values()))
    rows0 = first["rows"]
    L = [f"# Transport scorecard — {scoring} scoring", "",
         "Summary over every result set registered in `sources/`. Full tables: the reports "
         "linked below. Method, contract and decisions: "
         "[docs/EVALUATION_PROTOCOL.md](../docs/EVALUATION_PROTOCOL.md).", "",
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
    if scoring == "standard":
        ranked = sorted(((s["combined"], set_name, s) for set_name, res in results.items()
                         for s in res["scores"]), key=lambda x: (x[0], x[1], x[2]["run"]))
        L += ["", "## Ranking", "",
              "Every run of every set by its combined A–C score (lower is better). Runs on "
              "different windows are compared with references of their own years.", "",
              "| rank | run | set | window | combined A–C | A | B | C | D |",
              "|---|---|---|---|---|---|---|---|---|"]
        for i, (_, set_name, s) in enumerate(ranked, 1):
            res = results[set_name]
            w = f"{res['window'][0]}–{res['window'][1]}"
            w += f" — {res['window_flag']}" if res["window_flag"] else ""
            L.append(f"| {i} | {s['run']} | {set_name} | {w} | **{s['combined']:.2f}** | "
                     f"{s['A']:.2f} | {s['B']:.2f} | {s['C']:.2f} | {s['D']:.2f} |")
    else:
        L += ["", f"Ranking: the standard scoring ([{SUMMARY['standard']}]({SUMMARY['standard']}))."]
    for set_name, res in results.items():
        names, w = res["runs"], res["window"]
        flag = f" — **{res['window_flag']}**" if res["window_flag"] else ""
        L += ["", f"## {set_name} ({w[0]}–{w[1]}){flag}", "",
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
        d = f"{set_name}/"
        L.append(f"- Produced {res['produced']}. Report: {link(d + f'scorecard_{scoring}.md')}; figures "
                 f"{link(d + f'score_table_{scoring}.png')}, {link(d + f'heatmap_{scoring}.png')}")
    if missing:
        L += ["", f"Registered in `sources/`, not yet scored: {', '.join(missing)}."]
    return "\n".join(L) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1].strip())
    ap.add_argument("sets", nargs="*", help="result sets to score (default: every set)")
    ap.add_argument("--scoring", choices=list(SE.SCORINGS), help="one scoring (default: both)")
    a = ap.parse_args()
    sets = a.sets or SE.set_names()
    scorings = [a.scoring] if a.scoring else list(SE.SCORINGS)
    produced = datetime.date.today().isoformat()
    for set_name in sets:
        spec = SE.load_set(set_name)
        names, window = list(spec["runs"]), spec["window"]
        P = {n: json.load(open(os.path.join(SE.OUTDIR, set_name, f"{n}.json"))) for n in names}
        R = json.load(open(os.path.join(SE.OUTDIR, set_name, "references.json")))
        out_dir = os.path.join(RESULTS, set_name)
        os.makedirs(out_dir, exist_ok=True)
        for scoring in scorings:
            rows, scores = score_all(P, names, scoring)
            res = dict(name=set_name, description=spec["description"], scoring=scoring,
                       produced=produced, runs=names, window=list(window),
                       window_flag=SE.window_flag(window),
                       tolerance=SE.TOL, a4_range=list(SE.A4_RANGE), abs_bound=SE.ABS_BOUND,
                       scoring_definition={k: list(v) for k, v in SE.SCORINGS[scoring].items()},
                       wstar_estimator="jcm-strat Psi* route, log-p interpolated; interim "
                                       "(AIDE-eval issue #5)",
                       rows=rows, scores=scores,
                       profiles={n: P[n]["profiles"] for n in names},
                       references=R)
            out_j = os.path.join(out_dir, f"scorecard_{scoring}.json")
            out_m = os.path.join(out_dir, f"scorecard_{scoring}.md")
            json.dump(res, open(out_j, "w"), indent=1)
            open(out_m, "w").write(markdown(res, P, R))
            print(f"{set_name}, {scoring}")
            for s in scores:
                print(f"  {s['run']:9s} A {s['A']:.2f}  B {s['B']:.2f}  C {s['C']:.2f}  combined "
                      f"{s['combined']:.2f}  passed {s['passed']}/{s['total']}  D {s['D']:.2f} "
                      f"({s['D_passed']}/{s['D_total']})")
            print(f"wrote {out_j}\nwrote {out_m}")
    write_summaries(scorings)


def write_summaries(scorings):
    """Both summaries over every registered set that has a scorecard on disk."""
    for scoring in scorings:
        results, missing = {}, []
        for set_name in SE.set_names():
            f = os.path.join(RESULTS, set_name, f"scorecard_{scoring}.json")
            if os.path.exists(f):
                results[set_name] = json.load(open(f))
            else:
                missing.append(set_name)
        if not results:
            continue
        out_s = os.path.join(RESULTS, SUMMARY[scoring])
        open(out_s, "w").write(summary_md(scoring, results, missing))
        print(f"wrote {out_s}")


if __name__ == "__main__":
    main()
