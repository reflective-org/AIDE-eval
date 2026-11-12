"""
port_check - Check the ported diagnostics against jcm-strat's own numbers.

Evidence, not a pipeline stage: it feeds nothing. It reads jcm-strat's phase-29 and
phase-31 results, computes the diagnostics of the same runs
through 01's functions, and writes the difference tables
validation_results/port_check.md:

  1. phase 29: the eight runs, 2006-2010, group A at phase 29's 2005-2009 age years
  2. phase 31: AL7h 2006-2010 and 2011-2019, age references of the window's years,
     A2-WACCM included (success_<tag>/success_metrics.json, a2_waccm.json, and the
     group D files P31_RUNS names)

  checks/port_check.py

The runs are defined here (P29_RUNS, P31_RUNS), not read from sources/, so the check
does not change when the scored result sets do.

Two jcm-strat files are compared against:
  relative_metrics.json  the scored value (% off the reference, or the absolute
                         difference), rounded to 0.1 % / 0.01 - the brief's comparison
  success_metrics.json   the unrounded model and reference values behind it - the
                         0.5 % agreement test is applied to these

A metric agrees when both its model value and its reference value are within 0.5 %
of jcm-strat's, or within the rounding jcm-strat applied before scoring (teq_score.py
rounds jets to 0.1 m/s and 0.1 deg and the polar-cap bias to 0.01 K), whichever is
larger. B1, B2, D4-w and D7 are scored on log-p-interpolated w* and Psi*; for them
the table checks the nearest-level values (jcm-strat's method) for agreement and
lists the change the interpolation makes separately.

In part 1 group A is recomputed at jcm-strat phase 29's age reference years
(P29_AGE_YEARS, 2005-2009); the scorecard and part 2 take them over the window.
"""

import argparse
import copy
import datetime
import functools
import importlib.util
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # scripts/
sys.path.insert(0, HERE)
import strat_eval_common as SE

RESULTS = os.path.join(os.path.dirname(HERE), "validation_results")
# jcm-strat's own results and the runs they were computed from. Model files are read
# in place.
JCM = "/data/JCM_stripped"
P29_OUT = f"{JCM}/jcm-strat-phase29/docs/outputs/29_lqbo10d"
P29_RELATIVE = os.path.join(P29_OUT, "relative", "relative_metrics.json")
P29_SUCCESS = os.path.join(P29_OUT, "success_0610", "success_metrics.json")
P29_AGE_YEARS = (2005, 2009)      # phase 29's age reference years (sweep_score.CLAMS_YEARS)
P29_RUNS = dict(name="jcm-strat-phase29", reader="jcm_npz", window=[2006, 2010], notes={},
                runs={n: [f"{JCM}/jcm-strat-phase29/cache/variability/{n}_0610.npz"]
                      for n in ("P20b", "AL7", "AL7h", "AL7t", "L0", "LAL7", "QBO10d", "LQBO10d")})
# phase 31: AL7h 2006-2010 and 2011-2019, age references of the window's years, and the
# phase-30 group D of AL7h 2006-2010 that phase 31 reports.
# port source: (window, model file, success_<tag> directory, group D file)
P31_OUT = f"{JCM}/jcm-strat-phase31/docs/outputs/31_al7h_validation"
P31_RUNS = {
    "AL7h": ([2006, 2010], f"{JCM}/jcm-strat-phase29/cache/variability/AL7h_0610.npz", "success_0610",
             f"{JCM}/jcm-strat-phase30/docs/outputs/30_level_thinning/variability.json"),
    "AL7h_1119": ([2011, 2019], f"{JCM}/jcm-strat-phase31/cache/variability/AL7h_1119.npz",
                  "success_1119", os.path.join(P31_OUT, "variability.json")),
}
REL_LIMIT = 0.005
ROUNDING = {"C1": 0.05, "C2": 0.005}          # half the last digit teq_score.py keeps
GROUP_NAME = {"A": "A age of air", "B": "B residual circulation",
              "C": "C polar vortex", "D": "D variability and relations"}

# Causes of the differences that remain, established in the port (see the PR and
# docs/TRANSPORT_SCORECARD.md section 7). Printed under each group table.
NOTES = {
    "A": [],
    "B": ["Nearest-level values differ by -0.003 % (w*) and +0.008 % (mass flux) in every run "
          "and in the reference: the port uses CESM's a = 6.37122e6 m and g = 9.80616 m s-2 "
          "(strat_eval_common) where jcm-strat uses 6.371e6 m and "
          "9.80665 m s-2; w* scales with 1/a, Psi* with a/g.",
          "The scored values change because w* and Psi* are interpolated linearly in log p to "
          "the exact target pressure. Nearest levels: model 95.9 / 72.2 / 49.3 / 30.4 / 9.86 hPa, "
          "WACCM6 h6 94.9 / 66.8 / 47.1 / 32.6 / 9.55 hPa for 100 / 70 / 50 / 30 / 10 hPa. The "
          "largest change is in the WACCM6 references at 100 hPa (w* +14 %, mass flux +14 %), "
          "where w* and Psi* grow fast below the nearest level."],
    "C": ["jcm-strat teq_score.py takes calendar-month means of the 6-hourly frames by their "
          "timestamp (00Z on day 1 to 18Z on the last day). The port takes them from the daily "
          "means, and jcm-strat daily_zm.py builds each day from 06Z, 12Z, 18Z and 00Z of the next "
          "day. Every month therefore loses the 00Z frame that opens it and gains the one that "
          "opens the next month: one frame in about 120. The ERA5 side is identical. The shift "
          "table above estimates the effect from the daily series; it has the sign and size of "
          "every C2 difference. Open as AIDE-eval issue #7 (group C daily-mean input).",
          "A jet latitude is a grid latitude (argmax). Where two grid rows carry nearly the same "
          "wind, the one-frame shift can move the maximum by one row (1.87 deg). The tie table "
          "above lists the wind at both rows for every latitude that differs."],
    "D": ["D5 and D6 differ by 0.008 %: v'T' = v'theta' (p/1000)^kappa uses kappa = R/c_p = "
          "0.28573 (strat_eval_common) where jcm-strat hard-codes 0.2857. "
          "D4-w and D7 differ by 0.003 % through a (group B note).",
          "D4-w and D7 are scored on w* interpolated in log p to 70 hPa. ERA5's CDS levels include "
          "70 hPa, so its reference does not change; the model's nearest level is 72.2 hPa, so its "
          "value does (table above)."],
}


def limit(mid, value):
    return max(REL_LIMIT * abs(value), ROUNDING.get(mid[:2], 0.0))


def agrees(mid, port, jcm):
    return abs(port - jcm) <= limit(mid, jcm)


def fmt(x, nd=2):
    return "n/a" if x is None or not np.isfinite(x) else f"{x:+.{nd}f}"


# Every run the check reads, with the window it is diagnosed over
SPECS = {n: P29_RUNS for n in P29_RUNS["runs"]}
SPECS.update({n: dict(name=f"jcm-strat-phase31-{n}", reader="jcm_npz", window=w, notes={},
                      runs={n: [path]})
              for n, (w, path, _, _) in P31_RUNS.items() if n not in SPECS})


def diagnose_all():
    """01's diagnostics of every run in SPECS against the references of its window."""
    m01, refs, P = module_01(), {}, {}
    for n, spec in SPECS.items():
        w = tuple(spec["window"])
        if w not in refs:
            refs[w] = m01.references(w)
        P[n] = m01.diagnose(spec, n, refs[w])
        print(f"diagnosed {n} ({w[0]}-{w[1]})")
    return P


def load(names):
    S = json.load(open(P29_SUCCESS))
    R = {r["id"]: r for r in json.load(open(P29_RELATIVE))["rows"]}
    J = {n: {m["id"]: m for m in S[n]["metrics"] + S[n].get("D_metrics", [])} for n in names}
    return J, R


def group_table(g, names, P, J, R):
    """Scored value jcm-strat -> port, and the raw-value agreement, one row per metric
    (the metrics jcm-strat phase 29 has)."""
    ids = [m["id"] for m in P[names[0]]["metrics"] if m["group"] == g and m["id"] in R]
    L = [f"### {GROUP_NAME[g]}", "",
         "Scored value, jcm-strat `relative_metrics.json` → port (% off the reference; absolute "
         "metrics in their own unit); for B1, B2, D4-w and D7 the port value is the nearest-level "
         "one, the scored change follows in the interpolation table. ✓ = model and reference values "
         "agree with `success_metrics.json` to 0.5 % or to jcm-strat's rounding.", "",
         "| metric | " + " | ".join(names) + " |", "|---|" + "---|" * len(names)]
    n_ok = n_all = 0
    worst = {}
    for mid in ids:
        cells = []
        for n in names:
            m = next(x for x in P[n]["metrics"] if x["id"] == mid)
            pm = m.get("model_nearest", m["model"])
            pr = m.get("ref_nearest", m["ref"])
            jm, jr = J[n][mid]["model"], J[n][mid]["ref"]
            ok = agrees(mid, pm, jm) and agrees(mid, SE.ref_value(pr), SE.ref_value(jr))
            val, _ = SE.score(mid, pm, pr)
            cells.append(f"{R[mid][n]:+g} → {val:+.{1 if mid not in SE.ABS_BOUND else 2}f} "
                         f"{'✓' if ok else '✗'}")
            n_ok += ok
            n_all += 1
            d = abs(pm / jm - 1) if jm else abs(pm - jm)
            worst[mid] = max(worst.get(mid, 0.0), d)
        L.append(f"| **{mid}** {m['label']} | " + " | ".join(cells) + " |")
    L += ["", f"{n_ok} of {n_all} agree. Largest relative difference of the model value "
          "per metric: " + ", ".join(f"{k} {100 * v:.3f} %" for k, v in worst.items()) + "."]
    return L


def raw_table(g, names, P, J):
    """Unrounded model and reference values, jcm-strat against the port."""
    ids = [m["id"] for m in P[names[0]]["metrics"] if m["group"] == g and m["id"] in J[names[0]]]
    L = ["", "Unrounded values: jcm-strat `success_metrics.json` / port (nearest level where "
         "the scored value is interpolated). Reference in the last column (first run's grid).", "",
         "| metric | " + " | ".join(names) + " | reference |", "|---|" + "---|" * (len(names) + 1)]
    for mid in ids:
        cells = []
        for n in names:
            m = next(x for x in P[n]["metrics"] if x["id"] == mid)
            cells.append(f"{J[n][mid]['model']:.4g} / {m.get('model_nearest', m['model']):.4g}")
        m0 = next(x for x in P[names[0]]["metrics"] if x["id"] == mid)
        jr, pr = J[names[0]][mid]["ref"], m0.get("ref_nearest", m0["ref"])
        ref = (f"{jr} / {pr}" if isinstance(jr, list) else f"{jr:.4g} / {pr:.4g}")
        L.append(f"| {mid} | " + " | ".join(cells) + f" | {ref} |")
    return L


def interpolation_table(g, names, P):
    """What log-p interpolation changes: reference and model, per metric."""
    ids = [m["id"] for m in P[names[0]]["metrics"] if m["group"] == g and "model_nearest" in m]
    if not ids:
        return []
    L = ["", "Change made by log-p interpolation to the exact target pressure (the scored "
         "values) against the nearest level (jcm-strat). Reference: nearest → interpolated "
         "(change). Runs: change of the model value, and of the scored value in % points.", "",
         "| metric | reference | " + " | ".join(names) + " |", "|---|---|" + "---|" * len(names)]
    for mid in ids:
        m0 = next(x for x in P[names[0]]["metrics"] if x["id"] == mid)
        rn, ri = m0["ref_nearest"], m0["ref"]
        cells = []
        for n in names:
            m = next(x for x in P[n]["metrics"] if x["id"] == mid)
            sn, _ = SE.score(mid, m["model_nearest"], m["ref_nearest"])
            si, _ = SE.score(mid, m["model"], m["ref"])
            cells.append(f"{100 * (m['model'] / m['model_nearest'] - 1):+.1f} % ({si - sn:+.1f} pp)")
        L.append(f"| {mid} | {rn:.4g} → {ri:.4g} ({100 * (ri / rn - 1):+.1f} %) | "
                 + " | ".join(cells) + " |")
    return L


def module(name):
    spec = importlib.util.spec_from_file_location(
        "m" + name[:2], os.path.join(HERE, name))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


@functools.lru_cache(maxsize=None)
def module_01():
    return module("01_fields.py")


def phase29_age(names, P):
    """Replace each run's group A by A against the 2005-2009 age references."""
    m01 = module_01()
    R = m01.references_A(SE.years_of(P29_AGE_YEARS))
    for n in names:
        F = SE.load_run(SPECS[n], n)
        rows, _ = m01.group_A(F, R)
        P[n]["metrics"] = rows + [m for m in P[n]["metrics"] if m["group"] != "A"]
    return R


def phase31_metrics(run):
    """jcm-strat phase 31's values for one port source: A-C from success_metrics.json
    (A2 against shifted CLaMS), A2-WACCM from a2_waccm.json, D from the group D file."""
    _, _, tag, dfile = P31_RUNS[run]
    S = json.load(open(os.path.join(P31_OUT, tag, "success_metrics.json")))["AL7h"]
    a2 = json.load(open(os.path.join(P31_OUT, tag, "a2_waccm.json")))
    D = json.load(open(dfile))["runs"]["AL7h"]["metrics"]
    J = {m["id"]: m for m in S["metrics"] + D}
    J["A2-WACCM"] = dict(model=a2["rmse_waccm"], ref=0.0, years=a2["years"])
    return J, [os.path.join(P31_OUT, tag, f) for f in ("success_metrics.json", "a2_waccm.json")] + [dfile]


def rel(port, jcm):
    """Port - jcm-strat: relative in % (absolute where jcm-strat's value is 0)."""
    if jcm == 0:
        return f"{port - jcm:+.2g}"
    return f"{100 * (port / jcm - 1):+.3f} %"


def phase31_table(run, P, J):
    """Every metric of one source against jcm-strat phase 31: model and reference value,
    nearest level for B, D4-w and D7."""
    w = P[run]["window"]
    L = ["", f"### AL7h {w[0]}-{w[1]} (port source `{run}`)", "",
         "| metric | model jcm-strat | model port | Δ | reference jcm-strat | reference port | Δ | agree |",
         "|---|---|---|---|---|---|---|---|"]
    count = {g: [0, 0] for g in "ABCD"}
    for m in P[run]["metrics"]:
        mid = m["id"]
        pm, pr = m.get("model_nearest", m["model"]), SE.ref_value(m.get("ref_nearest", m["ref"]))
        jm, jr = J[mid]["model"], SE.ref_value(J[mid]["ref"])
        ok = agrees(mid, pm, jm) and agrees(mid, pr, jr)
        count[m["group"]][0] += ok
        count[m["group"]][1] += 1
        L.append(f"| **{mid}** | {jm:.5g} | {pm:.5g} | {rel(pm, jm)} | {jr:.5g} | {pr:.5g} | "
                 f"{rel(pr, jr)} | {'✓' if ok else '✗'} |")
    L += ["", "Agreeing: " + ", ".join(f"{g} {a} / {b}" for g, (a, b) in count.items()) + "."]
    return L, count


def phase31(P):
    """Part 2: AL7h against jcm-strat phase 31, both windows."""
    L = ["", "## Part 2 - jcm-strat phase 31: AL7h 2006-2010 and 2011-2019", "",
         "Age references (A1-A4, A2-WACCM) of the window's own years, as phase 31 "
         "(`SCORE_REF_YEARS`) and the committed scorecard; ERA5 and CLaMS of the same days for C "
         "and D; WACCM6 h6 1996-2014 for B. Phase 31 scored with absolute bounds, so only the "
         "values are compared, with the criterion of part 1. Δ = port / jcm-strat − 1 (absolute "
         "where jcm-strat's value is 0).", ""]
    for run in P31_RUNS:
        J, files = phase31_metrics(run)
        L += ["", f"jcm-strat files: " + ", ".join(f"`{f}`" for f in files)]
        tab, _ = phase31_table(run, P, J)
        L += tab
        L += c_frame_shift([run], P, {run: J})
    L += ["", "Causes of the remaining differences:", "",
          "- A, B, D: as in part 1 (CESM a, g and κ; ≤ 0.01 %). A2-WACCM reproduces "
          "`a2_waccm.json` to 1e-5 %.",
          "- C: the one-frame-per-month shift of part 1 (AIDE-eval issue #7); the shift table "
          "estimates each C2 difference from the daily series. AL7h 2006-2010 is the same run "
          "and the same model values as in part 1; phase 31 changed only its age references.",
          "- D, 2006-2010: phase 31 reports AL7h's group D from phase 30 "
          "(`30_level_thinning/variability.json`), whose AL7h entry is identical to phase 29's."]
    return L


def c_frame_shift(names, P, J):
    """Group C: the expected port - jcm-strat difference of the winter polar-cap T from the
    one-frame-per-month shift, and the wind at both rows of a jet latitude that differs.

    The month's mean changes by (x[00Z, day 1 of next month] - x[00Z, day 1]) / n_frames.
    The 00Z frames are not in the daily files; the mean of the two daily values either
    side of each stands in for it (the first and last day of the record at the edges)."""
    m01 = module_01()
    L = ["", "Expected effect of the one-frame-per-month shift on C2 against the observed "
         "difference (port - jcm-strat; jcm-strat rounds C2 to 0.01 K, so the observed value "
         "carries +-0.005 K) [K]:", "",
         "| run | C2-N observed | C2-N expected | C2-S observed | C2-S expected |",
         "|---|---|---|---|---|"]
    ties = []
    for n in names:
        F = SE.load_run(SPECS[n], n)
        p, lat, days = F["plev"], F["lat"], F["days"]
        k50 = SE.nearest_level(p, 50.0)
        ym = days.astype("datetime64[M]")
        keys = np.unique(ym)
        mon = np.array([int(str(k)[5:7]) for k in keys])
        cells = []
        for h, (months, sign) in m01.WINTER.items():
            cap = sign * lat >= 60
            w = np.cos(np.deg2rad(lat[cap]))
            tc = (F["T"][:, k50][:, cap] * w).sum(1, dtype="f8") / w.sum()
            sh = []
            for k in keys[np.isin(mon, months)]:
                i = np.where(ym == k)[0]
                a = tc[i[0]] if i[0] == 0 else 0.5 * (tc[i[0] - 1] + tc[i[0]])
                b = tc[i[-1]] if i[-1] + 1 >= days.size else 0.5 * (tc[i[-1]] + tc[i[-1] + 1])
                sh.append((b - a) / (4 * i.size))
            obs = next(x for x in P[n]["metrics"] if x["id"] == f"C2-{h}")["model"] - J[n][f"C2-{h}"]["model"]
            cells += [f"{obs:+.3f}", f"{np.mean(sh):+.3f}"]
        L.append(f"| {n} | " + " | ".join(cells) + " |")
        keys_m, U = m01.monthly(F["Uzm"], days)
        sel = np.array([k[1] for k in keys_m])
        for mid in [x["id"] for x in P[n]["metrics"] if x["id"].startswith("C1") and x["id"].endswith("-lat")]:
            pm = next(x for x in P[n]["metrics"] if x["id"] == mid)["model"]
            jm = J[n][mid]["model"]
            if agrees(mid, pm, jm):
                continue
            lev, h = mid[6:8], mid[9]
            months, _ = m01.WINTER[h]
            um = U[np.isin(sel, months)].mean(0)[SE.nearest_level(p, float(lev))]
            jj = int(np.argmin(np.abs(lat - jm)))
            jp = int(np.argmin(np.abs(lat - pm)))
            ties.append(f"| {n} | {mid} | {lat[jj]:.2f}: {um[jj]:.3f} | {lat[jp]:.2f}: {um[jp]:.3f} | "
                        f"{um[jp] - um[jj]:+.4f} |")
    if ties:
        L += ["", "Jet latitudes that differ: the port's winter-mean wind at jcm-strat's row and at "
              "its own row [m/s]:", "",
              "| run | metric | jcm-strat row | port row | difference |", "|---|---|---|---|---|"] + ties
    return L


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1].strip())
    ap.parse_args()
    names = list(P29_RUNS["runs"])
    window = P29_RUNS["window"]
    P_all = diagnose_all()
    J, R = load(names)
    P = {n: copy.deepcopy(P_all[n]) for n in names}
    phase29_age(names, P)                         # group A against 2005-2009
    groups = [g for g in "ABCD" if any(m["group"] == g for m in P[names[0]]["metrics"])]
    L = ["# Port check: transport diagnostics against jcm-strat phases 29 and 31", "",
         f"Generated by `scripts/checks/port_check.py` on {datetime.date.today()}, from the "
         "diagnostics of `01_fields.py` computed in the check. Do not edit by hand.", "",
         f"- jcm-strat: `{P29_RELATIVE}` and `{P29_SUCCESS}`",
         f"- Runs: {', '.join(names)}; window {window[0]}-{window[1]}",
         f"- Group A against CLaMS and WACCM6 age references of {P29_AGE_YEARS[0]}-"
         f"{P29_AGE_YEARS[1]}, jcm-strat phase 29's years (recomputed here; the scorecard "
         f"uses the window's own years)",
         f"- Groups ported: {', '.join(groups)}",
         "- Agreement: |port - jcm-strat| <= max(0.5 % of the jcm-strat value, the rounding "
         "jcm-strat applied before scoring), for both the model and the reference value", "",
         "## Part 1 - jcm-strat phase 29: eight runs, 2006-2010", ""]
    for g in groups:
        L += group_table(g, names, P, J, R)
        L += interpolation_table(g, names, P)
        L += raw_table(g, names, P, J)
        if g == "C":
            L += c_frame_shift(names, P, J)
        if NOTES[g]:
            L += ["", "Causes of the remaining differences:", ""] + [f"- {x}" for x in NOTES[g]]
        L.append("")
    L += phase31({n: P_all[n] for n in P31_RUNS})
    os.makedirs(RESULTS, exist_ok=True)
    out = os.path.join(RESULTS, "port_check.md")
    open(out, "w").write("\n".join(L) + "\n")
    print("\n".join(L))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
