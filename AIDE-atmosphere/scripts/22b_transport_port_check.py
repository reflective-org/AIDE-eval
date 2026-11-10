"""
22b - Check the ported transport diagnostics against jcm-strat's own numbers.

Evidence, not a pipeline stage (like 01b-01d): it feeds nothing. It reads the
per-run output of 22 and jcm-strat's phase-29 result, and writes the difference
table validation_results/transport/port_check.md.

  22b_transport_port_check.py [NAME ...]     default: every run in transport_common.SOURCES

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
"""

import argparse
import datetime
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import aide_val_common as C
import transport_common as T

RESULTS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "validation_results", "transport")
REL_LIMIT = 0.005
ROUNDING = {"C1": 0.05, "C2": 0.005}          # half the last digit teq_score.py keeps
GROUP_NAME = {"A": "A age of air", "B": "B residual circulation",
              "C": "C polar vortex", "D": "D variability and relations"}

# Causes of the differences that remain, established in the port (see the PR and
# docs/TRANSPORT_SCORECARD.md section 6). Printed under each group table.
NOTES = {
    "A": [],
    "B": [],
    "C": [],
    "D": [],
}


def limit(mid, value):
    return max(REL_LIMIT * abs(value), ROUNDING.get(mid[:2], 0.0))


def agrees(mid, port, jcm):
    return abs(port - jcm) <= limit(mid, jcm)


def fmt(x, nd=2):
    return "n/a" if x is None or not np.isfinite(x) else f"{x:+.{nd}f}"


def load(names):
    P = {n: json.load(open(os.path.join(C.OUTDIR, f"22_transport__{n}.json"))) for n in names}
    S = json.load(open(T.JCM_STRAT_SUCCESS))
    R = {r["id"]: r for r in json.load(open(T.JCM_STRAT_RELATIVE))["rows"]}
    J = {n: {m["id"]: m for m in S[n]["metrics"] + S[n].get("D_metrics", [])} for n in names}
    return P, J, R


def group_table(g, names, P, J, R):
    """Scored value jcm-strat -> port, and the raw-value agreement, one row per metric."""
    ids = [m["id"] for m in P[names[0]]["metrics"] if m["group"] == g]
    L = [f"## {GROUP_NAME[g]}", "",
         "Scored value, jcm-strat `relative_metrics.json` → port (% off the reference; absolute "
         "metrics in their own unit). ✓ = model and reference values agree with "
         "`success_metrics.json` to 0.5 % or to jcm-strat's rounding.", "",
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
            ok = agrees(mid, pm, jm) and agrees(mid, T.ref_value(pr), T.ref_value(jr))
            val, _ = T.score(mid, pm, pr)
            cells.append(f"{R[mid][n]:+g} → {val:+.{1 if mid not in T.ABS_BOUND else 2}f} "
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
    ids = [m["id"] for m in P[names[0]]["metrics"] if m["group"] == g]
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
            sn, _ = T.score(mid, m["model_nearest"], m["ref_nearest"])
            si, _ = T.score(mid, m["model"], m["ref"])
            cells.append(f"{100 * (m['model'] / m['model_nearest'] - 1):+.1f} % ({si - sn:+.1f} pp)")
        L.append(f"| {mid} | {rn:.4g} → {ri:.4g} ({100 * (ri / rn - 1):+.1f} %) | "
                 + " | ".join(cells) + " |")
    return L


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1].strip())
    ap.add_argument("names", nargs="*")
    a = ap.parse_args()
    names = a.names or list(T.SOURCES)
    P, J, R = load(names)
    groups = [g for g in "ABCD" if any(m["group"] == g for m in P[names[0]]["metrics"])]
    L = ["# Port check: transport diagnostics against jcm-strat phase 29", "",
         f"Generated by `scripts/22b_transport_port_check.py` on {datetime.date.today()} from "
         "`output/22_transport__<run>.json`. Do not edit by hand.", "",
         f"- jcm-strat: `{T.JCM_STRAT_RELATIVE}` and `{T.JCM_STRAT_SUCCESS}`",
         f"- Runs: {', '.join(names)}; window {T.WINDOW[0]}-{T.WINDOW[1]}",
         f"- Groups ported so far: {', '.join(groups)}",
         "- Agreement: |port - jcm-strat| <= max(0.5 % of the jcm-strat value, the rounding "
         "jcm-strat applied before scoring), for both the model and the reference value", ""]
    for g in groups:
        L += group_table(g, names, P, J, R)
        L += interpolation_table(g, names, P)
        L += raw_table(g, names, P, J)
        if NOTES[g]:
            L += ["", "Causes of the remaining differences:", ""] + [f"- {x}" for x in NOTES[g]]
        L.append("")
    os.makedirs(RESULTS, exist_ok=True)
    out = os.path.join(RESULTS, "port_check.md")
    open(out, "w").write("\n".join(L) + "\n")
    print("\n".join(L))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
