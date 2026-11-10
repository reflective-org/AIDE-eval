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
import importlib.util
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
    "B": ["Nearest-level values differ by -0.003 % (w*) and +0.008 % (mass flux) in every run "
          "and in the reference: the port uses CESM's a = 6.37122e6 m and g = 9.80616 m s-2 "
          "(aide_val_common, protocol section 4.6) where jcm-strat uses 6.371e6 m and "
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
          "every C2 difference.",
          "A jet latitude is a grid latitude (argmax). Where two grid rows carry nearly the same "
          "wind, the one-frame shift can move the maximum by one row (1.87 deg). The tie table "
          "above lists the wind at both rows for every latitude that differs."],
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


def module_22():
    spec = importlib.util.spec_from_file_location(
        "t22", os.path.join(os.path.dirname(os.path.abspath(__file__)), "22_transport_fields.py"))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def c_frame_shift(names, P, J):
    """Group C: the expected port - jcm-strat difference of the winter polar-cap T from the
    one-frame-per-month shift, and the wind at both rows of a jet latitude that differs.

    The month's mean changes by (x[00Z, day 1 of next month] - x[00Z, day 1]) / n_frames.
    The 00Z frames are not in the daily files; the mean of the two daily values either
    side of each stands in for it (the first and last day of the record at the edges)."""
    m22 = module_22()
    L = ["", "Expected effect of the one-frame-per-month shift on C2 against the observed "
         "difference (port - jcm-strat; jcm-strat rounds C2 to 0.01 K, so the observed value "
         "carries +-0.005 K) [K]:", "",
         "| run | C2-N observed | C2-N expected | C2-S observed | C2-S expected |",
         "|---|---|---|---|---|"]
    ties = []
    for n in names:
        F = T.load_source(n, range(T.WINDOW[0], T.WINDOW[1] + 1))
        p, lat, days = F["plev"], F["lat"], F["days"]
        k50 = T.nearest_level(p, 50.0)
        ym = days.astype("datetime64[M]")
        keys = np.unique(ym)
        mon = np.array([int(str(k)[5:7]) for k in keys])
        cells = []
        for h, (months, sign) in m22.WINTER.items():
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
        keys_m, U = m22.monthly(F["Uzm"], days)
        sel = np.array([k[1] for k in keys_m])
        for mid in [x["id"] for x in P[n]["metrics"] if x["id"].startswith("C1") and x["id"].endswith("-lat")]:
            pm = next(x for x in P[n]["metrics"] if x["id"] == mid)["model"]
            jm = J[n][mid]["model"]
            if agrees(mid, pm, jm):
                continue
            lev, h = mid[6:8], mid[9]
            months, _ = m22.WINTER[h]
            um = U[np.isin(sel, months)].mean(0)[T.nearest_level(p, float(lev))]
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
        if g == "C":
            L += c_frame_shift(names, P, J)
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
