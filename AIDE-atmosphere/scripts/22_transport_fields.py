"""
22 - Read a source's transport fields and compute the transport diagnostics.

The first stage of the transport scorecard (docs/TRANSPORT_SCORECARD.md). A
source is read through transport_common.READERS into the field contract of
section 1 of that document, then reduced to the metrics of groups A-D, each
beside its reference. Scoring is 23's; this script only measures.

  22_transport_fields.py [NAME ...]          diagnose the named sources (default:
                                             every registered source), one at a time
  22_transport_fields.py --inspect NAME      shapes, units and ranges of the
                                             contract fields, nothing written

The model-independent references (CLaMS, WACCM6) are computed once per call and
written beside the per-source results.

Output: output/22_transport__references.json
        output/22_transport__<NAME>.json      one per source
"""

import argparse
import json
import os
import sys
import warnings

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import aide_val_common as C
import transport_common as T

YEARS = range(T.WINDOW[0], T.WINDOW[1] + 1)
AGE_YEARS = list(range(T.AGE_REF_YEARS[0], T.AGE_REF_YEARS[1] + 1))


def row(group, mid, label, model, ref, **extra):
    return dict(group=group, id=mid, label=label, model=float(model),
                ref=(list(ref) if isinstance(ref, tuple) else float(ref)), **extra)


# --------------------------------------------------------------------- group A
# Age of air from the entry-age clock AOA150 (jcm-strat sweep_score.score, "age"
# block, and success_metrics.metrics A1-A4). The model clock is the mean over
# every day of the window. Levels are the nearest model level to 55 hPa, as in
# jcm-strat; A4 interpolates the tropical profile in log p.
def references_A():
    pc, latc, ac = T.clams_age(AGE_YEARS)
    pw, latw, aw = T.waccm_age(AGE_YEARS)
    # CLaMS counts from the surface; its own tropical age at 150 hPa is the offset
    # that turns it into an entry age comparable with AOA150
    offset = T.band_abs(latc, ac[T.nearest_level(pc, 150.0)], 0, 10)
    kc, kw = T.nearest_level(pc, 55.0), T.nearest_level(pw, 55.0)
    tc, tw = np.abs(latc) <= 10, np.abs(latw) <= 10
    with warnings.catch_warnings():              # all-fill columns below the surface
        warnings.simplefilter("ignore", RuntimeWarning)
        prof_c = np.nanmean(ac[:, tc], axis=1) - offset
        prof_w = np.nanmean(aw[:, tw], axis=1)
    return dict(
        years=AGE_YEARS, clams_offset_150=offset,
        waccm_55_tropics=T.band_abs(latw, aw[kw], 0, 10),
        clams_55_tropics=T.band_abs(latc, ac[kc], 0, 10),
        clams_55_5070=T.band_abs(latc, ac[kc], 50, 70),
        profile_clams_p=pc.tolist(), profile_clams_entry_tropics=prof_c.tolist(),
        profile_waccm_p=pw.tolist(), profile_waccm_entry_tropics=prof_w.tolist(),
        _clams=(pc, latc, ac))


def group_A(F, R):
    p, lat = F["plev"], F["lat"]
    age = F["AOA150"].mean(axis=0, dtype="f8")                     # (p, lat), yr
    pc, latc, ac = R["_clams"]
    off = R["clams_offset_150"]

    k55 = T.nearest_level(p, 55.0)
    trop, ext = T.band_abs(lat, age[k55], 0, 10), T.band_abs(lat, age[k55], 50, 70)

    # A2: model on the CLaMS grid against CLaMS minus its 150 hPa tropical age,
    # 100-5 hPa, |lat| <= 80
    am = T.on_grid(age, p, lat, pc, latc)
    dom = ((pc >= 5.0) & (pc <= 100.0))[:, None] & (np.abs(latc) <= 80.0)[None, :]
    rmse = T.weighted_rmse(am - (ac - off), latc, dom)

    prof = np.nanmean(age[:, np.abs(lat) <= 10], axis=1)
    at = lambda lv: float(np.interp(np.log(lv), np.log(p), prof))
    rows = [
        row("A", "A1", "entry age 55 hPa tropics [yr]", trop, R["waccm_55_tropics"]),
        row("A", "A2", "2-D age RMSE vs CLaMS [yr]", rmse, 0.0),
        row("A", "A3", "age step tropics -> 50-70 deg, 55 hPa [yr]", ext - trop,
            R["clams_55_5070"] - R["clams_55_tropics"]),
        row("A", "A4", "transit 70 -> 10 hPa [yr]", at(10.0) - at(70.0), T.A4_RANGE),
    ]
    return rows, dict(age_p=p.tolist(), age_tropics=prof.tolist(), age_level_55=float(p[k55]))


# --------------------------------------------------------------------- driver
def references():
    return dict(A=references_A())


def public(d):
    """Drop the in-memory fields (keys starting with _) before writing JSON."""
    return {k: (public(v) if isinstance(v, dict) else v)
            for k, v in d.items() if not k.startswith("_")}


def diagnose(name, R):
    spec = T.SOURCES[name]
    F = T.load_source(name, YEARS)
    rows, prof = [], {}
    for fn, key in ((group_A, "A"),):
        r, p = fn(F, R[key])
        rows += r
        prof.update(p)
    return dict(name=name, long_name=spec["long_name"], note=spec.get("note", ""),
                reader=spec["reader"], path=spec.get("path", ""),
                independent=bool(spec["independent"]),
                days=[str(F["days"][0]), str(F["days"][-1])], n_days=int(F["days"].size),
                n_levels=int(F["plev"].size), n_lats=int(F["lat"].size),
                metrics=rows, profiles=prof)


def inspect(name):
    """Print what the reader hands the diagnostics, so a new source can be checked
    against the contract before anything is scored."""
    spec = T.SOURCES[name]
    F = T.load_source(name, YEARS)
    p, lat, days = F["plev"], F["lat"], F["days"]
    print(f"{name}: {spec['long_name']}  reader {spec['reader']}  {spec['path']}")
    print(f"  days  {days.size}  {days[0]} .. {days[-1]}")
    print(f"  plev  {p.size} levels  {p.min():.4g} .. {p.max():.4g} hPa, ascending "
          f"{bool(np.all(np.diff(p) > 0))}; inside 5-150 hPa: {int(((p >= 5) & (p <= 150)).sum())}")
    print(f"  lat   {lat.size}  {lat.min():.3f} .. {lat.max():.3f}, ascending "
          f"{bool(np.all(np.diff(lat) > 0))}")
    for k in T.CONTRACT + T.OPTIONAL:
        if k not in F:
            print(f"  {k:7s} absent")
            continue
        x = F[k]
        print(f"  {k:7s} {str(x.shape):18s} {str(x.dtype):8s} [{T.UNITS[k]:5s}] "
              f"min {np.nanmin(x):11.4g}  mean {np.nanmean(x, dtype='f8'):11.4g}  "
              f"max {np.nanmax(x):11.4g}  NaN {int(np.isnan(x).sum())}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1].strip())
    ap.add_argument("names", nargs="*", help="sources to diagnose (default: all)")
    ap.add_argument("--inspect", metavar="NAME", help="print the contract fields of a source")
    a = ap.parse_args()
    if a.inspect:
        inspect(a.inspect)
        return
    names = a.names or list(T.SOURCES)
    unknown = [n for n in names if n not in T.SOURCES]
    if unknown:
        raise SystemExit(f"unknown source(s) {unknown}; registered: {list(T.SOURCES)}")
    os.makedirs(C.OUTDIR, exist_ok=True)
    R = references()
    out = os.path.join(C.OUTDIR, "22_transport__references.json")
    json.dump(dict(window=list(T.WINDOW), refs={k: T.REFS[k] for k in T.REFS},
                   **public(R)), open(out, "w"), indent=1)
    print(f"wrote {out}")
    for n in names:
        res = diagnose(n, R)
        out = os.path.join(C.OUTDIR, f"22_transport__{n}.json")
        json.dump(res, open(out, "w"), indent=1)
        print(f"{n}: " + "  ".join(f"{m['id']} {m['model']:.3g}" for m in res["metrics"]))
        print(f"wrote {out}")


if __name__ == "__main__":
    main()
