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
import xarray as xr

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


# --------------------------------------------------------------------- group B
# Residual circulation (jcm-strat sweep_score.wstar_fields / waccm_wstar and
# success_metrics B1, B2): Psi* and w* from the time means of [v], [v'theta'] and
# [theta] over the window, model and WACCM6 through the same functions. Scored at
# the exact target pressure, linear in log p; the nearest-level values jcm-strat
# scores are kept beside them for the port check. INTERIM w* estimator, AIDE-eval
# issue #5 (docs/TRANSPORT_SCORECARD.md).
W_LEVELS = (100.0, 70.0, 50.0, 30.0, 10.0)
FLUX_LEVELS = (100.0, 70.0)


def circulation(V, VTH, TH, p_hpa, lat):
    """Time-mean fields (time, p, lat) -> p [hPa, ascending], Psi* [kg/s], w* [mm/s]."""
    mean = lambda x: np.nanmean(x, axis=0, dtype="f8")
    p, _, psi = T.tem_streamfunction(mean(V), mean(VTH), mean(TH), p_hpa * 100.0, lat)
    return p / 100.0, psi, T.residual_w(p, psi, lat)


def circulation_metrics(p, psi, w, lat):
    out = {}
    for lv in W_LEVELS:
        out[f"B1-{lv:g}"] = (T.tropical_wstar(p, w, lat, lv),
                             T.tropical_wstar(p, w, lat, lv, nearest=True))
    for lv in FLUX_LEVELS:
        out[f"B2-{lv:g}"] = (T.upward_mass_flux(p, psi, lat, lv) / 1e9,
                             T.upward_mass_flux(p, psi, lat, lv, nearest=True) / 1e9)
    return out


def references_B():
    seg = T.REFS["waccm_tem_segment"]
    V, VTH, TH = (C.load_h6(seg, v, T.WACCM_ILEV) for v in ("Vzm", "VTHzm", "THzm"))
    yr = np.array([t.year for t in V.time.values])
    keep = (yr >= T.WACCM_TEM_YEARS[0]) & (yr <= T.WACCM_TEM_YEARS[1])
    lat = V.lat.values.astype("f8")
    p, psi, w = circulation(V.values[keep], VTH.values[keep], TH.values[keep],
                            V.ilev.values.astype("f8"), lat)
    return dict(years=list(T.WACCM_TEM_YEARS), segment=seg, n_days=int(keep.sum()),
                n_levels=int(p.size), top_hPa=float(p[0]),
                metrics={k: dict(interp=v[0], nearest=v[1])
                         for k, v in circulation_metrics(p, psi, w, lat).items()},
                nearest_hPa={f"{lv:g}": float(p[T.nearest_level(p, lv)]) for lv in W_LEVELS},
                profile_waccm_p=p.tolist(),
                profile_waccm_wstar_tropics=T.tropical_profile(w, lat).tolist())


def group_B(F, R):
    p, lat = F["plev"], F["lat"]
    p, psi, w = circulation(F["Vzm"], F["VTHzm"], F["THzm"], p, lat)
    M = circulation_metrics(p, psi, w, lat)
    rows = []
    for mid, (vi, vn) in M.items():
        lv = mid.split("-")[1]
        label = (f"w* {lv} hPa [mm/s]" if mid.startswith("B1")
                 else f"mass flux {lv} hPa [1e9 kg/s]")
        r = R["metrics"][mid]
        rows.append(row("B", mid, label, vi, r["interp"], model_nearest=vn, ref_nearest=r["nearest"],
                        nearest_hPa=float(p[T.nearest_level(p, float(lv))])))
    return rows, dict(wstar_p=p.tolist(), wstar_tropics=T.tropical_profile(w, lat).tolist())


# --------------------------------------------------------------------- group C
# Polar vortex against ERA5 for the same months (jcm-strat teq_score.py): monthly
# means of the model's daily u and T, ERA5 monthly zonal means interpolated onto
# the model's latitudes (linear) and levels (linear in log p, clamped to ERA5's
# range), then DJF (north) and JJA (south) means over the window's months. At the
# model level nearest 50 and 10 hPa, as jcm-strat: ERA5 sits on the same level,
# so the comparison is like for like. The reference depends on the model grid and
# is computed per source.
WINTER = {"N": ((12, 1, 2), 1.0), "S": ((6, 7, 8), -1.0)}


def monthly(x, days):
    """Calendar-month means of a daily (day, p, lat) field -> keys, (month, p, lat)."""
    ym = days.astype("datetime64[M]")
    keys = np.unique(ym)
    return ([(int(str(k)[:4]), int(str(k)[5:7])) for k in keys],
            np.stack([x[ym == k].mean(axis=0, dtype="f8") for k in keys]))


def era5_monthly(var, keys, p, lat):
    """ERA5 monthly zonal mean of `var` on (p, lat) for each (year, month) in keys."""
    out = {}
    for y in sorted({k[0] for k in keys}):
        with xr.open_dataset(T.REFS["era5_zm_monthly"].format(year=y)) as ds:
            tname = [c for c in ds.coords if c.startswith("time") or c == "valid_time"][0]
            x = ds[var].rename({tname: "time"}).interp(lat=lat)
            lp = np.log(ds.level.values.astype(float))
            o = np.argsort(lp)
            lpo = np.log(np.clip(p, np.exp(lp[o][0]), np.exp(lp[o][-1])))
            for m in range(1, 13):
                if (y, m) in keys:
                    prof = x.sel(time=x.time.dt.month == m).mean("time").values[o]
                    out[(y, m)] = np.stack([np.interp(lpo, lp[o], prof[:, j])
                                            for j in range(lat.size)], axis=1)
    return np.stack([out[k] for k in keys])


def jet(u, lat, sign):
    """Winter jet maximum poleward of 30 deg and its (grid) latitude."""
    m = sign * lat >= 30
    i = int(np.argmax(u[m]))
    return float(u[m][i]), float(lat[m][i])


def group_C(F, R):
    p, lat = F["plev"], F["lat"]
    keys, U = monthly(F["Uzm"], F["days"])
    _, TT = monthly(F["T"], F["days"])
    EU, ET = era5_monthly("uzm", keys, p, lat), era5_monthly("tzm", keys, p, lat)
    mon = np.array([k[1] for k in keys])
    k50, k10 = T.nearest_level(p, 50.0), T.nearest_level(p, 10.0)
    jets, caps = [], []
    for h, (months, sign) in WINTER.items():
        k = np.isin(mon, months)
        um, eu, tm, et = U[k].mean(0), EU[k].mean(0), TT[k].mean(0), ET[k].mean(0)
        name = {"N": "NH", "S": "SH"}[h]
        for lev, kk in (("50", k50), ("10", k10)):
            (mj, ml), (ej, el) = jet(um[kk], lat, sign), jet(eu[kk], lat, sign)
            jets.append((f"jet{lev}_{h}", f"{name} {lev} hPa", mj, ml, ej, el))
        cap = sign * lat >= 60
        w = np.cos(np.deg2rad(lat[cap]))
        caps.append(row("C", f"C2-{h}", f"polar-cap T bias 50 hPa {h} [K]",
                        np.sum((tm[k50] - et[k50])[cap] * w) / np.sum(w), 0.0,
                        model_T=float(np.sum(tm[k50][cap] * w) / np.sum(w)),
                        era5_T=float(np.sum(et[k50][cap] * w) / np.sum(w))))
    rows = []
    for key in ("jet50_N", "jet10_N", "jet50_S", "jet10_S"):
        _, lab, mj, ml, ej, el = next(j for j in jets if j[0] == key)
        rows.append(row("C", f"C1-{key}", f"jet {lab} [m/s]", mj, ej))
        rows.append(row("C", f"C1-{key}-lat", f"jet {lab} latitude [deg]", ml, el))
    return rows + caps, dict(level_50=float(p[k50]), level_10=float(p[k10]), n_months=len(keys))


# --------------------------------------------------------------------- driver
def references():
    return dict(A=references_A(), B=references_B(), C={})


def public(d):
    """Drop the in-memory fields (keys starting with _) before writing JSON."""
    return {k: (public(v) if isinstance(v, dict) else v)
            for k, v in d.items() if not k.startswith("_")}


def diagnose(name, R):
    spec = T.SOURCES[name]
    F = T.load_source(name, YEARS)
    rows, prof = [], {}
    for fn, key in ((group_A, "A"), (group_B, "B"), (group_C, "C")):
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
