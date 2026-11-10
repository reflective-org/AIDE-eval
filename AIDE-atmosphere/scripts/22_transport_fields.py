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


# --------------------------------------------------------------------- group D
# Variability and two relations against ERA5 (CDS daily zonal means) and CLaMS for
# the same days (jcm-strat variability_score.py, quantities and errors). Every
# anomaly is taken against the dataset's own smoothed day-of-year climatology.
# Levels are the nearest in log p, as jcm-strat's Z.ilev. The monthly tropical w*
# of D4-w and D7 goes through the group-B estimator, interpolated to 70 hPa
# (nearest-level values kept for the port check).
HEMI = {"N": dict(lat=60.0, band=(55, 65), win=(12, 1, 2), rel=(1, 2, 3), hf=(45, 75),
                  cap=(60, 90), sign=1.0),
        "S": dict(lat=-60.0, band=(-65, -55), win=(6, 7, 8), rel=(7, 8, 9, 10), hf=(-75, -45),
                  cap=(-90, -60), sign=-1.0)}


class Daily:
    """Daily zonal means (day, lev, lat) of one source, jcm-strat variability_score.Z."""

    def __init__(self, days, p, lat, fields):
        self.days, self.p, self.lat, self.f = days, p, lat, fields

    @classmethod
    def from_contract(cls, F):
        return cls(F["days"], F["plev"], F["lat"],
                   dict(u=F["Uzm"], v=F["Vzm"], T=F["T"], th=F["THzm"], vth=F["VTHzm"],
                        aoa150=F["AOA150"]))

    @classmethod
    def era5(cls, years):
        zs = [dict(np.load(T.REFS["era5_daily"].format(year=y))) for y in years]
        cat = lambda k: np.concatenate([z[k] for z in zs])
        return cls(cat("days").astype("datetime64[D]"), zs[0]["p_hpa"], zs[0]["lat"],
                   {k: cat(k) for k in ("u", "v", "T", "th", "vth")})

    def ilev(self, p):
        return int(np.argmin(np.abs(np.log(self.p / p))))

    def band(self, var, p, lo, hi):
        x = self.f[var][:, self.ilev(p)]
        m = (self.lat >= lo) & (self.lat <= hi)
        w = np.cos(np.deg2rad(self.lat[m]))
        return (x[:, m] * w).sum(1) / w.sum()

    def at_lat(self, var, p, lat):
        x = self.f[var][:, self.ilev(p)]
        return np.array([np.interp(lat, self.lat, r) for r in x])

    def vT(self, p, lo, hi):
        """v'T' from v'theta' on the pressure level."""
        return self.band("vth", p, lo, hi) / (C.P_REF_HPA / self.p[self.ilev(p)]) ** C.KAPPA

    def months(self):
        return (self.days.astype("datetime64[Y]").astype(int) + 1970,
                self.days.astype("datetime64[M]").astype(int) % 12 + 1)


def doy(days):
    """Day of year 0..364, 29 February folded onto 28 February."""
    d = (days - days.astype("datetime64[Y]")).astype(int)
    leap = (days.astype("datetime64[Y]").astype(int) + 1970) % 4 == 0
    return np.where(leap & (d >= 59), d - 1, d).clip(0, 364)


def anomaly(x, days, smooth=31):
    """x minus its own day-of-year climatology, smoothed by a 31-day running mean."""
    k = doy(days)
    clim = np.array([np.nanmean(x[k == i]) if np.any(k == i) else np.nan for i in range(365)])
    clim = np.where(np.isnan(clim), np.nanmean(clim), clim)
    pad = np.concatenate([clim[-smooth:], clim, clim[:smooth]])
    clim = np.convolve(pad, np.ones(smooth) / smooth, "same")[smooth:-smooth]
    return x - clim[k]


def monthly_tem_wstar(z, plev=70.0):
    """Monthly tropical w* (15S-15N) from monthly means of the daily [v], [v'theta'],
    [theta]: months, interpolated to plev, at the nearest level."""
    ym = z.days.astype("datetime64[M]")
    t, wi, wn = [], [], []
    for m in np.unique(ym):
        k = ym == m
        if k.sum() < 25:
            continue
        mean = lambda v: z.f[v][k].mean(0, dtype="f8")
        p, _, psi = T.tem_streamfunction(mean("v"), mean("vth"), mean("th"), z.p * 100.0, z.lat)
        w = T.residual_w(p, psi, z.lat)
        t.append(m)
        wi.append(T.tropical_wstar(p / 100.0, w, z.lat, plev))
        wn.append(T.tropical_wstar(p / 100.0, w, z.lat, plev, nearest=True))
    return np.array(t), np.array(wi), np.array(wn)


def deseason_detrend(t_month, x):
    m = t_month.astype(int) % 12
    clim = np.array([x[m == i].mean() for i in range(12)])
    a = x - clim[m]
    tt = np.arange(a.size)
    return a - np.polyval(np.polyfit(tt, a, 1), tt)


def clams_monthly_age(years, plev=55.0):
    """CLaMS monthly AGE [yr], 10S-10N, linear in log p to plev."""
    t, a = [], []
    for y in years:
        with T.clams_member(y) as d:
            A = np.ma.filled(np.ma.asarray(d.variables["AGE"][:], dtype="f8"), np.nan)
            press = np.asarray(d.variables["press"][:], dtype="f8")
            lat = np.asarray(d.variables["lat"][:], dtype="f8")
        A = np.where(A > -1e20, A, np.nan)
        lp = np.log(press)
        o = np.argsort(lp)
        m = np.abs(lat) <= 10
        w = np.cos(np.deg2rad(lat[m]))
        for mi in range(12):
            r = np.array([np.interp(np.log(plev), lp[o], A[mi][o, j]) for j in range(lat.size)])
            t.append(np.datetime64(f"{y}-{mi + 1:02d}"))
            a.append(np.nansum(r[m] * w) / np.nansum(w * np.isfinite(r[m])))
    return np.array(t), np.array(a)


def ssw_dates(u60n, days, min_sep_days=20):
    """Major NH SSW central dates (jcm-strat strat_compare.ssw_dates): the first
    easterly day Nov-Mar, >= 20 days after the previous event, and not a final
    warming (10 consecutive westerly days before 30 April)."""
    t, u = np.asarray(days).astype("datetime64[D]"), np.asarray(u60n)
    events, last = [], None
    for i in range(1, len(t)):
        month = int(str(t[i])[5:7])
        if u[i] < 0 <= u[i - 1] and month in (11, 12, 1, 2, 3):
            if last is not None and (t[i] - last).astype(int) < min_sep_days:
                continue
            yr = int(str(t[i])[:4]) + (1 if month >= 11 else 0)
            sel = (t > t[i]) & (t <= np.datetime64(f"{yr}-04-30"))
            run_start, ok = None, False
            for v, tt in zip(u[sel], t[sel]):
                if v > 0:
                    run_start = tt if run_start is None else run_start
                    if (tt - run_start).astype(int) >= 10:
                        ok = True
                        break
                else:
                    run_start = None
            if ok:
                events.append(t[i])
                last = t[i]
    return events


def quantities(z, age=None):
    """Every D quantity of one source. age = (months, tropical 55 hPa age) for a
    reference, None to take the source's own AOA150."""
    q = {}
    y, mo = z.months()
    for h, H in HEMI.items():
        u = z.at_lat("u", 10.0, H["lat"])
        # D1: winter means of u 55-65 deg, complete winters (DJF labelled by its January)
        wy = np.where(mo == 12, y + 1, y) if h == "N" else y
        win = np.isin(mo, H["win"])
        ub = z.band("u", 10.0, *H["band"])
        vals = [ub[win & (wy == yy)].mean() for yy in np.unique(wy[win])
                if (win & (wy == yy)).sum() >= 85]
        q[f"D1-{h}"] = float(np.std(vals, ddof=1)) if len(vals) > 1 else np.nan
        q[f"winters-{h}"] = len(vals)
        # D2: daily spread of u at 60 deg within winter
        q[f"D2-{h}"] = float(np.std(anomaly(u, z.days)[win]))
        # D5: daily spread of the 45-75 deg heat flux at 100 hPa, poleward positive
        hfa = anomaly(H["sign"] * z.vT(100.0, *H["hf"]), z.days)
        q[f"D5-{h}"] = float(np.std(hfa[win]))
        # D6: 50 hPa polar-cap T anomaly on the trailing 40-day mean heat-flux anomaly
        Ta = anomaly(z.band("T", 50.0, *H["cap"]), z.days)
        c = np.concatenate([[0.0], np.cumsum(hfa)])
        x = np.full(hfa.size, np.nan)
        x[39:] = (c[40:] - c[:-40]) / 40
        k = np.isin(mo, H["rel"]) & np.isfinite(x)
        q[f"D6-{h}"] = float(np.polyfit(x[k], Ta[k], 1)[0])
        q[f"D6r-{h}"] = float(np.corrcoef(x[k], Ta[k])[0, 1])
        if h == "N":                                                  # D3
            q["SSW-N"] = [str(e) for e in ssw_dates(u, z.days)]
        # D4: final warming - first easterly crossing of the 5-day-smoothed u after
        # which it stays below 5 m/s (Black & McDaniel), Mar-Jul N / Sep-Jan S
        us = np.full(u.size, np.nan)
        us[2:-2] = np.convolve(u, np.ones(5) / 5, "valid")
        fw = []
        for yy in np.unique(y):
            start = np.datetime64(f"{yy}-03-01") if h == "N" else np.datetime64(f"{yy}-09-01")
            end = np.datetime64(f"{yy}-08-01") if h == "N" else np.datetime64(f"{yy + 1}-02-01")
            k = (z.days >= start) & (z.days < end)
            if k.sum() < 120:
                continue
            s, dd = us[k], z.days[k]
            ok = np.isfinite(s)
            s, dd = s[ok], dd[ok]
            cross = [i for i in np.where((s[:-1] >= 0) & (s[1:] < 0))[0] if np.all(s[i + 1:] < 5.0)]
            if cross:
                fw.append((int(yy), int((dd[cross[0] + 1] - start).astype(int))))
        q[f"fw-{h}"] = fw
    # D4-w: annual-cycle amplitude of tropical w* at 70 hPa; D7: age vs w*
    tm, w_i, w_n = monthly_tem_wstar(z, 70.0)
    if age is None:
        ym = z.days.astype("datetime64[M]")
        a55 = z.band("aoa150", 55.0, -10, 10)
        at, av = tm, np.array([a55[ym == m].mean() for m in tm])
    else:
        at, av = age
        keep = np.isin(at, tm)
        at, av = at[keep], av[keep]
    aa = deseason_detrend(at, av)
    for tag, w in (("", w_i), ("_nearest", w_n)):
        mi = tm.astype(int) % 12
        clim = np.array([w[mi == i].mean() for i in range(12)])
        q[f"D4w{tag}"] = float(clim.max() - clim.min())
        q[f"w70_clim{tag}"] = clim.tolist()
        wa = deseason_detrend(tm, w)
        w3 = np.convolve(wa, np.ones(3) / 3, "full")[:wa.size]
        w3[:2] = np.nan
        ww = np.array([w3[np.where(tm == m)[0][0]] for m in at])
        k = np.isfinite(ww)
        q[f"D7{tag}"] = float(np.polyfit(ww[k], aa[k], 1)[0])
        q[f"D7r{tag}"] = float(np.corrcoef(ww[k], aa[k])[0, 1])
    return q


def references_D():
    years = list(range(T.WINDOW[0], T.WINDOW[1] + 1))
    q = quantities(Daily.era5(years), age=clams_monthly_age(years))
    return dict(years=years, q=q)


def group_D(F, R):
    m, r = quantities(Daily.from_contract(F)), R["q"]
    rows = []
    for h in "NS":
        rows.append(row("D", f"D1-{h}", f"year-to-year std of winter u 55-65{h} 10 hPa [m/s]",
                        m[f"D1-{h}"], r[f"D1-{h}"]))
    for h in "NS":
        rows.append(row("D", f"D2-{h}", f"daily std of winter u 60{h} 10 hPa [m/s]",
                        m[f"D2-{h}"], r[f"D2-{h}"]))
    rows.append(row("D", "D3", "NH major SSWs in the window [number]",
                    len(m["SSW-N"]), len(r["SSW-N"]), dates=m["SSW-N"], ref_dates=r["SSW-N"]))
    for h in "NS":
        dm, dr = dict(map(tuple, m[f"fw-{h}"])), dict(map(tuple, r[f"fw-{h}"]))
        yy = sorted(set(dm) & set(dr))                                # common springs
        a = float(np.mean([dm[y] for y in yy])) if yy else np.nan
        b = float(np.mean([dr[y] for y in yy])) if yy else np.nan
        rows.append(row("D", f"D4-fw{h}", f"final warming {h}, days after "
                        f"{'1 Mar' if h == 'N' else '1 Sep'}", a, b, springs=yy))
    rows.append(row("D", "D4-w", "annual-cycle amplitude of tropical w* 70 hPa [mm/s]",
                    m["D4w"], r["D4w"], model_nearest=m["D4w_nearest"], ref_nearest=r["D4w_nearest"]))
    for h in "NS":
        rows.append(row("D", f"D5-{h}", f"daily std of 100 hPa heat flux {h} [K m/s]",
                        m[f"D5-{h}"], r[f"D5-{h}"]))
    for h in "NS":
        rows.append(row("D", f"D6-{h}", f"heat flux -> polar-cap T slope {h} [K per K m/s]",
                        m[f"D6-{h}"], r[f"D6-{h}"], r_model=m[f"D6r-{h}"], r_ref=r[f"D6r-{h}"]))
    rows.append(row("D", "D7", "age vs w* slope [yr per mm/s]", m["D7"], r["D7"],
                    model_nearest=m["D7_nearest"], ref_nearest=r["D7_nearest"],
                    r_model=m["D7r"], r_ref=r["D7r"]))
    return rows, dict(w70_clim=m["w70_clim"])


# --------------------------------------------------------------------- driver
def references():
    return dict(A=references_A(), B=references_B(), C={}, D=references_D())


def public(d):
    """Drop the in-memory fields (keys starting with _) before writing JSON."""
    return {k: (public(v) if isinstance(v, dict) else v)
            for k, v in d.items() if not k.startswith("_")}


def diagnose(name, R):
    spec = T.SOURCES[name]
    F = T.load_source(name, YEARS)
    rows, prof = [], {}
    for fn, key in ((group_A, "A"), (group_B, "B"), (group_C, "C"), (group_D, "D")):
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
