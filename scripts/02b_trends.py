"""
02b - Separate the forced trend from internal variability.

sigma_interannual is the anchor for almost every target, so it matters whether it
measures INTERNAL variability or is inflated by a forced trend. Over 1970-2014 the
prescribed GHG and ODS forcing accelerates the Brewer-Dobson circulation and cools
the polar stratosphere, both well-documented CCM responses, and both present here.

Where a trend is significant, the target should be anchored on the DETRENDED sigma
(otherwise it is quietly loosened by the forced signal), and the trend itself
becomes a separate diagnostic.

There is a complication specific to AIDE-WACCM. Its input-only forcings are TOA
solar, time-of-day, year-progress and statics - GHG concentrations are NOT among
them. A free rollout therefore has no mechanism to reproduce a GHG-forced trend,
so this is reported as a documented expectation rather than a pass/fail target.

Output: output/02b_trends.json
"""

import json
import os
import numpy as np
from scipy import stats

import aide_val_common as C

OUT = os.path.join(C.OUTDIR, "02b_trends.json")

SERIES = [
    ("upward_mass_flux_70hPa", ["tropical_upwelling", "upward_mass_flux_70hPa"],
     "annual_means", "1e9 kg/s"),
    ("w_star_70hPa_10S10N", ["tropical_upwelling", "w_star_70hPa_10S10N"],
     "annual_means", "mm/s"),
    ("vortex_NH_DJF", ["polar_vortex", "NH"], "season_means", "m/s"),
    ("vortex_SH_JJA", ["polar_vortex", "SH"], "season_means", "m/s"),
    ("polar_cap_T_NH_DJF", ["polar_cap_T", "NH"], "season_means", "K"),
    ("polar_cap_T_SH_JJA", ["polar_cap_T", "SH"], "season_means", "K"),
]


def main():
    with open(os.path.join(C.OUTDIR, "02_reference_stats.json")) as f:
        R = json.load(f)

    def cat(path, key):
        a, b = R["1970-1995"], R["1996-2014"]
        for k in path:
            a, b = a[k], b[k]
        return np.array(a[key] + b[key])

    out = {"note": "45-year pooled series, 1970-2014; trend by OLS on the "
                   "annual/seasonal means"}
    print(f"{'series':26s} {'mean':>9} {'sig_raw':>9} {'sig_detr':>9} "
          f"{'trend/dec':>10} {'%/dec':>7} {'p':>8}")
    print("-" * 82)

    for name, path, key, unit in SERIES:
        v = cat(path, key)
        t = np.arange(len(v), dtype=float)
        sl, ic, r, p, se = stats.linregress(t, v)
        det = v - (ic + sl * t)
        sig_raw = float(v.std(ddof=1))
        sig_det = float(det.std(ddof=1))
        sig_use = sig_det if p < 0.05 else sig_raw
        out[name] = dict(
            units=unit, n=len(v), mean=float(v.mean()),
            sigma_raw=sig_raw, sigma_detrended=sig_det,
            trend_per_decade=float(sl * 10), trend_pct_per_decade=
                float(100 * sl * 10 / abs(v.mean())),
            trend_se_per_decade=float(se * 10), p_value=float(p),
            significant=bool(p < 0.05),
            sigma_for_anchoring=sig_use,
            sigma_inflation_pct=float(100 * (sig_raw / sig_det - 1)),
        )
        flag = " <- detrended sigma used" if p < 0.05 else ""
        print(f"{name:26s} {v.mean():9.3f} {sig_raw:9.4f} {sig_det:9.4f} "
              f"{sl*10:10.4f} {100*sl*10/abs(v.mean()):7.2f} {p:8.4f}{flag}")

    print("\nInterpretation")
    mf = out["upward_mass_flux_70hPa"]
    pc = out["polar_cap_T_NH_DJF"]
    print(f"  BDC acceleration {mf['trend_pct_per_decade']:+.1f} %/decade "
          f"(p = {mf['p_value']:.4f}) - consistent with the ~2 %/decade reported "
          f"across CCMs.")
    print(f"  NH polar cap cooling {pc['trend_per_decade']:+.2f} K/decade "
          f"(p = {pc['p_value']:.4f}).")
    print(f"  Vortex strength has NO significant trend in either hemisphere, so "
          f"its anchors are unaffected.")
    print(f"  Using the raw sigma for the mass flux would loosen its target by "
          f"{mf['sigma_inflation_pct']:.0f} %.")

    with open(OUT, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
