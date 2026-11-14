# Evaluation protocol — stratospheric dynamics

The standard evaluation of stratospheric dynamics and transport in this repository: one
scorecard of age of air (A), residual circulation (B), polar vortex (C) and variability (D),
each metric scored as its relative deviation from an observation-based or model reference.

| | |
|---|---|
| Standard window | 2006–2019, daily (§2) |
| Headline score | **standard** scoring: jcm-strat phase 31's scored set (the former "reduced" scoring), 20 metrics (§6) |
| Full listing | **extensive** scoring: the 34 metrics of the jcm-strat phase-29 scorecard (§6) |
| Tolerance | ±25 % of the reference; absolute bounds for A2, A2-WACCM, C1 latitudes, C2, D4 dates |
| Evaluate a model | one `sources/<set>.json`, then `scripts/evaluate.py <set>` (§9) |
| Results | [validation_results/](../validation_results/): [SCORING_STANDARD.md](../validation_results/SCORING_STANDARD.md) (ranking), [SCORING_EXTENSIVE.md](../validation_results/SCORING_EXTENSIVE.md), one folder per set (§10) |
| Open issues | w* estimator interim, AIDE-eval #5 (§7); group C daily-mean input, AIDE-eval #7 (§8) |
| Code, data flow | [ARCHITECTURE.md](ARCHITECTURE.md); decisions [DECISIONS.md](DECISIONS.md); deferred [DEFERRED.md](DEFERRED.md); status [PROGRESS.md](PROGRESS.md) |

Tables and bullets state what was measured. Inferences are in labelled `Interpretation`
blocks.

---

## 1. Purpose and scope

- **Purpose.** Score a model's stratospheric circulation and tracer transport against
  references of the same years, with one fixed metric set, tolerance and code path, so that
  models and model versions are comparable.
- **Origin.** The jcm-strat relative scorecard (phase 29) and phase 31's scored metric set,
  ported into this repository. The diagnostics are recomputed here from daily zonal means;
  the port reproduces jcm-strat's own values (§10.3).
- **Scope.** Zonal-mean diagnostics, 150–5 hPa, both hemispheres. No chemistry, no
  tropospheric metric, no longitudinal structure.
- **Not in scope.** The two-tier CESM-anchored protocol that this document replaced on
  2026-11-12 (§12).

---

## 2. Standard window

- **2006-01-01 to 2019-12-31**, 14 years, daily.
- References A, C and D are taken over the scored window's own years (§5); B is WACCM6 h6
  1996–2014 for every window.
- A result set on any other window is scored with references of its own years and flagged
  **"non-standard window (N of 14 standard years)"** — N = years inside 2006–2019 — in its
  reports, both summaries, the ranking and its figures.
- Reasons for 2006–2019: the ERA5 copies held here cover 2006–2019; CLaMS v3.1 covers
  2004–2023 and the WACCM6 entry age 1970–2019; jcm-strat phase 31 extended AL7h to 2019.

| Set | Window | Flag |
|---|---|---|
| `jcm-strat-AL7h-2006-2019` | 2006–2019 | standard |
| `jcm-strat-phase29-2006-2010` | 2006–2010 | non-standard window (5 of 14 standard years) — no later years exist for these runs |

---

## 3. Input field contract

A model supplies daily zonal means on pressure levels for **every day of 2006–2019**.

| Name | Units | Definition | Metrics that read it |
|---|---|---|---|
| `Uzm` | m s⁻¹ | zonal-mean zonal wind | C1, D1–D4 |
| `Vzm` | m s⁻¹ | zonal-mean meridional wind | B1, B2, D4-w, D7 |
| `THzm` | K | zonal-mean potential temperature, reference 1000 hPa | B1, B2, D4-w, D7; C2 and D6 through T |
| `VTHzm` | K m s⁻¹ | zonal-mean **eddy** heat flux v'θ' | B1, B2, D4-w, D5, D6, D7 |
| `AOA150` | yr | entry age: an age clock reset to zero wherever p > 150 hPa | A1–A4, A2-WACCM, D7 |
| `T` | K | zonal-mean temperature — optional; if absent `T = THzm (p/1000)^κ` | C2, D6 |
| `Wzm` | m s⁻¹ | log-pressure w̄ = −Hω/p, H = 7 km — optional, read by no metric | — |

- **Shape** `(time, plev, lat)`, every field on one common set of levels and latitudes.
- **`plev`** in hPa, any order (sorted ascending on reading). Pressure per level must not
  depend on surface pressure.
- **`lat`** in degrees, global, any order.
- **`time`** daily, each value labelled by the calendar day of its own mean; any CF
  calendar (`noleap` read through `cftime`). The record must cover 1 January of the first
  year to 31 December of the last; the reader refuses anything shorter.
- **`VTHzm` is the eddy flux**, formed from instantaneous 3-D fields before any time
  averaging, not the product of the zonal means.
- **`AOA150` in years** (365.25 days).
- **Levels** must reach well above 10 hPa: Ψ* is integrated from the top supplied level
  (jcm-strat 0.01 hPa, CESM h6 4.5 × 10⁻⁶ hPa). The metrics read 5–150 hPa.
- **Fill values** above 1e20 in absolute value are masked on reading.

### Supplying netCDF

- One or more files (globs allowed) with variables named as above, dims `time`, `plev`,
  `lat`, CF time units. Read by the `netcdf` reader.
- Register the files in one JSON (§9); check them with
  `01_fields.py --inspect <set>` (shapes, units, ranges, levels in 5–150 hPa, NaN counts).
- The jcm-strat runs are read by the `jcm_npz` reader from jcm-strat's `daily_zm.py` files
  (`u`→`Uzm`, `v`→`Vzm`, `th`→`THzm`, `vth`→`VTHzm`, `aoa150`/365.25→`AOA150`, `T` kept).
  A list of files is one record: identical grid, consecutive days, no gap or overlap.

---

## 4. Metrics

Window means unless stated. "Nearest level" is the model level nearest the target pressure,
as jcm-strat takes it; w* and Ψ* are interpolated in log p (§7). Reference values for the
two scored windows, from `output/<set>/references.json` and the per-run JSON (C references
are on the jcm-strat grid).

| ID | Quantity | Reference | 2006–2019 | 2006–2010 | Scored |
|---|---|---|---|---|---|
| **A1** | entry age 10°S–10°N (cos-weighted), nearest level to 55 hPa | WACCM6 REF-D1 entry age | 1.103 yr | 1.099 yr | relative |
| **A2** | cos-weighted RMSE of AOA150 against CLaMS minus CLaMS' 10°S–10°N age at 150 hPa, model on the CLaMS grid, 100–5 hPa, \|lat\| ≤ 80° | 0 (offset 0.090 / 0.091 yr) | — | — | absolute < 0.4 yr; extensive only |
| **A2-WACCM** | cos-weighted RMSE of AOA150 against the WACCM6 REF-D1 2-D entry age, model on the WACCM6 grid, 100–5 hPa, \|lat\| ≤ 80°, no offset | 0 | — | — | absolute < 0.4 yr; standard only |
| **A3** | age step 10°S–10°N → 50–70° (both hemispheres), nearest level to 55 hPa | CLaMS | 2.83 yr | 2.84 yr | relative |
| **A4** | transit 70 → 10 hPa of the 10°S–10°N age profile, linear in log p | range 2.3–3.2 yr (CLaMS 3.0, WACCM6 2.5) | 2.75 | 2.75 | relative to the midpoint |
| **B1-p** | w* 15°S–15°N (cos-weighted) at p = 100, 70, 50, 30, 10 hPa | WACCM6 h6 1996–2014, same estimator | 0.456, 0.228, 0.198, 0.271, 0.463 mm s⁻¹ | same | relative |
| **B2-p** | tropical upward mass flux: max − min of Ψ* over \|lat\| ≤ 60° at p = 100, 70 hPa | WACCM6, same | 12.39, 6.53 × 10⁹ kg s⁻¹ | same | relative |
| **C1-jetP_H** | winter (DJF N, JJA S) maximum of monthly-mean u poleward of 30°, nearest level to P = 50, 10 hPa | ERA5 monthly zonal means, same months, model grid | 50 hPa N 17.5, S 47.7; 10 hPa N 27.6, S 73.9 m s⁻¹ | 16.4, 46.7; 25.2, 71.9 | relative |
| **C1-…-lat** | grid latitude of that maximum | ERA5, same | — | — | absolute ± 4° on \|lat\| |
| **C2-H** | winter polar-cap (60–90°) T at the nearest level to 50 hPa, minus ERA5 | 0 | — | — | absolute ± 2 K |
| **D1-H** | std over winters of winter-mean u, 10 hPa, 55–65° | ERA5 CDS daily, same days | N 9.16, S 2.88 m s⁻¹ | 6.51, 2.16 | relative |
| **D2-H** | std of daily u anomalies at 60°, 10 hPa, within winter | ERA5 | N 16.67, S 5.83 m s⁻¹ | 14.99, 5.25 | relative |
| **D3** | NH major SSWs (u 60°N 10 hPa reversal Nov–Mar, ≥ 20 days apart, not a final warming) | ERA5 | 12 | 8 | relative |
| **D4-fwH** | mean final-warming date: first easterly crossing of the 5-day-smoothed u at 60°, 10 hPa after which u stays below 5 m s⁻¹ (Mar–Jul N, Sep–Jan S), springs common to model and ERA5 | ERA5 | N 40.7, S 80.4 days after 1 Mar / 1 Sep | 44.4, 89.8 | absolute ± 10 days |
| **D4-w** | annual-cycle amplitude of monthly w* 15°S–15°N at 70 hPa | ERA5 | 0.173 mm s⁻¹ | 0.185 | relative |
| **D5-H** | std of daily v'T' anomalies, 100 hPa, 45–75°, within winter | ERA5 | N 11.88, S 7.42 K m s⁻¹ | 11.47, 7.90 | relative |
| **D6-H** | slope of 50 hPa polar-cap T anomaly on the trailing 40-day 100 hPa heat-flux anomaly (Jan–Mar N, Jul–Oct S) | ERA5 | N 1.026, S 1.135 K per K m s⁻¹ | 0.972, 0.766 | relative |
| **D7** | slope of deseasonalised, detrended monthly 10°S–10°N age at 55 hPa on the trailing 3-month w* anomaly at 70 hPa | CLaMS age with ERA5 w* | −4.41 yr per mm s⁻¹ | −4.52 | relative |

- H = N or S. D anomalies are taken against each dataset's own 31-day-smoothed
  day-of-year climatology.
- D3 uses jcm-strat's detector (`strat_compare.ssw_dates`), ported.
- A2-WACCM: the WACCM6 REF-D1 clock starts at its tropical base point (0.47°N, 103 hPa),
  close to the 150 hPa reset of AOA150, so no offset is subtracted. A2 subtracts CLaMS'
  150 hPa tropical age, which does not put CLaMS on a 150 hPa clock.

---

## 5. References

| Reference | Location | Years | Metrics |
|---|---|---|---|
| CLaMS v3.1 driven by ERA5, zonal-mean AGE | `/data/CLaMS/CLaMS_v3/clams_v3.1_era5_zm_lat.zip` (2004–2023), read in memory | scored window | A2–A4, D7 |
| WACCM6 REF-D1 entry age (base 0.47°N, 103 hPa) | `/data/CESM2_REFD1_AOA/AoA_waccm6_refd1.04_AOA1mf_1970-2019_ba_0_100.0_ck_0_50.0.nc` | scored window | A1, A2-WACCM |
| CESM2.1.5-WACCM6 histSST h6 (`Vzm`, `VTHzm`, `THzm`) | `/data/cesm2.1.5_output/histSST`, through `strat_eval_common.load_h6` | **1996–2014 for every window** (6935 days, 61 interfaces from 4.5 × 10⁻⁶ hPa) | B1, B2 |
| ERA5 monthly zonal means (`uzm`, `tzm`, 25 levels) | `AIDE-atmosphere/data/era5/era5_zm_monthly_<y>.nc` | scored window | C1, C2 |
| ERA5 CDS daily zonal means (u, v, T, θ, v'θ', 17 levels 1–250 hPa) | `AIDE-atmosphere/data/era5/era5_<y>.npz` | scored window | D1–D7 |

- Every path is in `strat_eval_common.REFS`.
- The ERA5 files are copies of jcm-strat's, made by `scripts/copy_era5_refs.sh` (2006–2019,
  28 files) and checked against `scripts/era5_refs.sha256`; `AIDE-atmosphere/data/` is
  gitignored.
- **B is not period-matched.** WACCM6 histSST h6 exists for 1970–1995 and 1996–2014; the
  reference is the whole 1996–2014 run. For 2006–2019 the overlap is 2006–2014 (9 years).
- CESM h6 conventions, established against the files (evidence: the archived scripts
  `01*`, §12): `VTHzm` is the eddy flux, not `Vzm·THzm`; `Wzm` is log-pressure
  w = −Hω/p; below-surface points are the sentinel `1e35` without `_FillValue`, masked with
  `abs(x) < 1e20`; `MSKtem` is a fractional count, not a 0/1 mask.
- Constants are CESM's: a = 6.37122 × 10⁶ m, g = 9.80616 m s⁻², κ = R/c_p = 287.058 /
  1004.64, H = 7 km.

---

## 6. Scoring

jcm-strat `relative_scores.py`:

- **Relative metrics**: rel = model / ref − 1, err = \|rel\| / 0.25. A4 against the midpoint
  2.75 yr of its range. A negative slope keeps its sign in the ratio (−40 % is 40 % too
  small).
- **Absolute metrics** keep their bound: A2 and A2-WACCM < 0.4 yr, C1 latitudes ± 4° on
  \|lat\|, C2 ± 2 K, D4 final-warming dates ± 10 days; err = \|difference\| / bound.
- **Status**: pass err ≤ 1, near err ≤ 1.5, fail above. A metric that cannot be formed fails
  (err 3).
- **Group score** = mean err of the group's scored metrics. **Combined** = mean of the A, B
  and C group scores. D is reported on its own.
- The tolerances are the scorecard's; they are not derived from interannual σ.

### Two scorings

| | **standard** (headline, ranking) | **extensive** (full listing) |
|---|---|---|
| Definition | jcm-strat phase 31's scored set = the former "reduced" scoring | the 34 metrics of the phase-29 scorecard |
| Scored | 20: A 4, B 7, C 6, D 3 | 34: A 4, B 7, C 10, D 13 |
| Age RMSE | A2-WACCM | A2 (shifted CLaMS) |
| Listed, not scored | A2 and phase 31's `make_phase31_report.py` UNSCORED: C1-jet10_N, C1-jet10_N-lat, C1-jet10_S, C1-jet10_S-lat, D1-N, D1-S, D2-N, D2-S, D3, D4-fwN, D4-fwS, D6-N, D6-S, D7 | — (A2-WACCM computed, not shown) |
| Passes counted | A–C of 17, D of 3 | A–C of 21, D of 13 |

- Standard scores the 50 hPa vortex, the residual circulation, the age on the 150 hPa clock,
  D4-w and D5; it leaves out the 10 hPa winds, vortex variability, final warmings, SSW count
  and the two relations D6, D7.
- jcm-strat phase 31 scored its set with absolute bounds (`success_metrics.py`); standard
  applies the relative scoring above to that set. Group scores differ from phase 31's; the
  values behind them agree (§10.3).
- The ranking (`SCORING_STANDARD.md`) orders every run of every set by its standard combined
  score. Runs on different windows are each scored against references of their own years,
  so a cross-window rank compares two differences, not two models on one reference.

---

## 7. w* estimator — interim (AIDE-eval issue #5)

- **Formula** (jcm-strat `strat_circulation.tem_streamfunction`, `residual_w`), from
  time-mean [v], [v'θ'], [θ]:
  v* = [v] − ∂/∂p([v'θ'] / ∂[θ]/∂p); Ψ* = (2πa cos φ / g) ∫ v* dp from the top;
  w* = −H ω*/p with ω* = −(g / (2πa² cos φ)) ∂Ψ*/∂φ. No w̄ needed. \|lat\| > 88° masked.
- **Interpolation**: w* and Ψ* linear in log p to the exact target pressure
  (`interp_level`), not the nearest level. The same for the monthly 70 hPa w* of D4-w and D7.
- **Model and reference** (WACCM6 h6) go through the same functions.
- **Interim**: the estimator is open as AIDE-eval issue #5. A change there changes B, D4-w
  and D7.

What interpolation changes against the nearest level (phase-29 runs, 2006–2010, from
`port_check.md`):

| Metric | WACCM6 nearest → interpolated | Nearest levels model / WACCM6 [hPa] | Model value change |
|---|---|---|---|
| B1-100 | 0.398 → 0.456 mm s⁻¹ (+14.3 %) | 95.9 / 94.9 | +0.7 to +12.8 % |
| B1-70 | 0.211 → 0.228 (+8.0 %) | 72.2 / 66.8 | +2.2 to +5.1 % |
| B1-50 | 0.201 → 0.198 (−1.2 %) | 49.3 / 47.1 | −1.0 to −0.4 % |
| B1-30 | 0.255 → 0.271 (+6.3 %) | 30.4 / 32.6 | −0.6 to +0.2 % |
| B1-10 | 0.469 → 0.463 (−1.1 %) | 9.86 / 9.55 | −1.7 to −0.3 % |
| B2-100 | 10.85 → 12.39 × 10⁹ kg s⁻¹ (+14.1 %) | 95.9 / 94.9 | +5.1 to +6.2 % |
| B2-70 | 6.05 → 6.53 (+7.9 %) | 72.2 / 66.8 | −3.1 to −2.4 % |
| D4-w | ERA5 unchanged (70 hPa is an ERA5 level) | 72.2 / 70 | 0.0 to +8.2 % |
| D7 | ERA5/CLaMS unchanged | 72.2 / 70 | −9.2 to −5.2 % |

---

## 8. Group C daily-mean input (AIDE-eval issue #7)

- C takes calendar-month means of the daily zonal means. jcm-strat `teq_score.py` averages
  the 6-hourly frames by timestamp; jcm-strat `daily_zm.py` builds each day from 06Z, 12Z,
  18Z and 00Z of the next day. Each month therefore loses the 00Z frame that opens it and
  gains the one that opens the next month (one frame in about 120); the ERA5 side is the
  same in both.
- Effect against jcm-strat (port check): C2 differs by ≤ 0.022 K, jets by ≤ 0.37 %; the
  shift estimated from the daily series reproduces every C2 difference within jcm-strat's
  0.01 K rounding. In one case (L0, 2006–2010) it moves the NH 50 hPa jet latitude by one
  grid row (45.70°N vs 47.56°N; the two rows differ by 0.0004 m s⁻¹).
- Open: keep daily-mean input for C (any model can supply it), or require 6-hourly input to
  match jcm-strat exactly.

---

## 9. How to evaluate a model

1. **Write the contract fields** (§3) for 2006-01-01 to 2019-12-31 to netCDF.
2. **Register one result set**: `AIDE-atmosphere/sources/<set>.json`

   ```json
   {
     "description": "MYMODEL v1, free running, 2006-2019 daily zonal means",
     "reader": "netcdf",
     "window": [2006, 2019],
     "runs": {"MYMODEL": "/path/to/mymodel_*.nc"},
     "notes": {"MYMODEL": "one line on the configuration (optional)"}
   }
   ```

   - `runs` maps a run name to a path, a glob, or a list of them (one continuous record).
     Several runs on one window can share a set.
   - `window` other than `[2006, 2019]` is scored and flagged non-standard.
3. **Check it**: from `AIDE-atmosphere/scripts/`, `PY=../../.AIDE-eval_env/bin/python`,
   `$PY 01_fields.py --inspect <set>`.
4. **Evaluate**: `bash copy_era5_refs.sh` once (copies and verifies the ERA5 files), then
   `$PY evaluate.py <set>` → `01_fields` → `02_scorecard` → `03_figures`, and both summaries
   rewritten over every registered set. Heavy runs under the memory cap:
   `bash /data/JCM_stripped/jcm-strat-phase29/scripts/memcap.sh -- $PY evaluate.py <set>`.
   One 14-year single-run set takes about 25 s on CPU, peak 2.6 GB (01_fields).
5. **Read** `validation_results/<set>/scorecard_standard.md` with §7, §8 and §11.

- A free-running model is compared with the observed years climatologically; the jcm-strat
  runs are nudged to ERA5 in the troposphere and in the tropical zonal-mean wind, so for
  them the same-day comparison of C and D is tighter than for a free run.

---

## 10. Results

### 10.1 AL7h, 2006–2019 (standard window)

jcm-strat AL7h (T63 L95, AL7 T_eq correction at half amplitude), the continuous p25_AL7h
chain: phase 29's 2006–2010 file and phase 31's 2011–2019 continuation read as one record.
**Its T_eq correction was fitted on 2005-04..2008-01, inside the window.** From
`validation_results/jcm-strat-AL7h-2006-2019/scorecard_<scoring>.json`:

| | standard | extensive |
|---|---|---|
| A age of air | 0.72 | 0.93 |
| B residual circulation | 0.95 | 0.95 |
| C polar vortex | 0.39 | 0.36 |
| **combined A–C** | **0.69** | **0.74** |
| D variability | 0.12 | 1.10 |
| passed A–C | 15 / 17 | 18 / 21 |
| passed D | 3 / 3 | 6 / 13 |

- Standard, fail: B1-100 (−69 %), B1-50 (+46 %). Every other scored metric passes; A2-WACCM
  0.32 yr.
- Extensive, also fail: A2 (0.65 yr), D1-N (−52 %), D1-S (−51 %), D6-N (−39 %), D6-S
  (−42 %), D7 (−68 %); near: D2-N (−27 %), D2-S (−33 %).

### 10.2 Seven phase-29 runs, 2006–2010 (non-standard window, 5 of 14 standard years)

P20b, AL7, AL7t, L0, LAL7, QBO10d, LQBO10d (QBO30 and QBO10dL76 excluded; AL7h scored in
§10.1). Only 2006–2010 exists for these runs.

Standard:

| | P20b | AL7 | AL7t | L0 | LAL7 | QBO10d | LQBO10d |
|---|---|---|---|---|---|---|---|
| A | 0.56 | 0.95 | 0.98 | 0.79 | 1.47 | 0.71 | 1.20 |
| B | 1.05 | 1.00 | 1.01 | 0.98 | 1.06 | 0.86 | 0.78 |
| C | 1.07 | 0.30 | 0.13 | 1.03 | 0.45 | 0.63 | 0.63 |
| **combined A–C** | **0.89** | **0.75** | **0.71** | **0.93** | **0.99** | **0.73** | **0.87** |
| D | 0.39 | 0.42 | 0.53 | 0.54 | 0.53 | 0.28 | 0.73 |
| passed A–C | 12 / 17 | 13 / 17 | 11 / 17 | 12 / 17 | 11 / 17 | 14 / 17 | 10 / 17 |

Extensive:

| | P20b | AL7 | AL7t | L0 | LAL7 | QBO10d | LQBO10d |
|---|---|---|---|---|---|---|---|
| A | 0.58 | 1.17 | 1.21 | 0.99 | 1.70 | 0.91 | 1.44 |
| B | 1.05 | 1.00 | 1.01 | 0.98 | 1.06 | 0.86 | 0.78 |
| C | 1.09 | 0.26 | 0.21 | 1.02 | 0.45 | 0.79 | 0.47 |
| **combined A–C** | **0.91** | **0.81** | **0.81** | **1.00** | **1.07** | **0.85** | **0.90** |
| D | 1.35 | 1.28 | 1.29 | 1.17 | 1.45 | 0.93 | 1.39 |
| passed A–C | 13 / 21 | 17 / 21 | 15 / 21 | 14 / 21 | 15 / 21 | 15 / 21 | 14 / 21 |

- Standard: B1-100 fails in every run (−50 to −72 %); B1-70, B1-10, B2-100, B2-70, both
  50 hPa jets, the SH 50 hPa jet latitude, D5-N and D5-S pass in every run.
- AL7h on 2006–2010 alone (PR #6 results, commit 34d7d6f): standard combined 0.72,
  extensive 0.79.

> **Interpretation** — On the standard ranking AL7h 2006–2019 (0.69) leads, followed by AL7t
> (0.71), QBO10d (0.73) and AL7 (0.75) on 2006–2010. The first is on 14 years and its own
> references, the others on 5; the 0.02–0.06 gaps are within what a change of window moved
> AL7h's own score (0.72 on 2006–2010, 0.66 on 2011–2019 in the 2026-11-11 results). The
> lower-branch deficit (B1-100, −50 to −72 %) and the excess at 50 hPa (+32 to +63 %) are
> common to every run.

### 10.3 Port check

`scripts/checks/port_check.py` → [validation_results/port_check.md](../validation_results/port_check.md).
Its runs and jcm-strat paths are defined in the script, independent of `sources/`.
Agreement: model and reference value each within max(0.5 % of jcm-strat's value, jcm-strat's
rounding); B, D4-w, D7 at the nearest level.

| Part | Runs, window, age references | A | B | C | D |
|---|---|---|---|---|---|
| 1 — phase 29 | eight runs, 2006–2010, 2005–2009 | 32 / 32 | 56 / 56 | 66 / 80 | 104 / 104 |
| 2 — phase 31 | AL7h 2006–2010, window years | 5 / 5 | 7 / 7 | 8 / 10 | 13 / 13 |
| 2 — phase 31 | AL7h 2011–2019, window years | 5 / 5 | 7 / 7 | 9 / 10 | 13 / 13 |

- A, B, D: largest difference 0.008 % (CESM a, g and κ against jcm-strat's 6.371e6 m,
  9.80665 m s⁻², 0.2857). A2-WACCM reproduces phase 31's `a2_waccm.json` to 1e-5 %.
- C: every disagreeing cell is the one-frame-per-month shift of §8.
- Fed the nearest-level values, the extensive group scores reproduce jcm-strat's
  `relative_scores.json` for the eight runs except C (L0 0.98 → 1.02, others ± 0.01), same
  cause.
- The check runs on the two jcm-strat files AL7h is scored from, separately; the
  2006–2019 record is their concatenation (§3).

---

## 11. Limitations

- **Short windows.** D1 is a standard deviation over 4–13 winters and D3 a count of a few
  events; their sampling error is comparable to the 25 % tolerance. The 2006–2010 set
  carries it more than the 2006–2019 one.
- **AL7h is not out of sample** on 2006–2019: its T_eq correction was fitted on
  2005-04..2008-01.
- **B is not period-matched** and uses its own estimator (§5, §7).
- **Nearest levels** in A (55 hPa), C (50, 10 hPa) and D (10, 50, 55, 100 hPa) make those
  metrics depend on the model's vertical grid. Jet latitudes are grid latitudes (1.87° on
  the jcm-strat grid).
- **References disagree with each other.** A4's range spans CLaMS 3.0 and WACCM6 2.5 yr; A1
  takes WACCM6, A2 and A3 CLaMS; shifted CLaMS and the WACCM6 entry age differ by about
  0.5 yr RMSE on the A2 domain, more than the 0.4 yr bound.
- **Standard ≠ phase 31's score.** Relative scoring on phase 31's set; phase 31's
  absolute-bound group scores are not reproduced.
- **Cross-window ranking** (§6) compares runs against different reference years.

---

## 12. History

- **Until 2026-11-12** the repository's protocol was a two-tier CESM-anchored evaluation
  (tier 1, 5-year screening; tier 2, 35-year validation) with thresholds derived from
  CESM2.1.5-WACCM6 1970–2014. Its scripts, results (CESM, ERA5, ERA-Interim, JRA-55,
  MERRA-2, jcm-strat phase 10) and this document's previous version are archived as git tag
  **`tier-protocol-final`** (commit 34d7d6f).
- 2026-11-10: transport scorecard ported from jcm-strat phase 29 beside the tiers (PR #6).
- 2026-11-11: reference years follow the window; extensive and reduced scorings; AL7h
  2011–2019.
- 2026-11-12: the scorecard replaces the tier protocol; "reduced" renamed **standard**;
  standard window 2006–2019; AL7h scored 2006–2019; one JSON per result set and
  `evaluate.py`.
