# Transport scorecard

A separate scorecard for stratospheric tracer transport — age of air (A), residual
circulation (B), polar vortex (C) and variability (D) — kept beside the two CESM tiers of
[EVALUATION_PROTOCOL.md](EVALUATION_PROTOCOL.md). It does not change any threshold or rule of
that protocol. It is the `jcm-strat` relative scorecard of phase 29, ported into this
repository: the diagnostics are recomputed here from daily zonal means, not copied as numbers.

| | |
|---|---|
| Scored so far | the eight jcm-strat phase-29 runs, 2006–2010, and AL7h 2011–2019 (jcm-strat phase 31, out of sample), each under two scorings → [validation_results/transport/](../validation_results/transport/) |
| Scorings | **extensive** (the 34 phase-29 metrics) and **reduced** (jcm-strat phase 31's scored set) — §4 |
| One-page summaries | [SCORING_EXTENSIVE.md](../validation_results/transport/SCORING_EXTENSIVE.md), [SCORING_REDUCED.md](../validation_results/transport/SCORING_REDUCED.md) — every set, one scoring each |
| Scripts | `scripts/transport_common.py`, `22_transport_fields.py`, `23_transport_scorecard.py`, `24_transport_figures.py`; evidence `22b_transport_port_check.py`; `copy_era5_refs.sh` |
| Agreement with jcm-strat | [port_check.md](../validation_results/transport/port_check.md), summarised in §7 |
| w* estimator | **interim** — AIDE-eval issue #5 is open (§5) |

---

## 1. Field contract

What a model has to supply. This extends the five fields of protocol §4.2 with the age tracer;
`Wzm` is optional here because no metric of this scorecard reads it.

| Name | Units | Definition | Metrics that read it |
|---|---|---|---|
| `Uzm` | m s⁻¹ | zonal-mean zonal wind | C1, D1–D4 |
| `Vzm` | m s⁻¹ | zonal-mean meridional wind | B1, B2, D4-w, D7 |
| `THzm` | K | zonal-mean potential temperature, reference 1000 hPa | B1, B2, D4-w, D7; C2 and D6 through T |
| `VTHzm` | K m s⁻¹ | zonal-mean **eddy** heat flux v'θ' | B1, B2, D4-w, D5, D6, D7 |
| `AOA150` | years | zonal-mean entry age: an age clock reset to zero wherever p > 150 hPa | A1–A4, A2-WACCM, D7 |
| `Wzm` | m s⁻¹ | log-pressure w̄ (protocol §4.2) — optional | none |
| `T` | K | zonal-mean temperature — optional; if absent, `T = THzm (p/1000)^κ` | C2, D6 |

- **Shape** `(time, plev, lat)`, daily means, every field on one common set of levels and
  latitudes.
- **`plev` in hPa**, any order; the reader sorts it ascending. Pressure must be known per level,
  independent of surface pressure, as in protocol §4.4.
- **`lat` in degrees**, global (−90° to +90°), any order; the reader sorts it ascending.
- **`time`** daily, each value labelled by the calendar day of its own mean. CF time with any
  calendar; a `noleap` calendar is read through `cftime`.
- **`VTHzm` is the eddy flux**, not the product of the zonal means — protocol §4.2 and
  convention 1 of [CLAUDE.md](../../CLAUDE.md). It must be formed from instantaneous fields
  before any time averaging; the jcm-strat runs form it within each 6-hourly snapshot and
  average four snapshots to a day.
- **`AOA150` is in years** (365.25 days). The jcm-strat files store the clock in days; their
  reader converts.
- **Levels**: the residual streamfunction Ψ* is integrated from the top supplied level down, so
  the levels have to reach well above 10 hPa (the jcm-strat runs reach 0.01 hPa, CESM h6
  4.5 × 10⁻⁶ hPa). The metrics read 5–150 hPa.
- **Record**: each source carries the window it is scored over (`window` in its `SOURCES`
  entry); every reference except B follows that window (§3). Scored so far: 2006–2010 (1826
  days) and 2011–2019 (3287 days).
- Fill values must be masked before any reduction (protocol §5, rule 4).

Two readers deliver the contract, registered in `transport_common.READERS` in the form of
`20_series_from_zonal_mean.py` (`reader(spec, years)`, one `SOURCES` entry per source):

| Reader | Reads | Mapping |
|---|---|---|
| `jcm_npz` | one jcm-strat `daily_zm.py` npz per run (`days`, `p_hpa`, `lat`, `(day, lev, lat)` float32) | `u`→`Uzm`, `v`→`Vzm`, `th`→`THzm`, `vth`→`VTHzm`, `aoa150`/365.25→`AOA150`, `T` kept |
| `netcdf` | one or more netCDF files (glob) carrying the contract names, dims `(time, plev, lat)` | none; values beyond 1e20 masked; `T` derived from `THzm` if absent |

---

## 2. Metrics

Window means over the scored window unless stated. "Nearest level" is the model level nearest
the target pressure, as jcm-strat takes it; only w* and Ψ* are interpolated (§5). Reference
values are given for 2006–2010 and 2011–2019 where they differ.

| ID | Quantity | Reference, 2006–2010 / 2011–2019 | Scored |
|---|---|---|---|
| **A1** | entry age 10S–10N (cos-weighted), nearest level to 55 hPa | WACCM6 REF-D1 entry age, same reduction: 1.099 / 1.105 yr | relative |
| **A2** | cos-weighted RMSE of AOA150 against CLaMS minus CLaMS' own 10S–10N age at 150 hPa (0.091 / 0.089 yr), model put on the CLaMS grid, 100–5 hPa, \|lat\| ≤ 80° | 0 | absolute, < 0.4 yr; extensive only |
| **A2-WACCM** | cos-weighted RMSE of AOA150 against the WACCM6 REF-D1 2-D entry age, model put on the WACCM6 grid, 100–5 hPa, \|lat\| ≤ 80°, no offset (jcm-strat phase 31 `a2_waccm_age.py`) | 0 | absolute, < 0.4 yr; reduced only |
| **A3** | age step from 10S–10N to 50–70° (both hemispheres), nearest level to 55 hPa | CLaMS step: 2.84 / 2.82 yr | relative |
| **A4** | transit 70 → 10 hPa: difference of the 10S–10N age profile, linear in log p | range 2.3–3.2 yr (CLaMS 3.0, WACCM6 2.5) | relative to the midpoint 2.75 yr |
| **B1-p** | w* 15S–15N (cos-weighted) at p = 100, 70, 50, 30, 10 hPa | WACCM6 h6 1996–2014, same estimator: 0.456, 0.228, 0.198, 0.271, 0.463 mm s⁻¹ (both windows) | relative |
| **B2-p** | tropical upward mass flux: max − min of Ψ* over \|lat\| ≤ 60° at p = 100, 70 hPa | WACCM6: 12.39, 6.53 × 10⁹ kg s⁻¹ (both windows) | relative |
| **C1-jetP_H** | winter (DJF north, JJA south) maximum of the monthly-mean u poleward of 30°, nearest level to P = 50, 10 hPa | ERA5 monthly zonal means of the same months on the model's grid | relative |
| **C1-…-lat** | grid latitude of that maximum | ERA5, same | absolute, ± 4° on \|lat\| |
| **C2-H** | winter polar-cap (60–90°, cos-weighted) T at the nearest level to 50 hPa, minus ERA5 | 0 | absolute, ± 2 K |
| **D1-H** | standard deviation over winters of the winter-mean u at 10 hPa, 55–65° | ERA5 CDS daily, same days: 6.51 / 2.16 and 10.49 / 2.62 m s⁻¹ | relative |
| **D2-H** | standard deviation of daily u anomalies at 60°, 10 hPa, within winter | 14.99 / 5.25 and 16.46 / 5.63 m s⁻¹ | relative |
| **D3** | NH major SSWs (u 60°N 10 hPa reversal Nov–Mar, ≥ 20 days apart, not a final warming) | 8 / 4 in ERA5 | relative |
| **D4-fwH** | mean final-warming date: first easterly crossing of the 5-day-smoothed u 60° 10 hPa after which u stays below 5 m s⁻¹ (Mar–Jul north, Sep–Jan south), springs common to model and ERA5 | ERA5: 44.4 / 89.8 and 38.7 / 75.2 days after 1 Mar / 1 Sep | absolute, ± 10 days |
| **D4-w** | annual-cycle amplitude (max − min of the monthly climatology) of w* 15S–15N at 70 hPa | 0.185 / 0.193 mm s⁻¹ | relative |
| **D5-H** | standard deviation of daily anomalies of v'T' at 100 hPa, 45–75°, within winter | 11.47 / 7.90 and 11.80 / 7.04 K m s⁻¹ | relative |
| **D6-H** | slope of the 50 hPa polar-cap T anomaly on the trailing 40-day mean 100 hPa heat-flux anomaly, Jan–Mar north, Jul–Oct south | 0.972 / 0.766 and 1.055 / 1.234 K per K m s⁻¹ | relative |
| **D7** | slope of the deseasonalised, detrended monthly 10S–10N age at 55 hPa on the trailing 3-month w* anomaly at 70 hPa | CLaMS age with ERA5 w*: −4.52 / −3.98 yr per mm s⁻¹ | relative |

- H = N or S; anomalies in D are taken against each dataset's own 31-day-smoothed day-of-year
  climatology, so model and reference lose their own seasonal cycle.
- D3 uses jcm-strat's detector (`strat_compare.ssw_dates`), ported, not
  `aide_val_common.detect_ssw_cp07`; the two implement Charlton–Polvani differently.
- **A2-WACCM** replaces A2 in the reduced scoring. The WACCM6 REF-D1 clock starts at its own
  tropical base point (0.47°N, 103 hPa), close to the model's 150 hPa reset, so no offset is
  subtracted; A2 subtracts CLaMS' 150 hPa tropical age (0.09 yr), which does not put CLaMS on a
  150 hPa clock. AL7h: 0.325 yr (2006–2010), 0.318 yr (2011–2019); bias −0.065 / −0.082 yr.

---

## 3. References

| Reference | Location | Years | Metrics |
|---|---|---|---|
| CLaMS v3.1 driven by ERA5, zonal-mean AGE | `/data/CLaMS/CLaMS_v3/clams_v3.1_era5_zm_lat.zip` (2004–2023), read in memory | the scored window | A2–A4, D7 |
| WACCM6 REF-D1 entry age (base 0.47°N, 103 hPa) | `/data/CESM2_REFD1_AOA/AoA_waccm6_refd1.04_AOA1mf_1970-2019_ba_0_100.0_ck_0_50.0.nc` (1970–2019) | the scored window | A1, A2-WACCM, figure |
| CESM2.1.5-WACCM6 histSST h6 (`Vzm`, `VTHzm`, `THzm`) | `aide_val_common.load_h6("1996-2014", …)`, interfaces from the model top to 652 hPa, 1e35 masked | **1996–2014 for every window** | B1, B2 |
| ERA5 monthly zonal means (`uzm`, `tzm`, 25 levels) | `AIDE-atmosphere/data/era5/era5_zm_monthly_<y>.nc` | the scored window | C1, C2 |
| ERA5 CDS daily zonal means (u, v, T, θ, v'θ', 17 levels 1–250 hPa, 2.5°) | `AIDE-atmosphere/data/era5/era5_<y>.npz` | the scored window | D1–D7 |

Reference years per window:

| Window | A (CLaMS, WACCM6 entry age) | B (WACCM6 h6 TEM) | C (ERA5 monthly) | D (ERA5 CDS daily; CLaMS for D7) |
|---|---|---|---|---|
| 2006–2010 (phase-29 set) | 2006–2010 | 1996–2014 | 2006–2010 | 2006–2010 |
| 2011–2019 (AL7h, phase 31) | 2011–2019 | 1996–2014 | 2011–2019 | 2011–2019 |

- Every path is in one block, `transport_common.REFS`; the years follow each source's `window`.
- The ERA5 files are copies made by `scripts/copy_era5_refs.sh` from jcm-strat's cache, only
  the 28 a metric reads (both kinds, 2006–2019); `scripts/era5_refs.sha256` is their manifest,
  and the script checks it. They are identical to the files jcm-strat phases 29 and 31 read
  (phase 31's `cache/era5_ref` and `cache/variability_cds` are symlinks to the copied
  directories). `AIDE-atmosphere/data/` is gitignored.
- **Period matching.** A, C, D and D7 use the scored window's own years. The phase-29 set moved
  from jcm-strat's age reference years 2005–2009 to 2006–2010 on 2026-11-11; the effect is in
  §7. B uses the whole 1996–2014 segment: WACCM6 histSST h6 has no other segment and ends in
  2014, so B is not period-matched for either window (protocol §5, rule 3, binds on the tiers,
  not here).

---

## 4. Scoring and tolerances

jcm-strat `relative_scores.py`, unchanged:

- **Relative metrics**: rel = model / ref − 1 (−40 % is 40 % too small, also for a negative
  slope), err = \|rel\| / 0.25. A4 is scored against the midpoint 2.75 yr of its range.
- **Absolute metrics** keep their bound: A2 and A2-WACCM (< 0.4 yr), C1 latitudes (± 4°, on
  \|lat\|), C2 (± 2 K), D4 final-warming dates (± 10 days); err = \|difference\| / bound.
- **Status**: pass err ≤ 1, near err ≤ 1.5, fail above. A metric that cannot be formed fails
  (err 3).
- **Group score** = mean err of the group's scored metrics. **Combined** = mean of the A, B and C
  group scores; D is reported on its own, so a combined score stays comparable with jcm-strat's
  earlier phases.
- These tolerances are the scorecard's. They are not derived from CESM interannual σ as the
  tier thresholds are, and they do not replace or modify any tier verdict.

### Two scorings

Every result set is scored twice (`transport_common.SCORINGS`); the diagnostics are the same.

| | extensive | reduced |
|---|---|---|
| Metrics scored | the 34 of the phase-29 scorecard: A 4, B 7, C 10, D 13 | 20: A 4, B 7, C 6, D 3 |
| Age RMSE | A2 (shifted CLaMS) | A2-WACCM (WACCM6 entry age) |
| Listed, not scored | — | A2 and jcm-strat phase 31's `make_phase31_report.py` UNSCORED: C1-jet10_N, C1-jet10_N-lat, C1-jet10_S, C1-jet10_S-lat, D1-N, D1-S, D2-N, D2-S, D3, D4-fwN, D4-fwS, D6-N, D6-S, D7 |
| Passes counted | A–C of 21, D of 13 | A–C of 17, D of 3 |

- The reduced set leaves out the winds at 10 hPa (jets, vortex variability, final warmings),
  the SSW count and the two relations D6 and D7; it scores the 50 hPa vortex, the residual
  circulation, the age on the 150 hPa clock, D4-w and D5.
- **jcm-strat phase 31 scored its set with absolute bounds** (`success_metrics.py`); the reduced
  scoring here applies the relative scoring above to the same set. Its group scores therefore
  differ from phase 31's; the values behind them agree (§7).
- Unscored metrics stay in the reports ("Listed, not scored") and are faded in the heatmaps.

---

## 5. The w* estimator — interim

- **Formula.** jcm-strat `strat_circulation.tem_streamfunction` and `residual_w`, ported into
  `transport_common.py`: from the time-mean [v], [v'θ'] and [θ],
  v* = [v] − ∂/∂p([v'θ'] / ∂[θ]/∂p),
  Ψ* = (2πa cos φ / g) ∫ v* dp from the top level,
  w* = −H ω*/p with ω* = −(g / (2πa² cos φ)) ∂Ψ*/∂φ. No w̄ is needed. Rows poleward of 88° are
  masked.
- **Interpolation.** w* and Ψ* are interpolated linearly in log p to the exact target pressure
  with `aide_val_common.interp_level`, not taken at the nearest level as jcm-strat does. The
  same applies to the monthly 70 hPa w* of D4-w and D7, so every w* in the scorecard comes from
  one estimator.
- **Model and reference travel the same code**: WACCM6 h6 1996–2014 is reduced by the same
  functions with the same interpolation.
- **Not the tier estimator.** The protocol pins `aide_val_common.tem_residual` (§5, rule 1)
  for the tiers; this scorecard pins its own. B numbers here are not comparable with the tier
  mass flux or w̄*, and must not be read against their thresholds.
- **Interim.** The choice is open as AIDE-eval issue #5 and will change if that issue settles
  on another estimator.

What the interpolation changes, against the nearest level (from `port_check.md`, phase-29 set):

| Metric | WACCM6 nearest → interpolated | Nearest levels (model / WACCM6) [hPa] | Change of the model value, 8 runs |
|---|---|---|---|
| B1-100 | 0.398 → 0.456 mm s⁻¹ (+14.3 %) | 95.9 / 94.9 | +0.7 to +12.8 % |
| B1-70 | 0.211 → 0.228 (+8.0 %) | 72.2 / 66.8 | +2.2 to +5.1 % |
| B1-50 | 0.201 → 0.198 (−1.2 %) | 49.3 / 47.1 | −1.0 to −0.4 % |
| B1-30 | 0.255 → 0.271 (+6.3 %) | 30.4 / 32.6 | −0.6 to +0.2 % |
| B1-10 | 0.469 → 0.463 (−1.1 %) | 9.86 / 9.55 | −1.7 to −0.3 % |
| B2-100 | 10.85 → 12.39 × 10⁹ kg s⁻¹ (+14.1 %) | 95.9 / 94.9 | +5.1 to +6.2 % |
| B2-70 | 6.05 → 6.53 (+7.9 %) | 72.2 / 66.8 | −3.1 to −2.4 % |
| D4-w | 0.185 → 0.185 (ERA5 has a 70 hPa level) | 72.2 / 70 | 0.0 to +8.2 % |
| D7 | −4.52 → −4.52 (same) | 72.2 / 70 | −9.2 to −5.2 % |

---

## 6. Running it

From `AIDE-atmosphere/scripts/`, `PY=../../.AIDE-eval_env/bin/python`. Heavy steps under the
jcm-strat memory cap (`memcap.sh`, 96 G); the whole chain takes about a minute on CPU (22 peaks
at 2.6 GB for the 2011–2019 source).

```bash
bash copy_era5_refs.sh                     # ERA5 copies into ../data/era5/, checked against the manifest
$PY 22_transport_fields.py                 # every registered source; or name them
                                           #   -> output/22_transport__references__<Y0>-<Y1>.json (one per window)
                                           #   -> output/22_transport__<run>.json
$PY 22b_transport_port_check.py            # evidence only -> validation_results/transport/port_check.md
$PY 23_transport_scorecard.py              # every set, both scorings; or: SET ... [--scoring extensive|reduced]
                                           #   -> transport_scorecard__<set>__<scoring>__<date>.md / .json
                                           #   -> SCORING_EXTENSIVE.md, SCORING_REDUCED.md (one page per
                                           #      scoring over every set; only when every set is scored)
$PY 24_transport_figures.py                # every scorecard of the newest date [--date YYYY-MM-DD]
```

- Result sets are `transport_common.SETS`: `jcm-strat-phase29-2006-2010` (eight runs) and
  `jcm-strat-phase31-AL7h-2011-2019` (source `AL7h_1119`). `transport_common.COMPARE` adds the
  AL7h in-sample / out-of-sample comparison to the second.
- Figures: score table and heatmap per set and scoring (`__<set>__<scoring>__<date>`); w* and
  age profiles per set, and the window comparison, stamped `__<set>__<date>`.
- `output/` is gitignored; `validation_results/transport/` is generated and committed. Never
  edit either by hand.
- Re-running the chain (copy, 22, 22b, 23, 24) from an empty `output/` reproduces every
  committed artefact bit for bit under the pinned `requirements.txt` (checked 2026-11-11 and
  2026-11-12, two runs each).

---

## 7. Port check

`22b_transport_port_check.py` compares every metric with jcm-strat's own numbers. Agreement:
model and reference value both within 0.5 % of jcm-strat's, or within the rounding jcm-strat
applied before scoring. B, D4-w and D7 are checked at the nearest level.

**Part 1 — phase 29, eight runs, 2006–2010**, against `relative_metrics.json` and
`success_metrics.json`, group A at phase 29's age reference years 2005–2009 (recomputed in 22b):

| Group | Cells agreeing | Largest model-value difference | Cause of what remains |
|---|---|---|---|
| A | 32 / 32 | < 0.001 % | — |
| B (nearest level) | 56 / 56 | 0.008 % | CESM a and g instead of jcm-strat's 6.371e6 m and 9.80665 m s⁻² |
| C | 66 / 80 | C2: ≤ 0.022 K; jets ≤ 0.37 % | one 6-hourly frame per month (below) |
| D | 104 / 104 | 0.008 % | κ = R/c_p instead of 0.2857 in v'T' |

**Part 2 — phase 31, AL7h**, age references of the window's years, against
`31_al7h_validation/success_<0610|1119>/success_metrics.json`, `a2_waccm.json` and group D
(phase 31 `variability.json` for 2011–2019; phase 30's for 2006–2010, as phase 31 reports it):

| Window | A (incl. A2-WACCM) | B | C | D |
|---|---|---|---|---|
| 2006–2010 | 5 / 5 | 7 / 7 | 8 / 10 | 13 / 13 |
| 2011–2019 | 5 / 5 | 7 / 7 | 9 / 10 | 13 / 13 |

- A2-WACCM reproduces `a2_waccm.json` to 1e-5 % in both windows; A, B, D differ by ≤ 0.008 %
  for the causes of part 1.
- **C (AIDE-eval issue #7).** jcm-strat's `teq_score.py` takes calendar-month means of the
  6-hourly frames by timestamp; the port takes them from the daily means, whose days run
  06Z–00Z (`daily_zm.py`). Each month loses its opening 00Z frame and gains the next month's.
  The shift estimated from the daily series matches every C2 difference within jcm-strat's
  0.01 K rounding: the 16 of part 1, and the three of part 2 (AL7h 2006–2010 C2-N +0.019 K
  observed / +0.018 K expected, C2-S −0.008 / −0.005 K; 2011–2019 C2-N +0.020 / +0.020 K). In
  part 1 it also moves L0's 50 hPa NH jet by one grid row (45.70°N instead of 47.56°N; the two
  rows differ by 0.0004 m s⁻¹). All disagreeing cells are these. Whether to keep daily-mean
  input for C is open as issue #7.
- Fed the nearest-level values, `23`'s extensive group scores reproduce jcm-strat's
  `relative_scores.json` for all eight runs, except C (L0 0.98 → 1.02, others ± 0.01) for the
  reason above.

**Part 3 — age references 2005–2009 → 2006–2010** (phase-29 set; only group A reads them):

| Reference | 2005–2009 | 2006–2010 |
|---|---|---|
| A1 WACCM6 entry age 55 hPa | 1.107 yr | 1.099 yr |
| A3 CLaMS age step | 2.793 yr | 2.839 yr |
| A2 CLaMS 150 hPa offset | 0.092 yr | 0.091 yr |

| Run | A extensive | combined extensive | A reduced | combined reduced |
|---|---|---|---|---|
| P20b | 0.567 → 0.577 | 0.903 → 0.906 | 0.536 → 0.557 | 0.888 → 0.895 |
| AL7 | 1.176 → 1.173 | 0.814 → 0.813 | 0.938 → 0.947 | 0.746 → 0.749 |
| AL7h | 0.894 → 0.890 | 0.789 → 0.788 | 0.707 → 0.716 | 0.712 → 0.715 |
| AL7t | 1.213 → 1.210 | 0.810 → 0.809 | 0.975 → 0.984 | 0.703 → 0.706 |
| L0 | 0.987 → 0.986 | 0.996 → 0.995 | 0.778 → 0.789 | 0.929 → 0.933 |
| LAL7 | 1.704 → 1.703 | 1.072 → 1.072 | 1.459 → 1.469 | 0.988 → 0.992 |
| QBO10d | 0.911 → 0.908 | 0.853 → 0.852 | 0.704 → 0.714 | 0.731 → 0.734 |
| LQBO10d | 1.437 → 1.437 | 0.898 → 0.898 | 1.188 → 1.200 | 0.865 → 0.869 |

- Largest change of any group A or combined score: 0.021 (P20b, A reduced). A1 moves by
  +0.4 to +0.7 % points, A3 by −1.1 to −1.4, A4 not at all; A2 falls by 0.013–0.018 yr, A2-WACCM
  changes by −0.001 to +0.005 yr. The one status change is QBO10d's A2 (fail → near). P20b's A2
  err goes 1.05 → 1.002 and stays near (PR #6 rounded this to "near → pass").
- Per-metric values for every run are in part 3 of `port_check.md`.

---

## 8. Result, 2006–2010

Transcribed from the two reports
[extensive](../validation_results/transport/transport_scorecard__jcm-strat-phase29-2006-2010__extensive__2026-11-12.md)
and
[reduced](../validation_results/transport/transport_scorecard__jcm-strat-phase29-2006-2010__reduced__2026-11-12.md);
the metric tables and figures are beside them.

Extensive:

| | P20b | AL7 | AL7h | AL7t | L0 | LAL7 | QBO10d | LQBO10d |
|---|---|---|---|---|---|---|---|---|
| A age of air | 0.58 | 1.17 | 0.89 | 1.21 | 0.99 | 1.70 | 0.91 | 1.44 |
| B residual circulation | 1.05 | 1.00 | 0.96 | 1.01 | 0.98 | 1.06 | 0.86 | 0.78 |
| C polar vortex | 1.09 | 0.26 | 0.52 | 0.21 | 1.02 | 0.45 | 0.79 | 0.47 |
| **combined A–C** | **0.91** | **0.81** | **0.79** | **0.81** | **1.00** | **1.07** | **0.85** | **0.90** |
| D variability | 1.35 | 1.28 | 1.11 | 1.29 | 1.17 | 1.45 | 0.93 | 1.39 |
| passed A–C | 13 / 21 | 17 / 21 | 17 / 21 | 15 / 21 | 14 / 21 | 15 / 21 | 15 / 21 | 14 / 21 |
| passed D | 4 / 13 | 6 / 13 | 6 / 13 | 5 / 13 | 6 / 13 | 5 / 13 | 7 / 13 | 4 / 13 |

Reduced:

| | P20b | AL7 | AL7h | AL7t | L0 | LAL7 | QBO10d | LQBO10d |
|---|---|---|---|---|---|---|---|---|
| A age of air | 0.56 | 0.95 | 0.72 | 0.98 | 0.79 | 1.47 | 0.71 | 1.20 |
| B residual circulation | 1.05 | 1.00 | 0.96 | 1.01 | 0.98 | 1.06 | 0.86 | 0.78 |
| C polar vortex | 1.07 | 0.30 | 0.47 | 0.13 | 1.03 | 0.45 | 0.63 | 0.63 |
| **combined A–C** | **0.89** | **0.75** | **0.72** | **0.71** | **0.93** | **0.99** | **0.73** | **0.87** |
| D variability | 0.39 | 0.42 | 0.32 | 0.53 | 0.54 | 0.53 | 0.28 | 0.73 |
| passed A–C | 12 / 17 | 13 / 17 | 15 / 17 | 11 / 17 | 12 / 17 | 11 / 17 | 14 / 17 | 10 / 17 |
| passed D | 3 / 3 | 3 / 3 | 3 / 3 | 3 / 3 | 3 / 3 | 3 / 3 | 3 / 3 | 2 / 3 |

- Extensive scores against PR #6 (2005–2009 age references): P20b A 0.57 → 0.58 and combined
  0.90 → 0.91, AL7 A 1.18 → 1.17 (1.176 → 1.173); every other group score is unchanged at two
  decimals.
- B1-100 fails in every run (−50 to −72 %), B1-50 fails or is near in every run (+32 to
  +63 %), B2-70 passes in every run.
- Ranking on the combined score: extensive AL7h (0.79), AL7 and AL7t (0.81) first, LAL7 (1.07)
  last; reduced AL7t (0.71), AL7h (0.72), QBO10d (0.73) first, LAL7 (0.99) last. The reduced
  scoring removes the 10 hPa jets and most of D, where every run is weak (D: 0.93–1.45
  extensive, 0.28–0.73 reduced).

> **Interpretation** — The interpolated WACCM6 reference at 100 hPa is 14 % higher than the
> nearest-level one, so the lower-branch deficit jcm-strat reported (w* at 100 hPa about half
> of WACCM6's) is larger here, not smaller. The three best runs are within 0.02 of each other in
> either scoring; their order is not resolved by this scorecard.

---

## 9. AL7h out of sample, 2011–2019

AL7h's T_eq correction was fitted to ERA5 extratropical T over 2005-04..2008-01. jcm-strat
phase 31 continued the same AL7h chain to 2019; 2011–2019 is out of sample, 2006–2010 overlaps
the fit. Transcribed from the comparison section of
[transport_scorecard__jcm-strat-phase31-AL7h-2011-2019__extensive__2026-11-12.md](../validation_results/transport/transport_scorecard__jcm-strat-phase31-AL7h-2011-2019__extensive__2026-11-12.md)
(the same table is in the reduced report) and its figure
`transport_window_comparison__jcm-strat-phase31-AL7h-2011-2019__2026-11-12.png`.

| | extensive 2006–2010 | extensive 2011–2019 | reduced 2006–2010 | reduced 2011–2019 |
|---|---|---|---|---|
| A age of air | 0.89 | 0.95 | 0.72 | 0.73 |
| B residual circulation | 0.96 | 0.94 | 0.96 | 0.94 |
| C polar vortex | 0.52 | 0.30 | 0.47 | 0.30 |
| **combined A–C** | **0.79** | **0.73** | **0.72** | **0.66** |
| D variability | 1.11 | 1.28 | 0.32 | 0.29 |
| passed A–C | 17 / 21 | 18 / 21 | 15 / 17 | 15 / 17 |
| passed D | 6 / 13 | 5 / 13 | 3 / 3 | 3 / 3 |

- Scored in the reduced set, no metric changes status between the windows; B1-100 (−68 → −69 %)
  and B1-50 (+46 → +45 %) fail in both. The largest changes of a reduced metric are D4-w (+6 →
  −20 %), D5-N (+14 → −0.5 %), C1-jet50_N (−13 → −8 %) and the NH 50 hPa jet latitude (−3.7° → 0°).
- In the metrics only the extensive scoring counts, 2011–2019 is further from ERA5/CLaMS in D1-N
  (−46 → −57 %), D1-S (−29 → −46 %), D6-S (−26 → −44 %) and D7 (−61 → −75 %), closer in
  C1-jet10_N (−27 → −17 %). D3 goes from 3 SSWs against ERA5's 8 to 6 against 4 (−62 → +50 %).
- Every reference except B changes with the window, so a change in % off the reference mixes
  the model's change and the reference's (e.g. ERA5 D1-N 6.5 → 10.5 m s⁻¹). Nine winters in
  2011–2019 against four or five in 2006–2010 also change the sampling error of D1 and D3.

> **Interpretation** — Out of sample, AL7h's reduced group scores stay within 0.03 of in sample
> except C, which improves (0.47 → 0.30); the combined score is lower (better) in 2011–2019 in
> both scorings, driven by C. On these metrics the T_eq correction fitted to 2005–2008 shows no
> loss on the 50 hPa vortex in 2011–2019 (C2-S +1.26 → +1.42 K, inside ± 2 K). Extensive A
> worsens by 0.06 (A2 against shifted CLaMS, A3) and extensive D by 0.17: the vortex variability
> and the two relations (D1, D2, D6, D7), all outside the reduced set, are weak in both windows
> and weaker in 2011–2019.

---

## 10. How to score a new model

1. **Write the contract fields** (§1) as daily zonal means on pressure levels to one or more
   netCDF files with dims `(time, plev, lat)` and the names `Uzm`, `Vzm`, `THzm`, `VTHzm`,
   `AOA150` (`T`, `Wzm` optional). Form `VTHzm` from instantaneous 3-D fields before any time
   averaging; give `AOA150` in years, reset to zero wherever p > 150 hPa.
2. **Register it** in `transport_common.SOURCES`:
   `"NAME": dict(reader="netcdf", path="/path/to/files_*.nc", window=(2006, 2010),
   long_name="…", note="…", independent=True, time_resolution="daily")`, and add a result set
   to `transport_common.SETS`: `"SET": ("NAME", …)` (every source of a set shares its window).
3. **Check it**: `$PY 22_transport_fields.py --inspect NAME` prints shapes, units, ranges and
   the number of levels in 5–150 hPa. The levels must reach well above 10 hPa (§1).
4. **Score it**: `bash copy_era5_refs.sh`, then `$PY 22_transport_fields.py NAME`,
   `$PY 23_transport_scorecard.py SET`, `$PY 24_transport_figures.py`.
5. **Read it with its limits**: §11 and the interim w* note of §5.

- The window must lie within 2006–2019: C and D compare the same days of ERA5 and CLaMS, and
  only those years' ERA5 files are copied (2004–2005 and 2020–2023 would need copies; §13).
- A free-running model is compared with the observed years climatologically; jcm-strat's runs
  are nudged to ERA5 in the troposphere and in the tropical zonal-mean wind, so for them the
  same-day comparison is tighter than for a free run.

---

## 11. Limitations

- **Short windows.** D1 is a standard deviation over 4–5 (2006–2010) or 8–9 (2011–2019) winters
  and D3 a count of a few events; their sampling error is comparable to the 25 % tolerance, and
  differs between the two windows.
- **Nearest levels** in A (55 hPa), C (50, 10 hPa) and D (10, 50, 55, 100 hPa) make those
  metrics depend on the model's vertical grid; only w* and Ψ* are interpolated.
- **Jet latitudes** are grid latitudes (1.87° on the jcm-strat grid).
- **B is not period-matched** and uses a different estimator from the tiers (§3, §5); for
  2011–2019 it compares with WACCM6 1996–2014.
- **References disagree with each other** — A4's range spans CLaMS 3.0 and WACCM6 2.5 yr; A1
  takes WACCM6's entry age while A2 and A3 take CLaMS; WACCM6's entry age and shifted CLaMS
  differ by 0.48 yr (2006–2010) and 0.58 yr (2011–2019) RMSE on the A2 domain (phase 31
  `a2_waccm.json`), so A2 and A2-WACCM measure against references that disagree by more than
  the 0.4 yr bound.
- **Reduced ≠ phase 31's score.** The reduced scoring uses relative scoring; phase 31's group
  scores (absolute bounds) are not reproduced and not comparable.

---

## 12. Decisions

Settled by the user before the port (2026-11-10):

- A separate scorecard beside the two CESM tiers; the protocol's thresholds unchanged.
- Port the diagnostics, not the cached numbers; score the eight phase-29 runs over 2006–2010
  (QBO30 and QBO10dL76 excluded).
- w*: the jcm-strat formula, tropical values interpolated linearly in log p to the exact target
  pressures, |lat| > 88° masked — **interim, AIDE-eval issue #5 stays open**. Mass flux: range
  of Ψ* within 60°, interpolated to 100 and 70 hPa. WACCM6 B references recomputed from h6
  1996–2014 with the same method.
- Age: entry age only (`aoa150`).
- The field contract of §1 and the two readers in the `READERS` pattern of PR #4.
- Scoring as jcm-strat `relative_scores.py` (§4).

Settled by the user on 2026-11-11:

- Score AL7h out of sample, 2011–2019 (jcm-strat phase 31 `AL7h_1119.npz`, read in place), beside
  AL7h 2006–2010 in sample.
- Reference years follow the scored window for every set; the phase-29 set moves to 2006–2010
  age references. The WACCM6 w* / Ψ* reference stays h6 1996–2014, the only h6 segment.
- Two scorings of every set: extensive (the 34 metrics, unchanged) and reduced (phase 31's
  scored set, A2-WACCM for A2, relative scoring).

Taken without confirmation — listed in the PR for review:

- Readers and `SOURCES` live in `transport_common.py`, in 20's form, not in 20 itself.
- The `netcdf` reader takes a glob, dims named `time`, `plev`, `lat`, derives `T` from `THzm`
  with κ = R/c_p, masks values beyond 1e20.
- CESM constants (a, g, κ = R/c_p) instead of jcm-strat's; ≤ 0.008 % on any metric (§7).
- WACCM6 h6 read from the model top to 652 hPa with the sentinel masked; jcm-strat read all
  71 interfaces unmasked. No metric level is affected.
- Ψ* is interpolated to the target pressure before its latitude range is taken.
- The D4-w and D7 w* is interpolated too, as one estimator for every w*.
- A, C and D stay at the nearest level, as jcm-strat.
- D3 keeps jcm-strat's SSW detector.
- Monthly means for C from the daily means (§7; issue #7).
- Port-check criterion and a separate evidence script `22b`.
- Stamped artefact names, the scorecard JSON committed beside the report.
- Heatmap colours an absolute metric as a fraction of its bound, the bound at 25 %.
- Each source carries its `window`; result sets in `transport_common.SETS`; the 2011–2019
  source is named `AL7h_1119` (jcm-strat's file stem) in the set `jcm-strat-phase31-AL7h-2011-2019`.
- The extensive scoring shows exactly the 34 metrics (A2-WACCM computed, not shown); the
  reduced scoring lists A2 as unscored, as phase 31's report keeps the old A2 as a faded row.
- The comparison is written into both 2011–2019 reports and drawn once; in sample orange,
  out of sample blue, as phase 31's figures.
- 23 takes set names and `--scoring`; 24 takes `--date`; stamps `__<set>__<scoring>__<date>`.
- 22b keeps the phase-29 check at 2005–2009 by recomputing group A there.

---

## 13. Deferred

- **ERA5 1990–2022 context for D** (`variability_score.py --long`): how much D statistics of a
  5- or 9-year window scatter in ERA5. Not ported; its 33 yearly files were not copied.
- **ERA5 TEM w* context** (jcm-strat `era5_tem_monthly` 1996–2001): shown by jcm-strat beside
  B, never scored. Not ported.
- **CDS download script** (jcm-strat `fetch_era5_daily_zm.py`): the ERA5 files are copies of
  jcm-strat's; regenerating them from the CDS is not part of this repository yet.
- **Windows outside 2006–2019**: need the ERA5 files of those years (and CLaMS covers only
  2004–2023, WACCM6 entry age 1970–2019).
- **Phase 31's absolute-bound scores** and its imposed / tuned / emergent classes of each metric
  (`make_phase31_report.py`): not ported; the reduced scoring takes only its metric set.
- **AL7h 2011–2019 profiles beside 2006–2010**: the profile figures show one window each; the
  window comparison is metric by metric only.
