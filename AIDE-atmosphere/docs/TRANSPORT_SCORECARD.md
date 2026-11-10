# Transport scorecard

A separate scorecard for stratospheric tracer transport — age of air (A), residual
circulation (B), polar vortex (C) and variability (D) — kept beside the two CESM tiers of
[EVALUATION_PROTOCOL.md](EVALUATION_PROTOCOL.md). It does not change any threshold or rule of
that protocol. It is the `jcm-strat` relative scorecard of phase 29, ported into this
repository: the diagnostics are recomputed here from daily zonal means, not copied as numbers.

| | |
|---|---|
| Scored so far | the eight jcm-strat phase-29 runs, 2006–2010 → [validation_results/transport/](../validation_results/transport/) |
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
| `AOA150` | years | zonal-mean entry age: an age clock reset to zero wherever p > 150 hPa | A1–A4, D7 |
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
- **Record**: the window scored here is 2006–2010 (1826 days), because the references are
  matched to it (§3).
- Fill values must be masked before any reduction (protocol §5, rule 4).

Two readers deliver the contract, registered in `transport_common.READERS` in the form of
`20_series_from_zonal_mean.py` (`reader(spec, years)`, one `SOURCES` entry per source):

| Reader | Reads | Mapping |
|---|---|---|
| `jcm_npz` | one jcm-strat `daily_zm.py` npz per run (`days`, `p_hpa`, `lat`, `(day, lev, lat)` float32) | `u`→`Uzm`, `v`→`Vzm`, `th`→`THzm`, `vth`→`VTHzm`, `aoa150`/365.25→`AOA150`, `T` kept |
| `netcdf` | one or more netCDF files (glob) carrying the contract names, dims `(time, plev, lat)` | none; values beyond 1e20 masked; `T` derived from `THzm` if absent |

---

## 2. Metrics

Window means over 2006–2010 unless stated. "Nearest level" is the model level nearest the
target pressure, as jcm-strat takes it; only w* and Ψ* are interpolated (§5).

| ID | Quantity | Reference | Scored |
|---|---|---|---|
| **A1** | entry age 10S–10N (cos-weighted), nearest level to 55 hPa | WACCM6 REF-D1 entry age, same reduction: 1.107 yr | relative |
| **A2** | cos-weighted RMSE of AOA150 against CLaMS minus CLaMS' own 10S–10N age at 150 hPa (0.092 yr), model put on the CLaMS grid, 100–5 hPa, \|lat\| ≤ 80° | 0 | absolute, < 0.4 yr |
| **A3** | age step from 10S–10N to 50–70° (both hemispheres), nearest level to 55 hPa | CLaMS step: 2.79 yr | relative |
| **A4** | transit 70 → 10 hPa: difference of the 10S–10N age profile, linear in log p | range 2.3–3.2 yr (CLaMS 3.0, WACCM6 2.5) | relative to the midpoint 2.75 yr |
| **B1-p** | w* 15S–15N (cos-weighted) at p = 100, 70, 50, 30, 10 hPa | WACCM6 h6, same estimator: 0.456, 0.228, 0.198, 0.271, 0.463 mm s⁻¹ | relative |
| **B2-p** | tropical upward mass flux: max − min of Ψ* over \|lat\| ≤ 60° at p = 100, 70 hPa | WACCM6: 12.39, 6.53 × 10⁹ kg s⁻¹ | relative |
| **C1-jetP_H** | winter (DJF north, JJA south) maximum of the monthly-mean u poleward of 30°, nearest level to P = 50, 10 hPa | ERA5 monthly zonal means of the same months on the model's grid | relative |
| **C1-…-lat** | grid latitude of that maximum | ERA5, same | absolute, ± 4° on \|lat\| |
| **C2-H** | winter polar-cap (60–90°, cos-weighted) T at the nearest level to 50 hPa, minus ERA5 | 0 | absolute, ± 2 K |
| **D1-H** | standard deviation over winters of the winter-mean u at 10 hPa, 55–65° | ERA5 CDS daily, same days: 6.51 / 2.16 m s⁻¹ | relative |
| **D2-H** | standard deviation of daily u anomalies at 60°, 10 hPa, within winter | 14.99 / 5.25 m s⁻¹ | relative |
| **D3** | NH major SSWs (u 60°N 10 hPa reversal Nov–Mar, ≥ 20 days apart, not a final warming) | 8 in ERA5 | relative |
| **D4-fwH** | mean final-warming date: first easterly crossing of the 5-day-smoothed u 60° 10 hPa after which u stays below 5 m s⁻¹ (Mar–Jul north, Sep–Jan south), springs common to model and ERA5 | ERA5: 44.4 days after 1 Mar, 89.8 after 1 Sep | absolute, ± 10 days |
| **D4-w** | annual-cycle amplitude (max − min of the monthly climatology) of w* 15S–15N at 70 hPa | 0.185 mm s⁻¹ | relative |
| **D5-H** | standard deviation of daily anomalies of v'T' at 100 hPa, 45–75°, within winter | 11.47 / 7.90 K m s⁻¹ | relative |
| **D6-H** | slope of the 50 hPa polar-cap T anomaly on the trailing 40-day mean 100 hPa heat-flux anomaly, Jan–Mar north, Jul–Oct south | 0.972 / 0.766 K per K m s⁻¹ | relative |
| **D7** | slope of the deseasonalised, detrended monthly 10S–10N age at 55 hPa on the trailing 3-month w* anomaly at 70 hPa | CLaMS age with ERA5 w*: −4.52 yr per mm s⁻¹ | relative |

- H = N or S; anomalies in D are taken against each dataset's own 31-day-smoothed day-of-year
  climatology, so model and reference lose their own seasonal cycle.
- D3 uses jcm-strat's detector (`strat_compare.ssw_dates`), ported, not
  `aide_val_common.detect_ssw_cp07`; the two implement Charlton–Polvani differently.

---

## 3. References

| Reference | Location | Years | Metrics |
|---|---|---|---|
| CLaMS v3.1 driven by ERA5, zonal-mean AGE | `/data/CLaMS/CLaMS_v3/clams_v3.1_era5_zm_lat.zip`, read in memory | 2005–2009 (A), 2006–2010 (D7) | A2–A4, D7 |
| WACCM6 REF-D1 entry age (base 0.47°N, 103 hPa) | `/data/CESM2_REFD1_AOA/AoA_waccm6_refd1.04_AOA1mf_1970-2019_ba_0_100.0_ck_0_50.0.nc` | 2005–2009 | A1, figure |
| CESM2.1.5-WACCM6 histSST h6 (`Vzm`, `VTHzm`, `THzm`) | `aide_val_common.load_h6("1996-2014", …)`, interfaces from the model top to 652 hPa, 1e35 masked | 1996–2014 | B1, B2 |
| ERA5 monthly zonal means (`uzm`, `tzm`, 25 levels) | `AIDE-atmosphere/data/era5/era5_zm_monthly_<y>.nc` | 2006–2010 | C1, C2 |
| ERA5 CDS daily zonal means (u, v, T, θ, v'θ', 17 levels 1–250 hPa, 2.5°) | `AIDE-atmosphere/data/era5/era5_<y>.npz` | 2006–2010 | D1–D7 |

- Every path is in one block, `transport_common.REFS`, with the years beside it.
- The ERA5 files are copies made by `scripts/copy_era5_refs.sh` from jcm-strat's cache, only
  the ten a metric reads; `scripts/era5_refs.sha256` is their manifest, and the script checks
  it. They are identical to the files jcm-strat phase 29 read. `AIDE-atmosphere/data/` is
  gitignored.
- **Period matching.** C, D and D7 use the same days as the model. A uses 2005–2009, the
  years jcm-strat's `sweep_score.py` uses; taking 2006–2010 instead changes no A group score by
  more than 0.01 (P20b's A2 err 1.05 → 1.00). B uses the whole 1996–2014 segment: WACCM6
  histSST ends in 2014, so B is not period-matched (protocol §5, rule 3, binds on the tiers,
  not here).

---

## 4. Scoring and tolerances

jcm-strat `relative_scores.py`, unchanged:

- **Relative metrics**: rel = model / ref − 1 (−40 % is 40 % too small, also for a negative
  slope), err = \|rel\| / 0.25. A4 is scored against the midpoint 2.75 yr of its range.
- **Absolute metrics** keep their bound: A2 (< 0.4 yr), C1 latitudes (± 4°, on \|lat\|), C2
  (± 2 K), D4 final-warming dates (± 10 days); err = \|difference\| / bound.
- **Status**: pass err ≤ 1, near err ≤ 1.5, fail above. A metric that cannot be formed fails
  (err 3).
- **Group score** = mean err of the group's metrics. **Combined** = mean of the A, B and C
  group scores; D is reported on its own, so a combined score stays comparable with jcm-strat's
  earlier phases.
- These tolerances are the scorecard's. They are not derived from CESM interannual σ as the
  tier thresholds are, and they do not replace or modify any tier verdict.

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

What the interpolation changes, against the nearest level (from `port_check.md`):

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
jcm-strat memory cap (`memcap.sh`, 96 G); the whole chain takes under two minutes on CPU.

```bash
bash copy_era5_refs.sh                     # ERA5 copies into ../data/era5/, checked against the manifest
$PY 22_transport_fields.py                 # every registered source; or name them
                                           #   -> output/22_transport__references.json
                                           #   -> output/22_transport__<run>.json
$PY 22b_transport_port_check.py            # evidence only -> validation_results/transport/port_check.md
$PY 23_transport_scorecard.py              # -> validation_results/transport/transport_scorecard__<stamp>.md / .json
$PY 24_transport_figures.py                # -> four PNGs with the same stamp
```

- `output/` is gitignored; `validation_results/transport/` is generated and committed. Never
  edit either by hand.
- The stamp is `__<set>__<production date>`, as `17` stamps its report; the default set name is
  `jcm-strat-phase29-2006-2010`.
- Re-running the chain reproduces every committed artefact bit for bit under the pinned
  `requirements.txt` (checked 2026-11-10).

---

## 7. Port check

`22b_transport_port_check.py` compares every metric of every run with jcm-strat's
`relative_metrics.json` (the scored value) and `success_metrics.json` (the unrounded model and
reference values). Agreement: both values within 0.5 % of jcm-strat's, or within the rounding
jcm-strat applied before scoring. B, D4-w and D7 are checked at the nearest level.

| Group | Cells agreeing | Largest model-value difference | Cause of what remains |
|---|---|---|---|
| A | 32 / 32 | < 0.001 % | — |
| B (nearest level) | 56 / 56 | 0.008 % | CESM a and g instead of jcm-strat's 6.371e6 m and 9.80665 m s⁻² |
| C | 66 / 80 | C2: ≤ 0.022 K; jets ≤ 0.37 % | one 6-hourly frame per month (below) |
| D | 104 / 104 | 0.008 % | κ = R/c_p instead of 0.2857 in v'T' |

- **C.** jcm-strat's `teq_score.py` takes calendar-month means of the 6-hourly frames by
  timestamp; the port takes them from the daily means, whose days run 06Z–00Z
  (`daily_zm.py`). Each month loses its opening 00Z frame and gains the next month's. The shift
  estimated from the daily series matches all 16 C2 differences within jcm-strat's 0.01 K
  rounding, and moves L0's 50 hPa NH jet by one grid row (45.70°N instead of 47.56°N; the two
  rows differ by 0.0004 m s⁻¹). All 14 disagreeing cells are these.
- Fed the nearest-level values, `23`'s group scores reproduce jcm-strat's `relative_scores.json`
  for all eight runs, except C (L0 0.98 → 1.02, others ± 0.01) for the reason above.

---

## 8. Result, 2006–2010

Transcribed from
[transport_scorecard__jcm-strat-phase29-2006-2010__2026-11-10.md](../validation_results/transport/transport_scorecard__jcm-strat-phase29-2006-2010__2026-11-10.md);
the full metric table and the four figures are beside it.

| | P20b | AL7 | AL7h | AL7t | L0 | LAL7 | QBO10d | LQBO10d |
|---|---|---|---|---|---|---|---|---|
| A age of air | 0.57 | 1.18 | 0.89 | 1.21 | 0.99 | 1.70 | 0.91 | 1.44 |
| B residual circulation | 1.05 | 1.00 | 0.96 | 1.01 | 0.98 | 1.06 | 0.86 | 0.78 |
| C polar vortex | 1.09 | 0.26 | 0.52 | 0.21 | 1.02 | 0.45 | 0.79 | 0.47 |
| **combined A–C** | **0.90** | **0.81** | **0.79** | **0.81** | **1.00** | **1.07** | **0.85** | **0.90** |
| D variability | 1.35 | 1.28 | 1.11 | 1.29 | 1.17 | 1.45 | 0.93 | 1.39 |
| passed A–C | 13 / 21 | 17 / 21 | 17 / 21 | 15 / 21 | 14 / 21 | 15 / 21 | 15 / 21 | 14 / 21 |
| passed D | 4 / 13 | 6 / 13 | 6 / 13 | 5 / 13 | 6 / 13 | 5 / 13 | 7 / 13 | 4 / 13 |

- Against jcm-strat's own table the changes are in B (interpolation, §5) and in L0's C (§7).
  B1-100 fails in every run (−50 to −72 %), B1-50 fails or is near in every run (+32 to
  +63 %), B2-70 passes in every run.

> **Interpretation** — The interpolated WACCM6 reference at 100 hPa is 14 % higher than the
> nearest-level one, so the lower-branch deficit jcm-strat reported (w* at 100 hPa about half
> of WACCM6's) is larger here, not smaller. The ranking of the runs on the combined score is
> the same as jcm-strat's: AL7h first, then AL7 and AL7t, LAL7 last.

---

## 9. How to score a new model

1. **Write the contract fields** (§1) as daily zonal means on pressure levels, 2006–2010, to
   one or more netCDF files with dims `(time, plev, lat)` and the names `Uzm`, `Vzm`, `THzm`,
   `VTHzm`, `AOA150` (`T`, `Wzm` optional). Form `VTHzm` from instantaneous 3-D fields before
   any time averaging; give `AOA150` in years, reset to zero wherever p > 150 hPa.
2. **Register it** in `transport_common.SOURCES`:
   `"NAME": dict(reader="netcdf", path="/path/to/files_*.nc", long_name="…", note="…",
   independent=True, time_resolution="daily")`.
3. **Check it**: `$PY 22_transport_fields.py --inspect NAME` prints shapes, units, ranges and
   the number of levels in 5–150 hPa. The levels must reach well above 10 hPa (§1).
4. **Score it**: `bash copy_era5_refs.sh`, then `$PY 22_transport_fields.py NAME`,
   `$PY 23_transport_scorecard.py NAME --name SET`, `$PY 24_transport_figures.py`.
5. **Read it with its limits**: §10 and the interim w* note of §5.

- A model that does not cover 2006–2010 cannot be scored as configured: C and D compare the
  same days of ERA5 and CLaMS. Changing `transport_common.WINDOW` needs the ERA5 files of the
  new years (deferred, §12).
- A free-running model is compared with the observed years climatologically; jcm-strat's runs
  are nudged to ERA5 in the troposphere and in the tropical zonal-mean wind, so for them the
  same-day comparison is tighter than for a free run.

---

## 10. Limitations

- **Five years.** D1 is a standard deviation over 4–5 winters and D3 a count of a few events;
  their sampling error is comparable to the 25 % tolerance.
- **Nearest levels** in A (55 hPa), C (50, 10 hPa) and D (10, 50, 55, 100 hPa) make those
  metrics depend on the model's vertical grid; only w* and Ψ* are interpolated.
- **Jet latitudes** are grid latitudes (1.87° on the jcm-strat grid).
- **B is not period-matched** and uses a different estimator from the tiers (§3, §5).
- **References disagree with each other** — A4's range spans CLaMS 3.0 and WACCM6 2.5 yr; A1
  takes WACCM6's entry age while A2 and A3 take CLaMS.

---

## 11. Decisions

Settled by the user before the port:

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

Taken in the port without confirmation — listed in the PR for review:

- Readers and `SOURCES` live in `transport_common.py`, in 20's form, not in 20 itself.
- The `netcdf` reader takes a glob, dims named `time`, `plev`, `lat`, derives `T` from `THzm`
  with κ = R/c_p, masks values beyond 1e20.
- A references stay on jcm-strat's 2005–2009 (effect of 2006–2010 in §3).
- CESM constants (a, g, κ = R/c_p) instead of jcm-strat's; ≤ 0.008 % on any metric (§7).
- WACCM6 h6 read from the model top to 652 hPa with the sentinel masked; jcm-strat read all
  71 interfaces unmasked. No metric level is affected.
- Ψ* is interpolated to the target pressure before its latitude range is taken.
- The D4-w and D7 w* is interpolated too, as one estimator for every w*.
- A, C and D stay at the nearest level, as jcm-strat.
- D3 keeps jcm-strat's SSW detector.
- Monthly means for C from the daily means (§7).
- Port-check criterion and a separate evidence script `22b`.
- Stamped artefact names, the scorecard JSON committed beside the report.
- Heatmap colours an absolute metric as a fraction of its bound, the bound at 25 %.

---

## 12. Deferred

- **ERA5 1990–2022 context for D** (`variability_score.py --long`): how much five-year D
  statistics scatter in ERA5. Not ported; its 33 yearly files were not copied.
- **ERA5 TEM w* context** (jcm-strat `era5_tem_monthly_1996–2001`): shown by jcm-strat beside
  B, never scored. Not ported.
- **CDS download script** (jcm-strat `fetch_era5_daily_zm.py`): the ERA5 files are copies of
  jcm-strat's; regenerating them from the CDS is not part of this repository yet.
- **Another window**: scoring years other than 2006–2010 needs the ERA5 and CLaMS files of
  those years.
