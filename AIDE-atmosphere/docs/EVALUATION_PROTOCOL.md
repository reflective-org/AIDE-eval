# Evaluation protocol

Two tiers:

- **Tier 1** — a 5-year regression check, run on every new model version
- **Tier 2** — a 35-year validation

Every number below is transcribed from `output/16_anchors_45yr.json`, produced by
`scripts/16_anchors_45yr.py`. **Both tiers are anchored on CESM 1970–2014, the whole
record** — 44 annual years, 42 DJF and 45 JJA seasons. Nothing here has been run against an
AIDE-WACCM rollout: these are the criteria, plus the machinery exercised on a window of
CESM's own output (§3.4).

The anchor uses every year the run affords, which makes σ the best estimated and the
variance windows the tightest available. It also leaves no CESM output held back to test the
thresholds with; §2 is that limitation and where the supporting evidence now lives.

Results and interpretation are kept apart. Tables and the bullets under them state what was
measured; every inference drawn from them is in a labelled `Interpretation` block.

---

## The two tiers at a glance

| | **Tier 1 — screening** | **Tier 2 — validation** |
|---|---|---|
| Length | 5 years | 35 years |
| Anchor | 1970–2014 (45 yr) | 1970–2014 (45 yr) |
| Runs on | every new model version | a configuration that passed tier 1 |
| Scores | **individual years** | the **rollout mean** and the **variance** |
| Question | is the configuration broken? | is the configuration correct? |
| Thresholds | loose, empirical, **provisional** | derived, not adjustable by inspection |
| Expected to change | yes — tighten as failure modes emerge | no |
| A pass means | continue | the configuration meets the validation targets |

The two tiers answer different questions. Tier 1 is a regression test: low cost, a low
false-alarm rate on a healthy model, and sensitive to a broken configuration. Tier 2 is the
scientific claim.

---

## 1. Tier 1 — 5-year screening

### What is checked

**Every individual year, against a ±3σ band around the CESM 1970–2014 mean.** A 5-year
rollout gives 5 values per diagnostic and all 5 are checked. This is not a check on the
5-year mean: a regression test has to see a single catastrophic year, which averaging hides.

σ is detrended within the anchor wherever the trend is significant at p < 0.05, marked `*`.
Three diagnostics trigger it.

| Diagnostic | Units | Anchor mean | σ | **Accept (any single year)** | Anchor n | Widest anchor year |
|---|---|---|---|---|---|---|
| Mass flux, 70 hPa | 10⁹ kg s⁻¹ | 9.4487 | 0.2128* | **8.810 – 10.087** | 44 | +3.45σ (1978) |
| w̄*, 10°S–10°N | mm s⁻¹ | 0.1917 | 0.0144* | **0.1485 – 0.2350** | 44 | +2.57σ (1992) |
| Vortex NH, DJF | m s⁻¹ | 24.8 | 5.6 | **8.0 – 41.6** | 42 | +2.17σ (2014) |
| Vortex SH, JJA | m s⁻¹ | −81.6 | 2.5 | **−89.0 – −74.2** | 45 | +3.05σ (2011) |
| Polar cap T, NH | K | 209.55 | 1.90* | **203.8 – 215.3** | 42 | +2.35σ (2011) |
| Polar cap T, SH | K | 192.13 | 2.02 | **186.1 – 198.2** | 45 | +2.75σ (1992) |

Plus one count, which needs no band: **0–7 major NH SSWs over the 5 winters**, against 3.4
expected at the anchor rate of 0.689/winter (31 events in 45 winters, 95% Poisson interval).

### The anchor against its own band

Scoring all 45 anchor years back against the band they define gives **3 exceedances in 262
year-checks**: mass flux in 1971 and 1978, SH vortex in 2011.

> **Interpretation** — A band built from a detrended σ but centred on a mid-trend mean will
> flag the ends of a trending record. Both mass-flux exceedances are in the 1970s, before
> most of the BDC acceleration (D9). The 3 are not evidence that ±3σ is too tight; they are
> the forced trend appearing in a statistic that has had the trend removed from its width but
> not from its centre. A model flagged in the same way needs its period checked before
> its physics.

### Why ±3σ

The band clears the whole anchor record apart from those 3 checks, and multiplicity is the
reason not to start tighter. Six diagnostics × 5 years is **30 simultaneous checks per
screening run**. For a model that is perfect and merely Gaussian:

| Band | False alarm per check | Chance of ≥1 flag per screening run |
|---|---|---|
| ±2.0σ | 4.55% | **75.3%** |
| ±2.5σ | 1.24% | 31.3% |
| ±3.0σ | 0.27% | **7.8%** |
| ±3.5σ | 0.05% | 1.4% |

> **Interpretation** — ±2σ is unusable as a hard gate: it would flag three runs in four on a
> flawless model. ±3σ is the loosest band that is still informative and the tightest that
> stays quiet. If tier 1 tightens, this is the number to move, and the sensible path is to
> keep ±3σ as the hard gate and add tighter bands as advisory flags that are logged rather
> than failed on.

> **Settled for now: ±3σ.** Agreed 2026-08-18, to be revisited once real failure modes have
> been seen.

### Advisory, not gating, at tier 1

- **The 5-year mean**, held to `max(0.5σ, 1.96σ/√5) = 0.877σ`. Logged because it catches a
  uniform drift that no individual year is extreme enough to trip. Not a gate: at n = 5 the
  detection branch dominates so heavily that the threshold is nearly 1σ, which would pass
  almost anything (appendix A).

### What tier 1 does not do

It does not test variance, the mechanism relations, the forced trend, or anything about the
distribution.

> **Interpretation** — A tier-1 pass is not evidence that a configuration is good, only that
> it is not obviously broken. Do not report it as a validation.

---

## 2. What this anchor cannot verify

The standing rule is that a threshold CESM cannot meet on a different sample of its own
output is not a threshold (D12). The operational anchor cannot satisfy that rule: it spans
1970–2014, so **no CESM year lies outside it** and no held-out sample exists.

The evidence therefore comes from the anchors this one replaced, archived in
[`../stale/`](../stale/README.md):

| Test | Anchor | Scored on | Result |
|---|---|---|---|
| Tier 1, ±3σ, per year | 1996–2014 | 1970–1995, as five 5-year screening runs | 28 of 30 block verdicts pass; both failures mass flux, both pre-1980 |
| Tier 2, mean and variance | 1980–2014 | 1970–1994 | 4 of 6 mean, 6 of 6 variance |

> **Interpretation** — What transfers from those runs is the *method*, not the numbers: a
> ±3σ per-year gate and a 0.5σ mean tolerance are the right size on CESM output that did not
> set them, and the diagnostics that failed did so for an identifiable physical reason rather
> than from noise. The thresholds in this document are the same construction on a longer
> sample. They inherit that support and cannot add to it — which is the price of using the
> whole record, and the reason the archived material is kept rather than deleted.

---

## 3. Tier 2 — 35-year validation

**Anchored on 1970–2014**, the same anchor as tier 1. 35 years is past the 15.4-year
crossover, so every mean tolerance is the full-strength **0.5σ**; running longer buys no
further tightening on the mean.

The anchor joins two separate CESM runs across the 1995/96 restart. σ is detrended where the
trend is significant at p < 0.05, marked `*`.

### 3.1 Mean

| Diagnostic | Units | Anchor mean | σ | Tolerance | **Accept (rollout mean)** | as % of mean |
|---|---|---|---|---|---|---|
| Mass flux, 70 hPa | 10⁹ kg s⁻¹ | 9.4487 | 0.2128* | ±0.1064 | **9.3423 – 9.5552** | 1.1% |
| w̄*, 10°S–10°N | mm s⁻¹ | 0.1917 | 0.0144* | ±0.0072 | **0.1845 – 0.1990** | 3.8% |
| Vortex NH, DJF | m s⁻¹ | 24.8 | 5.6 | ±2.8 | **22.0 – 27.6** | 11.3% |
| Vortex SH, JJA | m s⁻¹ | −81.6 | 2.5 | ±1.2 | **−82.8 – −80.3** | 1.5% |
| Polar cap T, NH | K | 209.55 | 1.90* | ±0.95 | **208.60 – 210.50** | 0.5% |
| Polar cap T, SH | K | 192.13 | 2.02 | ±1.01 | **191.12 – 193.14** | 0.5% |

The DJF rows rest on 42 seasons and the JJA rows on 45; the annual rows on 44 years. One
DJF season is lost at each record edge and one annual year at the restart.

### 3.2 Variance

| Metric | 95% acceptance window | Width |
|---|---|---|
| **Interannual σ ratio** (rollout σ ÷ anchor σ) | **0.68 – 1.32** | 0.64 |
| **Daily DJF σ ratio**, u at 60°N | **0.88 – 1.12** | 0.24 |

The daily series gives 6.8 effective samples per winter — 4013 DJF days, 14-day
decorrelation time — against one for the annual series.

> **Interpretation** — Both windows are tighter than the same construction on the 35-year
> anchor gave (0.66–1.34 and 0.87–1.13), because the reference σ now rests on the whole
> record. That is the main practical gain from the longer anchor. The daily ratio remains the
> sharper test and is what will catch an emulator that has smoothed away its own weather.
> Detrend both sides before forming the interannual ratio: over the pooled record the
> mass-flux σ is 54% larger raw than detrended, so an undetrended ratio measures trend, not
> variability.

### 3.3 Counts and relations

| Check | Anchor | **Accept at 35 years** |
|---|---|---|
| Major NH SSW, count | 0.689/winter (31 in 45) | **15 – 34**, expect 24.1 |
| R1 wave → vortex, slope | −0.936 | resolvable to ±0.36 |
| R2 thermal wind, slope | −2.062 | resolvable to ±0.57 |

### 3.4 Running a validation

```bash
$PY scripts/17_validate.py FIRST_YEAR LAST_YEAR     # default 1996 2014
$PY scripts/18_validation_figures.py
```

Writes [`../validation_results/validation_result.md`](../validation_results/)
and four figures beside it, one per test family. Climate-model data enters through
`17_validate.climate_model_series`, which is the one function to replace when the model is a
model rollout rather than a window of CESM's own record.

The committed result scores **CESM 1996–2014**: tier 1 5/6, tier 2 mean 5/6, variance 7/7,
SSW pass, mechanism 1/2.

> **Interpretation** — That model lies inside the anchor, so it is a self-consistency
> check on the machinery and not a validation. Its two mean-side failures are both predicted:
> mass flux at +1.04σ is D9 with the sign reversed, 1996–2014 being the post-acceleration
> half of the record, and R1 is the relation D10 demoted. Any pass there is weaker evidence
> than the corresponding failure, and the report labels itself accordingly.

---

## 4. Input requirements

What a model has to supply to be scored against the thresholds above. Everything in this
section is what `scripts/aide_val_common.py` reads and what the diagnostics consume. The
CESM values quoted throughout are the reference implementation, not a constraint on the
climate model's own grid (§5, rule 2).

### 4.1 Minimum request

**Tier 1 needs 5 years, tier 2 needs 35 years.** Both need the same fields, at the same
resolution, on the same grid. Two forms of the request, depending on what the archive holds:

| | **A zonal-mean TEM tape** | **A standard pressure-level archive** |
|---|---|---|
| Fields | ū, v̄, w̄ (log-pressure), θ̄, v'θ' | `ua`, `va`, `wap`, `ta` (CF/CMIP names) |
| Shape | (time, level, lat), already zonally averaged | (time, level, lat, lon) daily 3-D — or zonal means **plus** 3-D `va` and `ta` |
| Temporal | daily means | daily means |
| Levels | pressure levels bracketing 10, 50, 70 and 100 hPa, with one level beyond each end | same |
| Latitude | global, −90° to +90° | same |
| Derivation | none | w̄ from `wap`, θ̄ from `ta`, v'θ' from 3-D `va` and `ta` (§4.6) |

- **v'θ' is the one field that cannot be recovered from zonal means.** An archive holding
  only zonal means must also hold the eddy heat flux itself, or the 3-D daily `va` and `ta`
  from which to form it.
- **θ̄ needs no 3-D data.** On a pressure surface the transform is linear in T, so
  `θ̄ = T̄ (1000/p)^κ` exactly.
- **`wap` is ω, not w.** It has to be converted (§4.6); using it, or a geometric vertical
  velocity, in place of log-pressure w̄ is an 11% error on all upwelling.
- CESM's `cam.h6` tape is the first column. §4.2–§4.6 give the exact definitions the two
  columns have to meet.

### 4.2 The five fields

| Symbol | CESM `cam.h6` name | Units | Definition | Diagnostics that need it |
|---|---|---|---|---|
| ū | `Uzm` | m s⁻¹ | zonal-mean zonal wind | vortex NH, vortex SH, SSW count, daily DJF σ |
| v̄ | `Vzm` | m s⁻¹ | zonal-mean meridional wind | v̄* (TEM); not scored on its own |
| w̄ | `Wzm` | m s⁻¹ | zonal-mean **log-pressure** vertical velocity, `w = −Hω/p` | w̄*, mass flux |
| θ̄ | `THzm` | K | zonal-mean potential temperature, reference 1000 hPa | polar cap T; `dθ/dz` in the TEM term |
| v'θ' | `VTHzm` | K m s⁻¹ | zonal-mean **eddy** meridional heat flux | TEM streamfunction ψ; R1 heat flux at 100 hPa |

- Nothing else is read. `UVzm` and `UWzm` are on the CESM tape and are not used by any tier-1
  or tier-2 diagnostic.
- **All five must sit on one common set of levels and one common latitude grid.**
  `tem_residual` differences them against each other level by level.
- `VTHzm` is the eddy flux `v'θ'`, not the product of the zonal means. Passing the full
  product where the eddy flux is expected reverses the sign of tropical w̄* (§4.6).

### 4.3 Temporal resolution

**Daily means are required. A monthly archive is not sufficient.**

- Two scored quantities cannot be formed from monthly output: the **SSW count** (needs the
  day the wind reverses) and the **daily DJF σ ratio**, which §3.2 identifies as the sharper
  of the two variance tests.
- Every other diagnostic is reduced from the same daily series. A monthly archive would
  support the mean tests of §3.1 and nothing else.
- CESM's calendar is `noleap` (365 days). A different calendar changes only the day counts
  in the completeness rules below.
- Each daily value is labelled by the calendar day of its own mean.

| Reduction | Rule | Completeness requirement |
|---|---|---|
| Annual mean | calendar-year mean of the daily series | ≥ 350 days, else the year is dropped |
| DJF season | Dec–Jan–Feb, labelled by the January year (Dec 1996 → season 1997) | ≥ 81 days (`0.9 × 30 × 3`) |
| JJA season | Jun–Jul–Aug, labelled by its own year | ≥ 81 days |
| SSW season | NH, Nov 1 – Mar 31, labelled by the January year | ≥ 100 days |
| Daily DJF σ | all DJF days pooled, no seasonal reduction first | decorrelation time τ = 14 days assumed |

- Record length: **5 years for tier 1, 35 years for tier 2** (§1, §3).
- Incomplete DJF seasons at both record edges are dropped, leaving one fewer DJF season
  than annual years. JJA lies inside a calendar year and loses none.

### 4.4 Vertical grid

| Diagnostic | Level or layer |
|---|---|
| Mass flux | 70 hPa |
| w̄*, 10°S–10°N | 70 hPa |
| Vortex NH, vortex SH, SSW | 10 hPa |
| Polar cap T, both hemispheres | layer 10–50 hPa |
| R1 heat flux (tier-2 mechanism) | 100 hPa |

- **Pressure levels, monotonically increasing**, with the pressure of each level known
  independently of surface pressure. On CESM's hybrid grid every level above 182 hPa has
  `hybi = 0`, so interface pressure is `hyai·P0` and no `PS` is needed.
- Interpolation onto 10, 70 and 100 hPa is **linear in ln p**. The level set must bracket
  each target, with at least one further level beyond each end of the 10–100 hPa range so
  the centred vertical derivative is defined at the brackets.
- The TEM term needs `dθ/dz` on log-pressure height `z = −H ln(p/1000)`, by centred
  differences over the supplied levels.
- Polar cap T is a mean over the levels falling inside 10–50 hPa, weighted by Δ ln p.
- CESM reference: 71 hybrid interface levels, of which the 41 spanning 0.099–652 hPa are
  loaded; 26 lie between 1 and 200 hPa, and 8 inside the 10–50 hPa polar cap layer. The
  brackets are 9.55 / 11.86 hPa for 10, 66.80 / 80.70 for 70, and 94.94 / 111.69 for 100.
- Fill values must be masked before any reduction (§5, rule 4). On CESM this is the `1e35`
  below-surface sentinel; on another model it is whatever that model writes.

### 4.5 Horizontal grid

- **Latitude must be global, −90° to +90°.** The mass-flux integral runs over every latitude
  where w̄* > 0 at 70 hPa, not over a prescribed tropical band.
- Latitude bands used: **10°S–10°N** (w̄*), **60–90°N** and **60–90°S** (polar cap T),
  **45–75°N** (R1 heat flux). Band means are **cos(lat)-weighted with inclusive endpoints**,
  not Gaussian-weighted.
- Point latitudes: **exactly ±60°**, reached by linear interpolation in latitude. CESM's
  0.9424° grid has no point at 60°.
- Zonal means only. No longitudinal information enters any diagnostic.
- CESM reference: 192 latitudes, 0.9424° spacing (f09).

> **Interpretation** — The grid above is what CESM supplies, not what a climate model must
> match. Rule 2 in §5 scores on the emulator's own grid with CESM reduced onto it, so a
> coarser model grid is admissible; an incomplete latitude range is not, because it changes
> the mass-flux integral itself.

### 4.6 Deriving the five fields from a standard archive

Most archives carry ω, T and 3-D winds rather than a zonal-mean TEM tape.

| Available (CF/CMIP name) | Derive | Constraint |
|---|---|---|
| `ua`, `va` | ū, v̄ — zonal mean, no transform | Must be on the same pressure levels as the rest |
| `wap` = ω (Pa s⁻¹) | `w = −Hω/p`, H = 7000 m | Geometric vertical velocity is not a substitute: using it is an 11% error on all upwelling |
| `ta` = T (K) | `θ = T (1000/p)^κ`, κ = 0.2857 | p in hPa. Linear in T on a pressure surface, so `θ̄ = T̄ (1000/p)^κ` |
| 3-D `va` and `ta` | `v'θ' = [vθ] − [v][θ]`, formed on the native 3-D grid **before** zonal averaging | Cannot be recovered from zonal means alone |

- Polar cap T is reconstructed from θ̄ by the reference path. An archive that carries `ta`
  can use it directly; the two are the same quantity.
- A coarse pressure archive can bracket 10, 50, 70 and 100 hPa and still degrade two
  diagnostics: `dθ/dz` in the TEM term is differenced over the levels supplied, and the
  polar cap average is taken over the levels inside 10–50 hPa, of which CESM has 8 (§4.4).

Constants, at CESM's own values:

| Constant | Value | Used in |
|---|---|---|
| a (Earth radius) | 6.37122 × 10⁶ m | mass flux, w̄* |
| g₀ | 9.80616 m s⁻² | mass-flux density |
| R | 287.058 J kg⁻¹ K⁻¹ | κ |
| c_p | 1004.64 J kg⁻¹ K⁻¹ | κ |
| κ = R/c_p | 0.2857 | θ ↔ T; computed from R and c_p, not hard-coded |
| H (scale height) | 7000 m | log-pressure w and z |
| p_ref | 1000 hPa | θ, log-pressure height |
| ρ at the flux surface | `p / (g₀H)` | mass flux (Butchart et al. 2010 / CCMVal convention) |

> **Interpretation** — These are CESM's constants, and the thresholds were derived with them.
> A rollout scored with different constants is being compared against a different quantity.
> This is D8 applied to the inputs rather than to the estimator: the two standard routes to
> w̄* differ by 10.8%, which is larger than the tier-2 mass-flux tolerance.

---

## 5. Rules that apply at both tiers

These are not negotiable per tier. They are what makes a comparison meaningful.

1. **Pin the estimator.** Model and truth go through `aide_val_common.tem_residual`. The
   two standard routes to w̄* differ by 10.8% on identical data — larger than the tier-2
   mass-flux tolerance itself (D8).
2. **Score on the emulator's grid**, with CESM reduced onto it, not the reverse.
3. **Detrend both sides, and period-match anything that overlaps CESM in time.** The
   anchor spans a real forced trend, so its mean sits mid-trend while its σ has the trend
   removed (D9). A free rollout with no GHG forcing cannot reproduce that trend, so an
   absolute upwelling comparison is the wrong test for one; score it against the CESM years
   it covers, or detrend both sides. This binds hardest on mass flux and w̄*, the two
   diagnostics whose σ is detrended in §1 and §3.1.
4. **Mask the `1e35` sentinel** with `abs(x) < 1e20` before any reduction. There is no
   `_FillValue`; an unguarded mean returns ~1e33.
5. **Report the rollout length with every number.** A threshold without its n is not a
   threshold.

---

## 6. Provenance

```bash
$PY scripts/02_reference_stats.py         # reads the h6 tape
                                          # writes output/02_reference_stats.json
$PY scripts/02b_trends.py                 # reads 02
                                          # writes output/02b_trends.json
$PY scripts/07_period_split.py            # reads the h6 tape
                                          # writes output/07_period_split.json
$PY scripts/16_anchors_45yr.py            # reads 07
                                          # writes output/16_anchors_45yr.json
                                          #   -> sections 1 and 3 of this document
$PY scripts/17_validate.py \
      --climate-model NAME 1996 2014      # reads 16 and 07; NAME defaults from the window
                                          # writes output/17_validation__<stamp>.json
                                          #   and validation_results/
                                          #     validation_result__<stamp>.md
$PY scripts/18_validation_figures.py      # reads 16, 07 and the newest 17 output
      [--climate-model NAME]              # writes validation_results/*__<stamp>.png
```

Those six scripts in that order, plus `aide_val_common.py` and `report_layout.py`. `16`,
`17` and `18` read only existing JSON and take seconds; if `output/` is empty, `02`, `02b`
and `07` rebuild it from the tape in about five minutes. See the pipeline order in
[../CLAUDE.md](../../CLAUDE.md).

Sections 1 and 3 change only when `16` changes; §3.4 only when `17` or `18` changes. Do not
edit the tables by hand.

Every scored result is stamped `__<climate model>__<production date>`, so the report, its
four figures and the JSON behind them carry the name of what was scored and the day it was
produced. `18` takes the stamp from `17`'s JSON rather than recomputing it, so a run that
crosses midnight cannot split its own figures from its report.

The split-anchor scripts that used to sit here, `14` and `15`, are in
[`../stale/`](../stale/README.md). They still run and still reproduce their JSON; §2 says
what they are still cited for.

---

## 7. Open for sign-off

Two choices here were made to produce a usable document and are not derived from the data.
Each is cheap to change.

**Settled, 2026-08-18:**

- **±3σ as the tier-1 gate**, to be revisited once real failure modes have been seen.

**Settled, 2026-08-19:**

- **Both tiers anchored on the whole 1970–2014 record**, so the thresholds are fixed and
  model-agnostic rather than matched to a CESM sub-period. §2 states the consequence: the
  anchor cannot be verified out of sample, and the evidence that the construction holds up
  is the archived split-anchor material.

**Still open:**

1. **Tier 1 gates on individual years only**, with the 5-year mean advisory. The brief
   specified individual years; the mean check was added here and can be dropped or promoted
   to a gate.
2. **SSW count and mechanism slopes are reported at tier 2** although the brief
   specified mean and variance. They are part of the existing target set and become
   resolvable at this length, so they are included and flagged rather than omitted.

---

## Appendix A · Where `0.5σ` and `1.96σ/√n` come from

No number in this appendix is transcribed from `output/`. It is the derivation of the rule,
not a result; the σ values it is applied to are in sections 1 and 3.

**σ is CESM's own interannual scatter** — the year-to-year variation in its DJF vortex or its
70 hPa mass flux under unchanged forcing. It is not an error bar on the model and not an
observational uncertainty. It is the unit every bias tolerance is measured in.

`max(0.5σ, 1.96σ/√n)` combines two statements that answer different questions:

| Branch | Question it answers | Depends on n? |
|---|---|---|
| `0.5σ` | What counts as an *acceptable* bias? | No |
| `1.96σ/√n` | What bias can an n-year run *resolve*? | Yes |

- **`0.5σ` is the scientific tolerance** — the emulator's mean must sit inside half of
  natural variability. Rollout length does not change what CESM's variability is, so this
  branch does not move with n.
- **`1.96σ/√n` is the detection floor.** An n-year mean carries sampling error σ/√n, so a
  threshold tighter than 1.96σ/√n lies inside the 95% noise band of the estimate itself: a
  model sitting on it cannot be distinguished from an unbiased one.

Taking the larger of the two demands no more accuracy than the physics requires and no more
than the run length can resolve. The branches cross where `0.5 = 1.96/√n`, i.e. **n =
(1.96/0.5)² = 15.4 years**: below that, detection binds; above it, tolerance binds. This is
why tier 1's 5-year mean threshold is 0.877σ and therefore advisory (§1), and why tier 2 at
35 years gets the full-strength 0.5σ on every diagnostic (§3.1).

## Appendix B · Decisions cited above

Condensed from the decision record kept during derivation; the numbering is retained so the
citations in the text resolve.

**D5 · The threshold rule is `target(n) = max(0.5σ, 1.96σ/√n)`.** Every bias tolerance is
the larger of half of natural variability and what an n-year mean can resolve at 95%
confidence. Derivation in appendix A.

**D7 · The interannual σ ratio is not usable at short n.** At n = 10 its 95% window is
0.43–1.57 — nothing realistic fails it, so a pass is uninformative. The daily DJF series
gives about 6.4 independent samples per winter (14-day decorrelation time) instead of one,
which is a usable test. §3.2 is where the interannual ratio becomes usable again at n = 35.

**D8 · A target is a statement about a *specific estimator* on a *specific grid*.** The two
standard routes to w̄* differ by 10.8% on identical data — larger than the tier-2 mass-flux
tolerance. Model and truth must go through `aide_val_common.tem_residual`.

**D9 · The upwelling targets are period-matched, not absolute.** Setting targets on
1996–2014 and scoring 1970–1995 fails both upwelling diagnostics: the earlier mass flux is
4.0% lower, and the trend is +0.31 ×10⁹ kg s⁻¹/decade in 1970–1995 (p < 10⁻⁴) but
statistically absent in 1996–2014 (p = 0.34). This is a protocol change, not a looser
number: AIDE-WACCM's input-only forcings are TOA solar, time-of-day, year-progress and
statics, with no GHGs, so a free rollout has no mechanism to reproduce a GHG-forced trend.
An absolute target would fail an unbiased emulator.

**D12 · A threshold CESM cannot meet on a different sample of its own output is not a
threshold.** The rule stands and is what the archived split anchors were built to satisfy —
the data that sets a threshold cannot also test it. The operational anchor deliberately
forgoes it in exchange for the whole record; §2 is that trade stated explicitly, and any new
diagnostic added here still owes the same held-out test before it is trusted.

**D13 · Reference σ comes from CESM, not from observations.** The right reference for
scoring an emulator *against this CESM run*, and the wrong one for any claim about realism —
see appendix C.

## Appendix C · Limitations that apply to every number here

- **Fixed SST.** No ENSO, so the interannual σ anchoring every tolerance is smaller than the
  real atmosphere's, and the tolerances derived from it are tighter.
- **Volcanic years.** 1970–1995 contains El Chichón and Pinatubo, part of why its upwelling
  variability is larger. The anchor spans them, so their variance is inside every σ here.
- **No held-out CESM sample.** The anchor spans the whole record, so the thresholds cannot be
  tested against CESM output that did not set them (§2, D12). Three of the anchor's own 262
  year-checks fall outside its ±3σ band.
- **The two segments are separate CESM runs** joined at a 1995/96 restart, not one
  continuous integration. Some of the 4% mass-flux step could be a restart discontinuity
  rather than a trend; the significant within-segment trend in 1970–1995 argues against
  that, but the two have not been separated. The tier-2 anchor spans this restart (§3).
- **SSW detection** is a local implementation of Charlton–Polvani, not a community-shared
  catalogue. The rate is consistent with the literature, but central dates have not been
  cross-checked independently.
- **Nothing here has been run against an actual emulator rollout.** Every number is
  CESM-vs-CESM.

---

## See also

- §4 — the fields, resolution and grid a model must supply to be scored
- [../validation_results/validation_result.md](../validation_results/) —
  a scored climate model, with the four figures beside it
- [../stale/README.md](../stale/README.md) — the split-anchor material §2 cites
- Appendix A — where `0.5σ` and `1.96σ/√n` come from
- Appendix B — the decisions cited above by number
- Appendix C — the limitations that apply to every number here
- [../stale/AIDE_WACCM_screening_1970-1995.pdf](../stale/AIDE_WACCM_screening_1970-1995.pdf)
  — the figure behind the §2 table, on the anchors it was built for
- `scripts/16_anchors_45yr.py` — every threshold on this page
- `scripts/17_validate.py`, `scripts/18_validation_figures.py` — the scored result and its
  figures
