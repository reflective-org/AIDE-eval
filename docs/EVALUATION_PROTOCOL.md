# Evaluation protocol

Two tiers:

- **Tier 1** — a 5-year regression check, run on every new model version.
- **Tier 2** — a 35-year validation, run once a configuration is a candidate.

Every number below is transcribed from `output/14_evaluation_tiers.json` and
`output/15_screen_out_of_sample.json`, produced by `scripts/14_evaluation_tiers.py` and
`scripts/15_screen_out_of_sample.py`. **Tier 1 is anchored on 1996–2014, tier 2 on
1980–2014**; the reasons are in each section. Nothing here has been run against an
AIDE-WACCM rollout: these are the criteria, plus both tiers applied to CESM's own earlier
output.

Results and interpretation are kept apart. Tables and the bullets under them state what was
measured; every inference drawn from them is in a labelled `Interpretation` block.

---

## The two tiers at a glance

| | **Tier 1 — screening** | **Tier 2 — validation** |
|---|---|---|
| Length | 5 years | 35 years |
| Anchor | 1996–2014 (19 yr) | 1980–2014 (35 yr) |
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

**Every individual year, against a ±3σ band around the CESM 1996–2014 mean.** A 5-year
rollout produces 5 values per diagnostic; all 5 are checked. The 5-year mean is not the
gate: a regression test has to detect a single extreme year, which averaging removes.

| Diagnostic | Units | Anchor mean | σ | **Accept (any single year)** | CESM anchor years span |
|---|---|---|---|---|---|
| Mass flux, 70 hPa | 10⁹ kg s⁻¹ | 9.670 | 0.220 | **9.010 – 10.330** | −1.82 / +1.62 σ |
| w̄*, 10°S–10°N | mm s⁻¹ | 0.1956 | 0.0127 | **0.1574 – 0.2337** | −1.79 / +1.45 σ |
| Vortex NH, DJF | m s⁻¹ | 25.0 | 5.7 | **7.9 – 42.1** | −2.07 / +2.09 σ |
| Vortex SH, JJA | m s⁻¹ | −81.2 | 2.7 | **−89.4 – −72.9** | −1.38 / +2.58 σ |
| Polar cap T, NH | K | 209.0 | 2.2 | **202.3 – 215.7** | −1.37 / +2.23 σ |
| Polar cap T, SH | K | 192.2 | 2.3 | **185.4 – 199.1** | −1.63 / +1.75 σ |

Plus one count, which needs no band: **1–8 major NH SSWs over the 5 winters** (expected 3.9
at the CESM rate of 0.79/winter, 95% Poisson interval). Only a rollout with no major
warmings, or with more than eight, falls outside the interval.

### Why ±3σ, and the limit on tightening it

- No CESM anchor year of any diagnostic reaches 2.6σ, so ±3σ covers the observed record
  with margin.
- The constraint on tightening is multiplicity: six diagnostics × 5 years is **30
  simultaneous checks per screening run**.

For a model that is unbiased and Gaussian:

| Band | False alarm per check | Chance of ≥1 flag per screening run |
|---|---|---|
| ±2.0σ | 4.55% | **75.3%** |
| ±2.5σ | 1.24% | 31.3% |
| ±3.0σ | 0.27% | **7.8%** |
| ±3.5σ | 0.05% | 1.4% |

> **Interpretation** — ±2σ is not usable as a hard gate: it flags three runs in four from an
> unbiased model, and a gate with that false-alarm rate is ignored. ±3σ is the loosest band
> that still resolves an excursion and the tightest that holds the per-run false-alarm rate
> below 10%. If tier 1 is tightened, ±3σ is the parameter to move; the path with the fewest
> side effects is to keep ±3σ as the hard gate and add tighter bands as advisory flags that
> are logged rather than failed on.

> **Settled for now: ±3σ.** Agreed 2026-08-18, to be revisited once real failure modes
> have been seen. Section 2 is the out-of-sample check on that choice.

### Advisory, not gating, at tier 1

- **The 5-year mean**, held to `max(0.5σ, 1.96σ/√5) = 0.877σ`. It detects a uniform drift
  that no individual year is extreme enough to trip. It is not a gate: at n = 5 the
  detection branch dominates and the threshold is close to 1σ, which admits most
  configurations (appendix A).

### What tier 1 does not do

It does not test variance, the mechanism relations, the forced trend, or the distribution. A
tier-1 pass is not evidence that a configuration is correct, only that it is not broken in
the six diagnostics checked. Do not report it as a validation.

### Period sensitivity of the upwelling diagnostics

- The bands above are anchored on 1996–2014.
- Over the **pooled** 1970–2014 record the mass flux reaches −4.34σ and w̄* −2.98σ. That
  excursion is the forced BDC trend, not natural scatter (D9), and the trend is concentrated
  in 1970–1995.
- The four vortex and temperature diagnostics are insensitive to the choice of record: NH
  vortex and NH polar cap T agree to within 0.01σ, and SH vortex and SH polar cap T widen by
  0.30σ and 0.63σ, on one side only.

> **Interpretation** — A screening band built from the pooled record would accept a rollout
> with a 1970s circulation. Score a rollout against the CESM years it covers, or detrend
> both sides.

---

## 2. Does the tier-1 gate survive out of sample?

The standing rule: a threshold CESM cannot meet on a different sample of its own output is
not a threshold (D12). The ±3σ gate was run against **CESM 1970–1995, treated as five
consecutive 5-year screening runs** — 150 individual year-checks across the six diagnostics,
with the five blocks spanning 1970–1994.

**Figure: [AIDE_WACCM_screening_1970-1995.pdf](AIDE_WACCM_screening_1970-1995.pdf)**
(3 pp — the six diagnostics year by year, the same result as block verdicts, then the
tier-2 test of section 3.3).

| Diagnostic | Offset of 1970–1995 from anchor | Years outside ±3σ | Worst year |
|---|---|---|---|
| Mass flux, 70 hPa | **−1.77σ** | **4 of 25** (1970, 1971, 1975, 1978) | −4.34σ (1978) |
| w̄*, 10°S–10°N | −0.53σ | 0 of 25 | −2.98σ (1974) |
| Vortex NH, DJF | −0.07σ | 0 of 24 | −2.08σ (1977) |
| Vortex SH, JJA | −0.26σ | 0 of 26 | −1.68σ (1983) |
| Polar cap T, NH | +0.42σ | 0 of 24 | +1.74σ (1985) |
| Polar cap T, SH | −0.07σ | 0 of 26 | +2.38σ (1992) |

Results:

- Four of six diagnostics produce no flag in 26 years of out-of-sample CESM.
- The vortex and temperature diagnostics fall inside the band for their whole 1970–1995
  record.
- Mass flux flags 4 of 25 years, all in the 1970s, the worst at −4.34σ in 1978; the period
  as a whole sits 1.77σ below the anchor.
- As screening-run verdicts rather than year-checks: **2 of 30 block verdicts fail**, both
  mass flux, both before 1980. The other 28 pass.

> **Interpretation** — The gate does not fire on a healthy model, which is the required
> behaviour for screening. The mass-flux failures are the outcome D9 predicts: the forced BDC
> trend, not a broken configuration, and the reason the protocol requires scoring a rollout
> against the CESM years it covers. Every flag the ±3σ gate raised on 26 years of CESM output
> pointed at forced physics, none at noise, so a tier-1 flag means "inspect this", not "this
> configuration is broken".

---

## 3. Tier 2 — 35-year validation

**Anchored on 1980–2014** — 35 years, the same length as the rollout it scores, and a
different anchor from tier 1, which stays on 1996–2014.

Two properties of that window, neither of them ideal:

- **It spans the 1995/96 restart.** 1970–1995 and 1996–2014 are two separate CESM runs
  joined at a restart, not one integration. The tier-1 anchor sits inside one run; this one
  does not.
- **It overlaps the test period below by 15 years** (1980–1994). Section 3.3 is therefore a
  consistency check, not the independent out-of-sample test that D12 asks for.

σ is detrended within the window where the trend is significant at p < 0.05, the same rule
`07_period_split.stats_of` applies. Only the mass flux triggers it here (marked `*`).

### 3.1 Mean

35 years is past the 15.4-year crossover, so the detection branch is not binding and every
mean tolerance is the full-strength **0.5σ** — the same target a 20-year or a 45-year
rollout would face. Running longer than 35 years gives no further tightening of the mean.

| Diagnostic | Units | Anchor mean | σ | Tolerance | **Accept (rollout mean)** | as % of mean |
|---|---|---|---|---|---|---|
| Mass flux, 70 hPa | 10⁹ kg s⁻¹ | 9.564 | 0.204* | ±0.102 | **9.462 – 9.666** | 1.1% |
| w̄*, 10°S–10°N | mm s⁻¹ | 0.1965 | 0.0137 | ±0.0068 | **0.1896 – 0.2033** | 3.5% |
| Vortex NH, DJF | m s⁻¹ | 25.7 | 5.4 | ±2.7 | **23.0 – 28.4** | 10.5% |
| Vortex SH, JJA | m s⁻¹ | −81.5 | 2.7 | ±1.3 | **−82.8 – −80.1** | 1.7% |
| Polar cap T, NH | K | 209.1 | 1.9 | ±1.0 | **208.1 – 210.1** | 0.5% |
| Polar cap T, SH | K | 192.3 | 2.1 | ±1.1 | **191.3 – 193.4** | 0.6% |

`*` detrended σ (mass-flux trend significant at p < 0.001 over this window). The DJF and JJA
rows rest on 34 and 35 seasons respectively — one DJF season is lost at the record edge.

Moving the anchor from 1996–2014 to 1980–2014 tightens five of the six tolerances — mass
flux from ±0.110 to ±0.102, NH polar cap T from ±1.12 to ±0.97 K — because the longer window
has the smaller σ. w̄* loosens, from ±0.0064 to ±0.0068.

### 3.2 Variance

At 35 years the interannual σ ratio becomes a usable test. D7 retired it at n = 10, where
the window was 0.43–1.57 and nothing realistic could fail it. At 35 years against a 35-year
anchor:

| Metric | 95% acceptance window | Width |
|---|---|---|
| **Interannual σ ratio** (rollout σ ÷ CESM σ) | **0.66 – 1.34** | 0.67 |
| **Daily DJF σ ratio**, u at 60°N | **0.87 – 1.13** | 0.26 |

- Both windows are tighter than against the 19-year anchor (0.60–1.40 and 0.84–1.16),
  because the reference σ is now estimated from 35 years rather than 19. That is the main
  practical gain from the longer anchor.
- **The daily ratio is the sharper test, at 35 years as at shorter lengths.** The daily DJF
  series gives about 6.4 independent samples per winter (14-day decorrelation time), so it
  constrains the variance about 2.6× more tightly than the annual series. Report both.
- **Detrend both sides before forming the interannual ratio.** A 35-year rollout spans long
  enough for the forced trend to inflate its raw σ: over the pooled record the mass-flux σ
  is 54% larger raw than detrended. The anchor σ in the table above is already detrended
  where the trend is significant, and the rollout must be treated the same way.

> **Interpretation** — The daily ratio is the test that detects an emulator that has smoothed
> away its own weather. An interannual ratio formed without detrending both sides measures
> trend, not variability.

### 3.3 Tested on 1970–1994

CESM 1970–1994 scored as if it were a rollout — 24 annual years, 23–25 seasons depending on
the diagnostic. Page 3 of
[AIDE_WACCM_screening_1970-1995.pdf](AIDE_WACCM_screening_1970-1995.pdf) is this table as a
figure.

| Diagnostic | Test mean | Offset | Mean verdict | σ ratio | Variance verdict |
|---|---|---|---|---|---|
| Mass flux, 70 hPa | 9.272 | −1.43σ | **FAIL** | 0.94 | pass |
| w̄*, 10°S–10°N | 0.1881 | −0.61σ | **FAIL** | 0.91 | pass |
| Vortex NH, DJF | 24.7 | −0.19σ | PASS | 1.07 | pass |
| Vortex SH, JJA | −81.8 | −0.13σ | PASS | 0.83 | pass |
| Polar cap T, NH | 210.0 | +0.48σ | PASS | 0.79 | pass |
| Polar cap T, SH | 192.1 | −0.10σ | PASS | 0.87 | pass |

Results:

- Plus the SSW count: **15 major NH warmings in 1970–1994** against an expected 16.4
  [9, 25] at the anchor rate of 0.66/winter — **pass**.
- Four of six pass the mean test. The two failures are the upwelling pair, both low: mass
  flux at −1.43σ is nearly three times the ±0.50σ tolerance, and w̄* at −0.61σ is just
  outside it.
- All six pass the variance test, spanning 0.79–1.07 inside a 0.66–1.34 window.
- 15 of the test period's 24 years also lie inside the anchor.

> **Interpretation** — Both mean failures are the forced BDC trend: 1970–1994 predates most
> of the acceleration, so its mean circulation is weaker than the 1980–2014 anchor's. This is
> not a defect of the anchor choice — the same failure appears against the 1996–2014 anchor,
> and it is the D9 finding restated at a different length. The variance result is the
> transferable one: CESM's year-to-year spread is stable across the whole record even where
> its mean level is not, so a variance target set on one period transfers to another and a
> mean target on the upwelling diagnostics does not. Read section 3.3 as a consistency check,
> not a validation: anchor and test period are not independent, so the passes are weaker
> evidence than the failures. Overlapping windows make agreement easier, so a diagnostic that
> still fails under overlap is failing on a real offset.

### What else becomes testable at 35 years

| | At 35 years | Note |
|---|---|---|
| Major NH SSW count | expect 27.6, accept **18–38** | rate resolved to ±37% |
| R1 wave → vortex slope | −1.345 ± 0.34 | resolvable |
| R2 thermal wind slope | −2.194 ± 0.45 | resolvable |
| Forced trend, mass flux | +0.1945 ± 0.0706 per decade | **resolvable** |
| Forced trend, w̄* | +0.0058 ± 0.0048 per decade | **resolvable** (marginally) |
| Forced trend, vortex NH | +0.598 ± 1.843 per decade | not resolvable |
| Forced trend, vortex SH | +0.351 ± 0.802 per decade | not resolvable |
| Forced trend, polar cap T NH | −0.562 ± 0.632 per decade | not resolvable |
| Forced trend, polar cap T SH | +0.300 ± 0.657 per decade | not resolvable |

Slope standard errors are scaled from the 44-year fits in `02b_trends.json` as `n^{-3/2}`.

> **Interpretation** — The upwelling trend becoming resolvable is the scientific reason to
> prefer 35 years over 20. It is also the quantity that made the targets fail out of sample
> (D9), so a 35-year rollout is the first length at which that failure can be diagnosed
> rather than absorbed by period-matching. The vortex and temperature trends stay within noise at any
> rollout length this project will run; do not report them as passes.

---

## 4. Input requirements

What a model has to supply to be scored against the thresholds above. Everything in this
section is what `scripts/aide_val_common.py` reads and what the diagnostics consume; the
CESM column is the reference implementation, not a constraint on the candidate model's own
grid (§5, rule 2).

### 4.1 The five fields

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
  product where the eddy flux is expected reverses the sign of tropical w̄* (§4.5).

### 4.2 Temporal resolution

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

### 4.3 Vertical grid

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

### 4.4 Horizontal grid

- **Latitude must be global, −90° to +90°.** The mass-flux integral runs over every latitude
  where w̄* > 0 at 70 hPa, not over a prescribed tropical band.
- Latitude bands used: **10°S–10°N** (w̄*), **60–90°N** and **60–90°S** (polar cap T),
  **45–75°N** (R1 heat flux). Band means are **cos(lat)-weighted with inclusive endpoints**,
  not Gaussian-weighted.
- Point latitudes: **exactly ±60°**, reached by linear interpolation in latitude. CESM's
  0.9424° grid has no point at 60°.
- Zonal means only. No longitudinal information enters any diagnostic.
- CESM reference: 192 latitudes, 0.9424° spacing (f09).

> **Interpretation** — The grid above is what CESM supplies, not what a candidate model must
> match. Rule 2 in §5 scores on the emulator's own grid with CESM reduced onto it, so a
> coarser model grid is admissible; an incomplete latitude range is not, because it changes
> the mass-flux integral itself.

### 4.5 Deriving the five fields from a standard archive

Most archives carry ω, T and 3-D winds rather than a zonal-mean TEM tape.

| Available | Derive | Constraint |
|---|---|---|
| ω (Pa s⁻¹) | `w = −Hω/p`, H = 7000 m | Geometric vertical velocity is not a substitute: using it is an 11% error on all upwelling |
| T (K) | `θ = T (1000/p)^κ`, κ = 0.2857 | p in hPa |
| 3-D v and θ | `v'θ' = [vθ] − [v][θ]`, formed on the native 3-D grid **before** zonal averaging | Cannot be recovered from zonal means alone |
| T archived directly | use it for polar cap T | The reference path reconstructs T from θ; a directly archived T is equivalent |

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

### 4.6 Minimum request, in one line

- **Tier 1** — 5 years of daily, global, zonal-mean ū, v̄, w̄ (log-pressure), θ̄ and v'θ' on a
  common set of pressure levels bracketing 10, 50, 70 and 100 hPa.
- **Tier 2** — the same fields, 35 years.

---

## 5. Rules that apply at both tiers

These are not negotiable per tier. They are what makes a comparison meaningful.

1. **Pin the estimator.** Model and truth go through `aide_val_common.tem_residual`. The
   two standard routes to w̄* differ by 10.8% on identical data — larger than the tier-2
   mass-flux tolerance itself (D8).
2. **Score on the emulator's grid**, with CESM reduced onto it, not the reverse.
3. **Period-match, or detrend both sides.** The upwelling targets are period-matched, not
   absolute (D9). This binds harder at tier 2, where 35 years spans a real trend.
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
$PY scripts/14_evaluation_tiers.py        # reads 07 and 02b
                                          # writes output/14_evaluation_tiers.json
$PY scripts/15_screen_out_of_sample.py    # reads 07 and 14
                                          # writes output/15_screen_out_of_sample.json
                                          #   and docs/AIDE_WACCM_screening_1970-1995.pdf
```

The pipeline is those five scripts in that order, plus `aide_val_common.py` and
`report_layout.py`. `14` and `15` read only existing JSON and take seconds; if `output/` is
empty, `02`, `02b` and `07` rebuild it from the tape in about five minutes. See the pipeline
order in [../CLAUDE.md](../CLAUDE.md).

If any number in this document changes, it changes because `14` or `15` changed. Do not
edit the tables by hand.

---

## 7. Open for sign-off

Two choices here were made to produce a usable document and are not derived from the data.
Each is cheap to change.

**Settled, 2026-08-18:**

- **±3σ as the tier-1 gate**, to be revisited once real failure modes have been seen.
  Section 2 shows how it behaves out of sample.
- **Tier 2 anchored on 1980–2014 and tested on 1970–1994.** Section 3 states the two
  consequences: the anchor spans the 1995/96 restart, and it overlaps the test period by
  15 years, so section 3.3 is a consistency check rather than an independent test.

**Still open:**

1. **Tier 1 gates on individual years only**, with the 5-year mean advisory. The brief
   specified individual years; the mean check was added here and can be dropped or promoted
   to a gate.
2. **SSW count, mechanism slopes and trends are reported at tier 2** although the brief
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
threshold.** Anchor and test periods are split rather than pooled, even though pooling gives
a longer record: the data that sets a threshold cannot also test it. Sections 2 and 3.3 are
this rule applied to the two tiers.

**D13 · Reference σ comes from CESM, not from observations.** The right reference for
scoring an emulator *against this CESM run*, and the wrong one for any claim about realism —
see appendix C.

## Appendix C · Limitations that apply to every number here

- **Fixed SST.** No ENSO, so the interannual σ anchoring every tolerance is smaller than the
  real atmosphere's, and the tolerances derived from it are tighter.
- **Volcanic years.** 1970–1995 contains El Chichón and Pinatubo, part of why its upwelling
  variability is larger. The tier-1 anchor is volcanically quiet, so its σ is at the low end.
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
- Appendix A — where `0.5σ` and `1.96σ/√n` come from
- Appendix B — the decisions cited above by number
- Appendix C — the limitations that apply to every number here
- [AIDE_WACCM_screening_1970-1995.pdf](AIDE_WACCM_screening_1970-1995.pdf) — the figure for
  section 2
- `scripts/14_evaluation_tiers.py`, `scripts/15_screen_out_of_sample.py` — every number on
  this page
