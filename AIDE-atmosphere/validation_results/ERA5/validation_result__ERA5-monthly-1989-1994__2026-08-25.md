# Tier 1 (restricted) — ERA5 monthly 1989-1994

| | |
|---|---|
| Source | ERA5 monthly means via CDS (reanalysis-era5-pressure-levels-monthly-means) |
| Scored | 1989–1994 |
| Anchor | CESM 1970-2014 |
| Produced | 2026-08-25 |
| Scope | 4 of 6 tier-1 diagnostics |

## Scored

| Diagnostic | Unit | Accept | Seasons | Outside | Worst | 5-yr blocks | Verdict |
|---|---|---|---|---|---|---|---|
| Vortex NH, DJF | m s⁻¹ | 7.97 – 41.56 | 5 | 0 | +2.65σ (1993) | 0/1 | PASS |
| Vortex SH, JJA | m s⁻¹ | -88.95 – -74.20 | 6 | 4 | +6.50σ (1989) | 1/1 | **FAIL** |
| Polar cap T, NH | K | 203.84 – 215.25 | 5 | 0 | +2.24σ (1992) | 0/1 | PASS |
| Polar cap T, SH | K | 186.07 – 198.19 | 6 | 0 | +1.82σ (1992) | 0/1 | PASS |

## Shape — advisory only

A tier-2 shape check on a tier-1-length sample. Below the n = 15.4 crossover the tolerance is set by what the sample resolves rather than by what the physics tolerates, so a pass carries no information and is not counted. A failure would still count.

| Diagnostic | Amplitude offset | Tolerance | Phase offset | Tolerance | Status |
|---|---|---|---|---|---|
| Vortex NH, DJF | +6.056 (+2.31σ) | ±2.096 | -0.19 mo | ±0.20 mo | **FAIL** |
| Vortex SH, JJA | -1.131 (-0.56σ) | ±1.613 | -0.09 mo | ±0.18 mo | advisory |
| Polar cap T, NH | +1.774 (+1.53σ) | ±0.928 | -0.35 mo | ±0.11 mo | **FAIL** |
| Polar cap T, SH | +0.336 (+0.31σ) | ±0.875 | -0.21 mo | ±0.13 mo | **FAIL** |

## Not evaluable from a monthly archive

| Diagnostic | Why |
|---|---|
| Mass flux, 70 hPa | needs v'θ′, which cannot be formed after time averaging |
| w̄*, 10°S–10°N | needs v'θ′ |
| Major NH SSW count | needs the day the wind reverses |
| Daily DJF σ ratio | needs the daily series |

## Artefacts

| Figure | Status |
|---|---|
| tier1_screening | 4 of 6 diagnostics, gated |
| tier1_sigma | the same four, in anchor σ |
| shape_seasonal | 4 of 6, advisory |
| tier2_mean, tier2_variance | not produced — tier 2 needs 35 years |
| counts_and_relations | not produced — the SSW count needs daily data |
| shape_daily_distribution | not produced — needs the daily series |
| shape_w_star_profile | not produced — needs v'θ′ |

## Caveats

- ERA5 is a reanalysis of the real atmosphere; the anchor is fixed-SST CESM with no ENSO. Differences here are not by themselves an ERA5 error.
- The polar-cap layer uses 4 ERA5 levels inside 10–50 hPa against CESM's 8.
- Seasonal means from monthly data are exact, not approximate: every reduction involved is linear, and the day weighting is preserved.

