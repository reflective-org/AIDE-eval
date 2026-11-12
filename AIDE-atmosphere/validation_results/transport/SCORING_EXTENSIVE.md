# Transport scorecard — extensive scoring

Summary over every result set, produced 2026-11-12. Full tables: the reports linked below. Method, contract and decisions: [docs/TRANSPORT_SCORECARD.md](../../docs/TRANSPORT_SCORECARD.md).

## What is scored

- **extensive**: the 34 metrics of the phase-29 scorecard, every one scored.
- A age of air: 4 scored — A1, A2, A3, A4
- B residual circulation: 7 scored — B1-100, B1-70, B1-50, B1-30, B1-10, B2-100, B2-70
- C polar vortex: 10 scored — C1-jet50_N, C1-jet50_N-lat, C1-jet10_N, C1-jet10_N-lat, C1-jet50_S, C1-jet50_S-lat, C1-jet10_S, C1-jet10_S-lat, C2-N, C2-S
- D variability: 13 scored — D1-N, D1-S, D2-N, D2-S, D3, D4-fwN, D4-fwS, D4-w, D5-N, D5-S, D6-N, D6-S, D7
- Score: rel = model / ref − 1, err = \|rel\| / 0.25; pass err ≤ 1, near ≤ 1.5, fail above. A2, A2-WACCM, C1 latitudes, C2 and D4 final warming keep absolute bounds; A4 against the midpoint 2.75 yr.
- Group score = mean err of its scored metrics; combined = mean of A, B, C; D separate.
- References follow each set's window (A, C, D); B is WACCM6 h6 1996–2014 for every set.
- w*: jcm-strat Ψ* route, log-p interpolated — interim, AIDE-eval issue #5. Group C input from daily means — AIDE-eval issue #7.

## jcm-strat-phase29-2006-2010 (2006–2010)

| | P20b | AL7 | AL7h | AL7t | L0 | LAL7 | QBO10d | LQBO10d |
|---|---|---|---|---|---|---|---|---|
| A age of air | 0.58 | 1.17 | 0.89 | 1.21 | 0.99 | 1.70 | 0.91 | 1.44 |
| B residual circulation | 1.05 | 1.00 | 0.96 | 1.01 | 0.98 | 1.06 | 0.86 | 0.78 |
| C polar vortex | 1.09 | 0.26 | 0.52 | 0.21 | 1.02 | 0.45 | 0.79 | 0.47 |
| **combined A-C** | **0.91** | **0.81** | **0.79** | **0.81** | **1.00** | **1.07** | **0.85** | **0.90** |
| D variability | 1.35 | 1.28 | 1.11 | 1.29 | 1.17 | 1.45 | 0.93 | 1.39 |
| passed A-C | 13 / 21 | 17 / 21 | 17 / 21 | 15 / 21 | 14 / 21 | 15 / 21 | 15 / 21 | 14 / 21 |
| passed D | 4 / 13 | 6 / 13 | 6 / 13 | 5 / 13 | 6 / 13 | 5 / 13 | 7 / 13 | 4 / 13 |

- Fail in every run: B1-100, D7
- Pass in every run: B1-70, B1-10, B2-100, B2-70, C1-jet50_N, C1-jet50_S, C1-jet50_S-lat, C1-jet10_S, C1-jet10_S-lat, D5-N, D5-S
- Report: [`transport_scorecard__jcm-strat-phase29-2006-2010__extensive__2026-11-12.md`](transport_scorecard__jcm-strat-phase29-2006-2010__extensive__2026-11-12.md); figures [`transport_score_table__jcm-strat-phase29-2006-2010__extensive__2026-11-12.png`](transport_score_table__jcm-strat-phase29-2006-2010__extensive__2026-11-12.png), [`transport_heatmap__jcm-strat-phase29-2006-2010__extensive__2026-11-12.png`](transport_heatmap__jcm-strat-phase29-2006-2010__extensive__2026-11-12.png)

## jcm-strat-phase31-AL7h-2011-2019 (2011–2019)

| | AL7h_1119 |
|---|---|
| A age of air | 0.95 |
| B residual circulation | 0.94 |
| C polar vortex | 0.30 |
| **combined A-C** | **0.73** |
| D variability | 1.28 |
| passed A-C | 18 / 21 |
| passed D | 5 / 13 |

- Fail: A2, B1-100, B1-50, D1-N, D1-S, D3, D6-N, D6-S, D7
- Near: D2-N, D2-S
- Report: [`transport_scorecard__jcm-strat-phase31-AL7h-2011-2019__extensive__2026-11-12.md`](transport_scorecard__jcm-strat-phase31-AL7h-2011-2019__extensive__2026-11-12.md); figures [`transport_score_table__jcm-strat-phase31-AL7h-2011-2019__extensive__2026-11-12.png`](transport_score_table__jcm-strat-phase31-AL7h-2011-2019__extensive__2026-11-12.png), [`transport_heatmap__jcm-strat-phase31-AL7h-2011-2019__extensive__2026-11-12.png`](transport_heatmap__jcm-strat-phase31-AL7h-2011-2019__extensive__2026-11-12.png)

### AL7h_1119 against AL7h: in sample 2006–2010 vs out of sample 2011–2019

| | 2006–2010 | 2011–2019 |
|---|---|---|
| A age of air | 0.89 | 0.95 |
| B residual circulation | 0.96 | 0.94 |
| C polar vortex | 0.52 | 0.30 |
| **combined A-C** | **0.79** | **0.73** |
| D variability | 1.11 | 1.28 |
| passed A-C | 17 / 21 | 18 / 21 |

- Scored metrics changing status: C1-jet10_N, D1-S, D2-N, D6-S
- Metric by metric: [`transport_scorecard__jcm-strat-phase31-AL7h-2011-2019__extensive__2026-11-12.md`](transport_scorecard__jcm-strat-phase31-AL7h-2011-2019__extensive__2026-11-12.md); figure [`transport_window_comparison__jcm-strat-phase31-AL7h-2011-2019__2026-11-12.png`](transport_window_comparison__jcm-strat-phase31-AL7h-2011-2019__2026-11-12.png)
