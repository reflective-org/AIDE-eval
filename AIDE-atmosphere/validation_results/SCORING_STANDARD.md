# Transport scorecard — standard scoring

Summary over every result set registered in `sources/`. Full tables: the reports linked below. Method, contract and decisions: [docs/EVALUATION_PROTOCOL.md](../docs/EVALUATION_PROTOCOL.md).

## What is scored

- **standard** (the headline score and ranking): jcm-strat phase 31's scored set, the former "reduced" scoring — A2-WACCM replaces A2; A2 and C1-jet10_N, C1-jet10_N-lat, C1-jet10_S, C1-jet10_S-lat, D1-N, D1-S, D2-N, D2-S, D3, D4-fwN, D4-fwS, D6-N, D6-S, D7 are listed but not scored. Phase 31 scored with absolute bounds; here the relative scoring applies to that set.
- A age of air: 4 scored — A1, A2-WACCM, A3, A4
- B residual circulation: 7 scored — B1-100, B1-70, B1-50, B1-30, B1-10, B2-100, B2-70
- C polar vortex: 6 scored — C1-jet50_N, C1-jet50_N-lat, C1-jet50_S, C1-jet50_S-lat, C2-N, C2-S
- D variability: 3 scored — D4-w, D5-N, D5-S
- Listed, not scored: A2, C1-jet10_N, C1-jet10_N-lat, C1-jet10_S, C1-jet10_S-lat, D1-N, D1-S, D2-N, D2-S, D3, D4-fwN, D4-fwS, D6-N, D6-S, D7
- Score: rel = model / ref − 1, err = \|rel\| / 0.25; pass err ≤ 1, near ≤ 1.5, fail above. A2, A2-WACCM, C1 latitudes, C2 and D4 final warming keep absolute bounds; A4 against the midpoint 2.75 yr.
- Group score = mean err of its scored metrics; combined = mean of A, B, C; D separate.
- References follow each set's window (A, C, D); B is WACCM6 h6 1996–2014 for every set.
- w*: jcm-strat Ψ* route, log-p interpolated — interim, AIDE-eval issue #5. Group C input from daily means — AIDE-eval issue #7.

## Ranking

Every run of every set by its combined A–C score (lower is better). Runs on different windows are compared with references of their own years.

| rank | run | set | window | combined A–C | A | B | C | D |
|---|---|---|---|---|---|---|---|---|
| 1 | AL7h | jcm-strat-AL7h-2006-2019 | 2006–2019 | **0.69** | 0.72 | 0.95 | 0.39 | 0.12 |
| 2 | AL7t | jcm-strat-phase29-2006-2010 | 2006–2010 — non-standard window (5 of 14 standard years) | **0.71** | 0.98 | 1.01 | 0.13 | 0.53 |
| 3 | QBO10d | jcm-strat-phase29-2006-2010 | 2006–2010 — non-standard window (5 of 14 standard years) | **0.73** | 0.71 | 0.86 | 0.63 | 0.28 |
| 4 | AL7 | jcm-strat-phase29-2006-2010 | 2006–2010 — non-standard window (5 of 14 standard years) | **0.75** | 0.95 | 1.00 | 0.30 | 0.42 |
| 5 | LQBO10d | jcm-strat-phase29-2006-2010 | 2006–2010 — non-standard window (5 of 14 standard years) | **0.87** | 1.20 | 0.78 | 0.63 | 0.73 |
| 6 | P20b | jcm-strat-phase29-2006-2010 | 2006–2010 — non-standard window (5 of 14 standard years) | **0.89** | 0.56 | 1.05 | 1.07 | 0.39 |
| 7 | L0 | jcm-strat-phase29-2006-2010 | 2006–2010 — non-standard window (5 of 14 standard years) | **0.93** | 0.79 | 0.98 | 1.03 | 0.54 |
| 8 | LAL7 | jcm-strat-phase29-2006-2010 | 2006–2010 — non-standard window (5 of 14 standard years) | **0.99** | 1.47 | 1.06 | 0.45 | 0.53 |

## jcm-strat-AL7h-2006-2019 (2006–2019)

| | AL7h |
|---|---|
| A age of air | 0.72 |
| B residual circulation | 0.95 |
| C polar vortex | 0.39 |
| **combined A-C** | **0.69** |
| D variability | 0.12 |
| passed A-C | 15 / 17 |
| passed D | 3 / 3 |

- Fail: B1-100, B1-50
- Near: none
- Produced 2026-11-12. Report: [`jcm-strat-AL7h-2006-2019/scorecard_standard.md`](jcm-strat-AL7h-2006-2019/scorecard_standard.md); figures [`jcm-strat-AL7h-2006-2019/score_table_standard.png`](jcm-strat-AL7h-2006-2019/score_table_standard.png), [`jcm-strat-AL7h-2006-2019/heatmap_standard.png`](jcm-strat-AL7h-2006-2019/heatmap_standard.png)

## jcm-strat-phase29-2006-2010 (2006–2010) — **non-standard window (5 of 14 standard years)**

| | P20b | AL7 | AL7t | L0 | LAL7 | QBO10d | LQBO10d |
|---|---|---|---|---|---|---|---|
| A age of air | 0.56 | 0.95 | 0.98 | 0.79 | 1.47 | 0.71 | 1.20 |
| B residual circulation | 1.05 | 1.00 | 1.01 | 0.98 | 1.06 | 0.86 | 0.78 |
| C polar vortex | 1.07 | 0.30 | 0.13 | 1.03 | 0.45 | 0.63 | 0.63 |
| **combined A-C** | **0.89** | **0.75** | **0.71** | **0.93** | **0.99** | **0.73** | **0.87** |
| D variability | 0.39 | 0.42 | 0.53 | 0.54 | 0.53 | 0.28 | 0.73 |
| passed A-C | 12 / 17 | 13 / 17 | 11 / 17 | 12 / 17 | 11 / 17 | 14 / 17 | 10 / 17 |
| passed D | 3 / 3 | 3 / 3 | 3 / 3 | 3 / 3 | 3 / 3 | 3 / 3 | 2 / 3 |

- Fail in every run: B1-100
- Pass in every run: B1-70, B1-10, B2-100, B2-70, C1-jet50_N, C1-jet50_S, C1-jet50_S-lat, D5-N, D5-S
- Produced 2026-11-12. Report: [`jcm-strat-phase29-2006-2010/scorecard_standard.md`](jcm-strat-phase29-2006-2010/scorecard_standard.md); figures [`jcm-strat-phase29-2006-2010/score_table_standard.png`](jcm-strat-phase29-2006-2010/score_table_standard.png), [`jcm-strat-phase29-2006-2010/heatmap_standard.png`](jcm-strat-phase29-2006-2010/heatmap_standard.png)
