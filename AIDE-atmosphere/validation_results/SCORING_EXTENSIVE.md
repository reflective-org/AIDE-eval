# Transport scorecard — extensive scoring

Summary over every result set registered in `sources/`. Full tables: the reports linked below. Method, contract and decisions: [docs/EVALUATION_PROTOCOL.md](../docs/EVALUATION_PROTOCOL.md).

## What is scored

- **extensive** (the full listing): the 34 metrics of the phase-29 scorecard, every one scored.
- A age of air: 4 scored — A1, A2, A3, A4
- B residual circulation: 7 scored — B1-100, B1-70, B1-50, B1-30, B1-10, B2-100, B2-70
- C polar vortex: 10 scored — C1-jet50_N, C1-jet50_N-lat, C1-jet10_N, C1-jet10_N-lat, C1-jet50_S, C1-jet50_S-lat, C1-jet10_S, C1-jet10_S-lat, C2-N, C2-S
- D variability: 13 scored — D1-N, D1-S, D2-N, D2-S, D3, D4-fwN, D4-fwS, D4-w, D5-N, D5-S, D6-N, D6-S, D7
- Score: rel = model / ref − 1, err = \|rel\| / 0.25; pass err ≤ 1, near ≤ 1.5, fail above. A2, A2-WACCM, C1 latitudes, C2 and D4 final warming keep absolute bounds; A4 against the midpoint 2.75 yr.
- Group score = mean err of its scored metrics; combined = mean of A, B, C; D separate.
- References follow each set's window (A, C, D); B is WACCM6 h6 1996–2014 for every set.
- w*: jcm-strat Ψ* route, log-p interpolated — interim, AIDE-eval issue #5. Group C input from daily means — AIDE-eval issue #7.

Ranking: the standard scoring ([SCORING_STANDARD.md](SCORING_STANDARD.md)).

## jcm-strat-AL7h-2006-2019 (2006–2019)

| | AL7h |
|---|---|
| A age of air | 0.93 |
| B residual circulation | 0.95 |
| C polar vortex | 0.36 |
| **combined A-C** | **0.74** |
| D variability | 1.10 |
| passed A-C | 18 / 21 |
| passed D | 6 / 13 |

- Fail: A2, B1-100, B1-50, D1-N, D1-S, D6-N, D6-S, D7
- Near: D2-N, D2-S
- Produced 2026-11-12. Report: [`jcm-strat-AL7h-2006-2019/scorecard_extensive.md`](jcm-strat-AL7h-2006-2019/scorecard_extensive.md); figures [`jcm-strat-AL7h-2006-2019/score_table_extensive.png`](jcm-strat-AL7h-2006-2019/score_table_extensive.png), [`jcm-strat-AL7h-2006-2019/heatmap_extensive.png`](jcm-strat-AL7h-2006-2019/heatmap_extensive.png)

## jcm-strat-phase29-2006-2010 (2006–2010) — **non-standard window (5 of 14 standard years)**

| | P20b | AL7 | AL7t | L0 | LAL7 | QBO10d | LQBO10d |
|---|---|---|---|---|---|---|---|
| A age of air | 0.58 | 1.17 | 1.21 | 0.99 | 1.70 | 0.91 | 1.44 |
| B residual circulation | 1.05 | 1.00 | 1.01 | 0.98 | 1.06 | 0.86 | 0.78 |
| C polar vortex | 1.09 | 0.26 | 0.21 | 1.02 | 0.45 | 0.79 | 0.47 |
| **combined A-C** | **0.91** | **0.81** | **0.81** | **1.00** | **1.07** | **0.85** | **0.90** |
| D variability | 1.35 | 1.28 | 1.29 | 1.17 | 1.45 | 0.93 | 1.39 |
| passed A-C | 13 / 21 | 17 / 21 | 15 / 21 | 14 / 21 | 15 / 21 | 15 / 21 | 14 / 21 |
| passed D | 4 / 13 | 6 / 13 | 5 / 13 | 6 / 13 | 5 / 13 | 7 / 13 | 4 / 13 |

- Fail in every run: B1-100, D7
- Pass in every run: B1-70, B1-10, B2-100, B2-70, C1-jet50_N, C1-jet50_S, C1-jet50_S-lat, C1-jet10_S, C1-jet10_S-lat, D5-N, D5-S
- Produced 2026-11-12. Report: [`jcm-strat-phase29-2006-2010/scorecard_extensive.md`](jcm-strat-phase29-2006-2010/scorecard_extensive.md); figures [`jcm-strat-phase29-2006-2010/score_table_extensive.png`](jcm-strat-phase29-2006-2010/score_table_extensive.png), [`jcm-strat-phase29-2006-2010/heatmap_extensive.png`](jcm-strat-phase29-2006-2010/heatmap_extensive.png)
