# Progress — AIDE-SAI validation

## Status (2026-11-26)

| Item | State |
|---|---|
| Plan | step 1 approved: CESM, 2030s, 3 years ([VALIDATION_PLAN.md](VALIDATION_PLAN.md)) |
| Injection rates | digitised, [`tilmes2025_fig2b.csv`](../validation_results/injection_rates/tilmes2025_fig2b.csv) |
| AIDE-SAI schedule support | PR #9 merged (891f782) |
| AIDE-SAI schedules | PR #12 open: corrected 2030–2034 rates; 2027–2034 file with spin-up |
| `ccmi2027_inj` (comparison run, 2027–2032) | started 2026-11-26 00:37 UTC, GPU 5 (tmux `train_cesm_ccmi2027_inj`) |
| CCMI reference reader, REMAP | not started (Task 2) |
| Comparison | not started (Task 4) |

## Timeline

| Date | Where | Change |
|---|---|---|
| 2026-11-24 | AIDE-eval, draft plan | Reference, comparability limits, tasks |
| 2026-11-24 | AIDE-SAI PR #9 | Yearly SO2 injection schedule (`INJ_SO2_SCHED_FILE`) and the CCMI 2030–2034 file |
| 2026-11-24 | AIDE-SAI `runs/` | `ccmi2030_inj` and `ccmi2030_ctrl` started (tmux `train_cesm_ccmi2030_*`) |
| 2026-11-25 | AIDE-SAI `runs/` | `ccmi2030_ctrl` interrupted after step 2165 and resumed from its checkpoint under 68e6a66, on GPU 1 (was GPU 0) |
| 2026-11-25 | AIDE-eval PR 1 | Scaffold, docs, Fig. 2b fetch and digitiser |
| 2026-11-25 | AIDE-eval PR 1 review | Fig. 2b read at the year vertex (first reading ~10 % high); WACCM6 entry-age path fixed |
| 2026-11-25 | AIDE-SAI PR #12, `runs/` | Corrected schedule; `ccmi2030_inj_tick` started |
| 2026-11-26 | AIDE-SAI PR #12, `runs/` | `ccmi2030_*` runs (no spin-up) stopped and deleted; 2027–2034 schedule; `ccmi2027_inj` started |
