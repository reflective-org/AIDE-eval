# Progress — AIDE-SAI validation

## Status (2026-11-25)

| Item | State |
|---|---|
| Plan | step 1 approved: CESM, 2030s, 3 years ([VALIDATION_PLAN.md](VALIDATION_PLAN.md)) |
| Injection rates | digitised, [`tilmes2025_fig2b.csv`](../validation_results/injection_rates/tilmes2025_fig2b.csv) |
| AIDE-SAI schedule support | PR #9 merged (891f782) |
| `ccmi2030_inj`, `ccmi2030_ctrl` | running; step 2308/4380 at 2026-11-25 20:10 UTC; rate step at h 8760 seen |
| CCMI reference reader, REMAP | not started (Task 2) |
| Comparison | not started (Task 4) |

## Timeline

| Date | Where | Change |
|---|---|---|
| 2026-11-24 | AIDE-eval, draft plan | Reference, comparability limits, tasks |
| 2026-11-24 | AIDE-SAI PR #9 | Yearly SO2 injection schedule (`INJ_SO2_SCHED_FILE`) and the CCMI 2030–2034 file |
| 2026-11-24 | AIDE-SAI `runs/` | `ccmi2030_inj` and `ccmi2030_ctrl` started (tmux `train_cesm_ccmi2030_*`) |
| 2026-11-25 | AIDE-SAI `runs/` | `ccmi2030_ctrl` interrupted after step 2165 and resumed from its checkpoint under 68e6a66 |
| 2026-11-25 | AIDE-eval PR 1 | Scaffold, docs, Fig. 2b fetch and digitiser |
