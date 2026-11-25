# Progress — AIDE-SAI validation

## Status (2026-11-25)

| Item | State |
|---|---|
| Plan | step 1 approved: CESM, 2030s, 3 years ([VALIDATION_PLAN.md](VALIDATION_PLAN.md)) |
| Injection rates | digitised, [`tilmes2025_fig2b.csv`](../validation_results/injection_rates/tilmes2025_fig2b.csv) |
| AIDE-SAI schedule support | PR #9 merged (891f782) |
| AIDE-SAI corrected rates | PR #12 open (3.81, 4.72, 5.60 Tg SO2/yr for 2030–2032) |
| `ccmi2030_inj_tick` (comparison run) | started 2026-11-25 ~22:15 UTC, GPU 6 |
| `ccmi2030_inj` (+10 % rate case) | running; step 2308/4380 at 2026-11-25 20:10 UTC; rate step at h 8760 seen |
| `ccmi2030_ctrl` | running; step 2308/4380 at 2026-11-25 20:10 UTC; resume seam at step 2165 (below) |
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

## Known issues in the runs

- **`ccmi2030_ctrl` has a resume seam at step 2165 (h 12 990).** The resume rebuilt the
  advection boundary values (`qfroz`) from hour 0 (January), so the redone step differs from
  the original: N/N0 12.47 → 11.84, M/M0 0.5289 → 0.5271, SO2 burden −4.6 %. Fixed for later
  resumes in AIDE-SAI PR #11 (4ac8d61). Kept by decision (2026-11-25); compare across this
  hour with care.
