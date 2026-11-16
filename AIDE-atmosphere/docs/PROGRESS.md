# Progress

What is scored, when, and through which PR. Repository `reflective-org/AIDE-eval`.

## Scored

| Result set | Runs | Window | Standard combined A–C | Produced |
|---|---|---|---|---|
| `jcm-strat-AL7h-2006-2019` | AL7h | 2006–2019 (standard) | 0.69 | 2026-11-12 |
| `jcm-strat-phase29-2006-2010` | P20b, AL7, AL7t, L0, LAL7, QBO10d, LQBO10d | 2006–2010 (non-standard, 5 of 14) | 0.71 (AL7t) … 0.99 (LAL7) | 2026-11-12 |

- Reports and ranking: [validation_results/SCORING_STANDARD.md](../validation_results/SCORING_STANDARD.md).
- No AIDE-WACCM emulator rollout scored yet.

## Timeline

| Date | PR / branch | Change |
|---|---|---|
| 2026-08-19 … 2026-08-26 | `reanalysis-validation`, PR #4 `jcm-rollout-validation` | Two-tier CESM protocol; four reanalyses and jcm-strat phase 10 scored against it (archived, tag `tier-protocol-final`) |
| 2026-11-10 | PR #6 `strat-transport-scoring` | Transport scorecard ported from jcm-strat phase 29: field contract, readers, groups A–D, port check, eight runs 2006–2010 |
| 2026-11-11 | PR #6 update | Reference years follow the window; extensive and reduced scorings; A2-WACCM; AL7h 2011–2019 |
| 2026-11-12 | PR #6 update | One-page summary per scoring |
| 2026-11-12 | PR #8 `strat-eval-standard` (stacked on PR #6) | Scorecard replaces the tier protocol; tier material deleted (tag); `strat_eval_common.py`; scripts 01–03, `checks/port_check.py`; `sources/<set>.json` and `evaluate.py`; standard / extensive; standard window 2006–2019 and non-standard flag; AL7h 2006–2019; seven-run set; per-set results; docs rewritten |

## PR and branch state (2026-11-14)

- PR #4 `jcm-rollout-validation` → `reanalysis-validation`: merged 2026-11-13.
- PR #6 `strat-transport-scoring` → `jcm-rollout-validation`: merged 2026-11-14.
- PR #8 `strat-eval-standard` → `strat-transport-scoring`: merged 2026-11-14.
- PR #9 `strat-eval-standard` → `main`: merged 2026-11-14. #4, #6, #8 were stacked and had
  not reached `main`; #9 brings all of it there. Decision status is unchanged by the merges
  ([DECISIONS.md](DECISIONS.md), "Unconfirmed — PR #6" and "Unconfirmed — PR #8").
- Branches deleted 2026-11-14 (local and GitHub): `publish-tier-protocol`,
  `feature-45yr-anchors`, `remove-trend-resolvability`, `reanalysis-validation`,
  `jcm-rollout-validation`, `strat-transport-scoring`. Every commit in them is in
  `strat-eval-standard`; only GitHub's merge commits of PRs #4, #6, #8 (no file changes)
  are not.
- Remaining branches: `main` (all current work), `dev` (PRs #1–#3, 2026-08), and
  `strat-eval-standard` (merged into `main`). How future AIDE components are organised in
  this repo (branches vs directories on `main`) is open, pending a team decision.
- Issues #5 (w* estimator) and #7 (group C input) open ([DEFERRED.md](DEFERRED.md)).

## Reproduction

- The chain (ERA5 manifest check, `evaluate.py`, `checks/port_check.py`) run twice from an
  empty `output/` and `validation_results/` on 2026-11-12: all 23 committed artefacts and
  the 10 output JSON identical (sha256). About 80 s on CPU.
- 2026-11-12, after consolidating into one checkout (`/data/AIDE-eval`, the
  worktree `/data/AIDE-eval-strat-transport` removed): two more runs from the main checkout
  reproduce every committed artefact bit for bit; result JSON carry repo-relative paths, so
  they no longer depend on the checkout location.
