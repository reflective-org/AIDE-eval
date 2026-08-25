1. Think Before Coding
Don't assume. Don't hide confusion. Surface tradeoffs.
Before implementing:
* State your assumptions explicitly. If uncertain, ask.
* If multiple interpretations exist, present them - don't pick silently.
* If a simpler approach exists, say so. Push back when warranted.
* If something is unclear, stop. Name what's confusing. Ask.
2. Simplicity First
Minimum code that solves the problem. Nothing speculative.
* No features beyond what was asked.
* No abstractions for single-use code.
* No "flexibility" or "configurability" that wasn't requested.
* No error handling for impossible scenarios.
* If you write 200 lines and it could be 50, rewrite it.
Ask yourself: "Would a senior engineer say this is overcomplicated?" If yes, simplify.
3. Surgical Changes
Touch only what you must. Clean up only your own mess.
When editing existing code:
* Don't "improve" adjacent code, comments, or formatting.
* Don't refactor things that aren't broken.
* Match existing style, even if you'd do it differently.
* If you notice unrelated dead code, mention it - don't delete it.
When your changes create orphans:
* Remove imports/variables/functions that YOUR changes made unused.
* Don't remove pre-existing dead code unless asked.
The test: Every changed line should trace directly to the user's request.
4. Goal-Driven Execution
Define success criteria. Loop until verified.
Transform tasks into verifiable goals:
* "Add validation" → "Write tests for invalid inputs, then make them pass"
* "Fix the bug" → "Write a test that reproduces it, then make it pass"
* "Refactor X" → "Ensure tests pass before and after"
For multi-step tasks, state a brief plan:
1. [Step] → verify: [check]
2. [Step] → verify: [check]
3. [Step] → verify: [check]


 1. You have to always make a plan and each plan needs to have a logical task and each task should have subtasks that are small enough for commits.
 2. Each big progress should be a PR that will be pushed to the repo for me to review
 3. do not make any assumptions without checking if it is ok. All assumptions must be checked with me.
 4. Please do not include Claude in the PRs and commits descriptions.
 5. Make various documentations from README, Architecture, Progress, Plans, Key decisions, Deffered, Features, etc.
 6. This is a scientific model and needs to be validated at any step. Visualizations are super helpful to understand issues and progress. Please make sure you make what make sense in terms of visulization.
 7. The rule is always asking questions even for small tweaks and decisions.
 8. The goal is scientific integrity, newer language, great documentation, collaboration, reproducibility, and transparency.
9. Be scientific in your writing, don't use weird sayings or adverbs
10. Be concise and prioritize bulletpoints over long text




---

# Project context — AIDE-atmosphere_validation

## What this repo is

A derivation suite for the **two-tier evaluation protocol**, not a model. It measures
CESM2.1.5-WACCM6 fixed-SST output and produces the thresholds an AIDE-WACCM emulator
rollout must meet for tropical upwelling and the polar vortex. Every number here is
CESM-vs-CESM; nothing has yet been run against an actual emulator rollout.

Authoritative outputs: `output/16_anchors_45yr.json` (the thresholds, both tiers anchored on
CESM 1970–2014) and `output/17_validation__<stamp>.json` (a scored climate model).
The protocol is `docs/EVALUATION_PROTOCOL.md`, and it is self-contained — thresholds,
derivation (appendix A), decisions (appendix B) and limitations (appendix C) in one file.
`validation_results/` holds a scored climate model as `validation_result__<stamp>.md` plus seven figures,
all generated from the same JSON, so they cannot drift apart.

The repo was trimmed to this protocol on 2026-08-19. Removed and **not retained**: the
target table and CSV, the observability and mechanism diagnostics, the explainer and 10-year
scorecard PDFs, and the separate decision/concepts/tolerance/deferred documents. Do not
re-derive any of it without asking.

Superseded but **kept**, in `stale/`: the split-anchor threshold scripts `14` and `15`, the
screening figure, and the 23 pp technical report. `14` and `15` still run
(`PYTHONPATH=../scripts`) and are cited by §2 of the protocol as the out-of-sample evidence
the 45-year anchor cannot re-earn. Do not treat any number in `stale/` as current, and do not
resurrect either script into the pipeline without asking.

Two things were restored after the trim:

- **`scripts/01*`** — the convention evidence (below). Runnable: they read only the tape.
- **`stale/AIDE_WACCM_validation_targets.pdf`** — the 23 pp technical report, the first PDF
  the project produced. **Frozen, not regenerable**: the script that built it needed JSON
  from five diagnostics no longer in the repo. It is a historical record of the 10/20-year
  target suite, not a live artefact, and its thresholds are superseded by the tier protocol.
  Do not cite its numbers as current; do not hand-edit it; this is the one artefact in the
  repo that does not match code that can rebuild it.

## Repository layout

Only `README.md`, `CLAUDE.md`, `LICENSE`, `requirements.txt`, `.gitignore` and the
gitignored `.AIDE-eval_env/` sit at the repo root. Everything else — `scripts/`, `docs/`,
`output/`, `logs/`, `stale/`, `validation_results/` — is under **`AIDE-atmosphere/`**, and
every bare path in this file is relative to it.

## Data and environment

| | |
|---|---|
| Data root | `/data/cesm2.1.5_output/histSST` (read-only input; never write here) |
| Stream | `cam.h6`, daily zonal-mean TEM tape: `Uzm Vzm Wzm THzm VTHzm UVzm UWzm` |
| Segments | 1970–1995 and 1996–2014 — two **separate runs** joined at a restart, not one integration |
| Interpreter | `.AIDE-eval_env/bin/python` at the repo root, built from `requirements.txt` (Python 3.10.12) |
| Output dir | `C.OUTDIR` in `scripts/aide_val_common.py`, derived from `__file__` as `<AIDE-atmosphere>/output`; `output/` and `logs/` are gitignored |

Reproduce in order (from `AIDE-atmosphere/scripts/`, `PY=../../.AIDE-eval_env/bin/python`, ~5 min total;
verified to rebuild every JSON, the report and all seven PNGs bit-for-bit from the tape):

```bash
$PY 02_reference_stats.py && $PY 02b_trends.py && $PY 07_period_split.py \
  && $PY 16_anchors_45yr.py && $PY 17_validate.py 1996 2014 \
  && $PY 18_validation_figures.py
```

`16`, `17` and `18` read only existing JSON and take seconds — iterate on a threshold, or
score another climate model, without re-reading the tape. `17` takes `--climate-model NAME`
and the period as arguments and stamps every artefact `__<climate model>__<production date>`;
its `climate_model_series` is the seam to replace for a model rollout. `01`, `01b`, `01c` and `01d` are evidence, not pipeline stages: they
feed nothing and are run standalone to re-confirm the conventions below.

## Data conventions that are already settled — do not re-derive, do not violate

Each was established against the files and each silently changes the answer by
10–1000% if assumed wrong. `aide_val_common.py` relies on all four. Evidence is in
`scripts/01*` and `logs/01*`; all four were re-run on 2026-08-19 and reproduce the numbers
quoted here.

1. **`VTHzm` is already the eddy flux `v'θ'`.** Do NOT subtract `Vzm·THzm`.
   Verified against raw 3-D h7 fields: 15% median error as-is, 1087% if subtracted.
   Subtracting it makes tropical w* come out downward. Same for `UVzm`, `UWzm`.
2. **`Wzm` is the log-pressure vertical velocity** `w = −Hω/p`, H = 7 km — not geometric.
   Continuity regression slope 1.011, r = 1.000. Getting this wrong is 11% on all upwelling.
3. **Below-surface points are the sentinel `1e35` with no `_FillValue`.** xarray will not
   mask them; an unguarded `.mean()` returns ~1e33. Mask with `abs(x) < 1e20`.
4. **`MSKtem` is not a 0/1 mask.** It is a `(time, lat, lon)` fractional count of
   above-surface interfaces and cannot be applied to the zonal-mean fields.

## Scientific protocol rules that constrain any new diagnostic

- **Pin the estimator.** The two standard routes to w̄* differ by **10.8%** on the same
  data — larger than the 1.1% tier-2 mass-flux tolerance. A target is a statement about a
  *specific estimator*, valid only when model and truth go through the identical code
  path: `aide_val_common.tem_residual`.
- **Score on the emulator's own grid**, with truth reduced onto that grid.
- **The upwelling targets are period-matched, not absolute.** They fail out of sample
  because the BDC trend is concentrated in 1970–1995. Compare a rollout against the CESM
  years it actually covers, or detrend both sides.
- **A target CESM itself cannot meet on a different sample of its own output is not a
  target.** Any new target gets the same train/test treatment as `07_period_split.py`.
- **The anchor has no held-out sample.** It spans 1970–2014, so no CESM output remains to
  test a new threshold against (D12). Any new diagnostic still owes that test — the archived
  split anchors in `stale/` are the template.
- **State the rollout length.** Thresholds are `max(0.5σ, 1.96σ/√n)`; the two branches
  cross at n = 15.4 yr. The interannual σ ratio and the forced BDC trend are **not
  testable** below ~15 yr — they are why tier 2 is 35 years, and they must not be reported
  as passes at tier-1 length.

## Working rules specific to this repo

- `output/` and `logs/` are **generated**. Never hand-edit a JSON, CSV or PDF —
  change the script and re-run it, so the artefact always matches the code that made it.
- The dependency pins in `requirements.txt` are part of the result. All 15 artefacts were
  verified bit-for-bit under them, and the same check was run against the older
  `/home/ubuntu/.../paradis_venv` interpreter — identical output. Changing a pin means
  re-running that check, not assuming it still holds.
- The report and figures are built from the JSON. If a number changes, re-run `17` and `18`
  in the same pass, or the tables and the figures disagree.
- Numbers in `docs/EVALUATION_PROTOCOL.md` §1 and §3 are transcribed from
  `output/16_anchors_45yr.json`. If a script changes a number, update the markdown in the
  same commit — and check the transcription, not just the prose.
- `validation_results/` is **generated and committed**, unlike `output/`. Never hand-edit the
  report or a figure; re-run `17` and `18`.
- Known limitations are in `docs/EVALUATION_PROTOCOL.md` appendix C — read it before
  treating any figure as fixed (fixed SST means no ENSO; the tier-2 anchor spans the
  1995/96 restart; SSW detection is a local Charlton–Polvani implementation, not a shared
  catalogue).
