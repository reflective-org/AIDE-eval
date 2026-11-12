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

The **standard evaluation of stratospheric dynamics** for AIDE models: one scorecard of age
of air (A), residual circulation (B), polar vortex (C) and variability (D), each metric
scored as its relative deviation from CLaMS, WACCM6 or ERA5 for the same years. It is the
jcm-strat relative scorecard (phase 29) with phase 31's scored set as the headline,
ported here; the diagnostics are recomputed from daily zonal means.

- Protocol: `docs/EVALUATION_PROTOCOL.md` (window, contract, metrics, references, scoring,
  open issues, how to evaluate a model, results, limitations).
- Code and data flow: `docs/ARCHITECTURE.md`. Decisions: `docs/DECISIONS.md` (settled and
  unconfirmed). Deferred and open: `docs/DEFERRED.md`. Status: `docs/PROGRESS.md`.
- **Standard window 2006–2019.** Any other window is scored against references of its own
  years and flagged "non-standard window (N of 14 standard years)".
- **Two scorings**: **standard** (= jcm-strat phase 31's set, the former "reduced"; the
  headline and ranking) and **extensive** (the 34 phase-29 metrics, the full listing).
- Scored: `jcm-strat-AL7h-2006-2019` (standard window) and `jcm-strat-phase29-2006-2010`
  (seven runs, non-standard). No AIDE-WACCM emulator rollout scored yet.
- The two-tier CESM-anchored protocol this repo held until 2026-11-12 is superseded; do not
  re-derive or cite its thresholds.

## Repository layout

Only `README.md`, `CLAUDE.md`, `LICENSE`, `requirements.txt`, `.gitignore` and the
gitignored `.AIDE-eval_env/` sit at the repo root. Everything else is under
**`AIDE-atmosphere/`**; bare paths below are relative to it.

| Path | Content |
|---|---|
| `sources/<set>.json` | one result set: description, reader, window, runs (path or list of paths), optional notes |
| `scripts/strat_eval_common.py` | shared constants, `REFS`, h6 loader, contract, set registry, readers, scoring, TEM math |
| `scripts/01_fields.py`, `02_scorecard.py`, `03_figures.py` | measure → score → draw; `evaluate.py <set> [...]` runs all three |
| `scripts/checks/port_check.py` | evidence against jcm-strat phases 29 and 31; outside the chain |
| `scripts/copy_era5_refs.sh`, `era5_refs.sha256` | ERA5 copies into `data/era5/` (gitignored) and their manifest |
| `output/<set>/` | per-run and reference JSON — generated, gitignored |
| `validation_results/` | `SCORING_STANDARD.md`, `SCORING_EXTENSIVE.md`, `port_check.md`, `<set>/` — generated, committed |

## Data and environment

| | |
|---|---|
| Interpreter | `.AIDE-eval_env/bin/python` at the repo root, from `requirements.txt` (Python 3.10.12) |
| Model input | paths in `sources/<set>.json`, read in place (jcm-strat files under `/data/JCM_stripped/jcm-strat-phase29`, `-phase31`; never modify them) |
| References | CLaMS `/data/CLaMS/CLaMS_v3/clams_v3.1_era5_zm_lat.zip`; WACCM6 entry age `/data/CESM2_REFD1_AOA/`; CESM2.1.5-WACCM6 h6 `/data/cesm2.1.5_output/histSST` (read-only); ERA5 copies in `data/era5/` |
| Memory | wrap heavy runs in `/data/JCM_stripped/jcm-strat-phase29/scripts/memcap.sh` (96 G per job); one 14-year run peaks at 2.6 GB |

Reproduce (from `scripts/`, `PY=../../.AIDE-eval_env/bin/python`, about 80 s):

```bash
bash copy_era5_refs.sh
$PY evaluate.py && $PY checks/port_check.py
```

## Data conventions that are already settled — do not re-derive, do not violate

Established against the CESM h6 files; each changes the answer by 10–1000 % if assumed
wrong. Evidence: scripts `01*` and their logs, archived with the tier protocol as git tag
`tier-protocol-final`.

1. **`VTHzm` is already the eddy flux `v'θ'`.** Do NOT subtract `Vzm·THzm`. Verified against
   raw 3-D h7 fields: 15 % median error as-is, 1087 % if subtracted. Same for `UVzm`, `UWzm`.
   The field contract requires the eddy flux from every model.
2. **`Wzm` is the log-pressure vertical velocity** `w = −Hω/p`, H = 7 km — not geometric.
3. **Below-surface points are the sentinel `1e35` with no `_FillValue`.** Mask with
   `abs(x) < 1e20` (`load_h6`, the `netcdf` reader).
4. **`MSKtem` is not a 0/1 mask.** It is a `(time, lat, lon)` fractional count of
   above-surface interfaces and cannot be applied to the zonal-mean fields.

## Scientific rules that constrain any change

- **One estimator for model and reference.** Every w* and Ψ* — model, WACCM6 and ERA5 — goes
  through `tem_streamfunction` / `residual_w` with log-p interpolation. This estimator is
  **interim** (issue #5); do not mix it with another route.
- **Reference years follow the scored window** (A, C, D). B is WACCM6 h6 1996–2014 for every
  window and is not period-matched; say so with any B number.
- **State the window and its flag with every score.** Short windows make D1 and D3 sampling
  errors comparable to the 25 % tolerance.
- **Changing a metric, tolerance, reference or scoring changes what a score means**: record it
  in `docs/DECISIONS.md`, regenerate every set, and update the protocol in the same commit.
- **Group C input is daily means** (issue #7); do not compare C2 with jcm-strat below 0.02 K.
- The port check must keep reproducing jcm-strat (`validation_results/port_check.md`); a code
  change that moves its agreement counts is a regression unless explained.

## Working rules specific to this repo

- `output/` is **generated** and gitignored; `validation_results/` is **generated and
  committed**. Never hand-edit a JSON, report or figure — change the script and re-run
  `evaluate.py` (and `checks/port_check.py`), so the artefact matches the code that made it.
- The dependency pins in `requirements.txt` are part of the result. The chain reproduces all
  23 committed artefacts bit for bit under them (checked 2026-11-12, two runs). Changing a
  pin means re-running that check.
- Numbers in `docs/EVALUATION_PROTOCOL.md` §4 and §10 are transcribed from
  `output/<set>/references.json` and `validation_results/<set>/scorecard_<scoring>.json`. If a
  script changes a number, update the markdown in the same commit and check the
  transcription.
- A new model is one `sources/<set>.json` and `evaluate.py <set>`; no code change.
- Known limitations are in `docs/EVALUATION_PROTOCOL.md` §11 — read them before treating a
  score as fixed.
