# AIDE-WACCM atmosphere validation

The **two-tier evaluation protocol** an AIDE-WACCM emulator rollout is scored against for
tropical upwelling (Brewer–Dobson circulation) and the stratospheric polar vortex,
measured from CESM2.1.5-WACCM6 fixed-SST output.

This repository derives thresholds. It is not a model, and it does not score one:
every number in it is CESM-vs-CESM, including the dry run of both tiers against CESM's own
earlier output. Scoring an actual emulator rollout is the next piece of work and has never
been done.

## Start here

| If you want | Read |
|---|---|
| The protocol — both tiers, every threshold, the reasoning | [docs/EVALUATION_PROTOCOL.md](AIDE-atmosphere/docs/EVALUATION_PROTOCOL.md) |
| A scored climate model, with figures | [validation_results/](AIDE-atmosphere/validation_results/) |
| How to score your own model | [docs/EVALUATION_PROTOCOL.md](AIDE-atmosphere/docs/EVALUATION_PROTOCOL.md) §3.5, and §4 for the inputs |
| Where the tolerances come from | [docs/EVALUATION_PROTOCOL.md](AIDE-atmosphere/docs/EVALUATION_PROTOCOL.md) appendix A |
| Why a choice was made | [docs/EVALUATION_PROTOCOL.md](AIDE-atmosphere/docs/EVALUATION_PROTOCOL.md) appendix B |
| What is knowingly not covered | [docs/EVALUATION_PROTOCOL.md](AIDE-atmosphere/docs/EVALUATION_PROTOCOL.md) appendix C, §7 |
| Settled data conventions and working rules | [CLAUDE.md](CLAUDE.md) |
| The evidence for those conventions | `scripts/01*`, `logs/01*` |
| Superseded material, kept | [stale/](AIDE-atmosphere/stale/README.md) |

The protocol document is self-contained: every threshold, its derivation, the decisions
behind it and the limitations on it are in that one file, and every number in it is
transcribed from `output/16_anchors_45yr.json`.

## The two tiers

Both tiers are anchored on **CESM 1970–2014, the whole record** — 44 annual years, 42 DJF
and 45 JJA seasons — so the thresholds are fixed and model-agnostic rather than matched to a
CESM sub-period.

**Tier 1 — 5-year screening.** Every individual year against a ±3σ band on the anchor mean,
six diagnostics, plus an SSW count of 0–7 over 5 winters. A regression test: fast, quiet on
a healthy model, loud on a broken one.

**Tier 2 — 35-year validation.** The rollout mean to ±0.5σ, the interannual σ ratio to
0.68–1.32 and the daily DJF σ ratio to 0.88–1.12, plus the SSW count, the two mechanism
slopes and the upwelling trend. This is the scientific claim.

**What the anchor gives up.** Spanning the whole record leaves no CESM output held back to
test the thresholds against, so §2 of the protocol cites the archived split-anchor runs in
[stale/](AIDE-atmosphere/stale/README.md) instead: on anchors that did not include the scored period, 28 of
30 tier-1 block verdicts and 4 of 6 tier-2 mean tests pass, the failures being the upwelling
pair and the forced BDC trend behind them. The method carries over; the numbers do not.

## Input requirements

What a model has to supply to be scored against the thresholds above. Everything in this
section is what `scripts/aide_val_common.py` reads and what the diagnostics consume. The
CESM values quoted throughout are the reference implementation, not a constraint on the
candidate model's own grid (§5, rule 2).

**Tier 1 needs 5 years, tier 2 needs 30 years.** Both need the same fields, at the same
resolution, on the same grid. Two forms of the request, depending on what the archive holds:

| | **A zonal-mean TEM tape** | **A standard pressure-level archive** |
|---|---|---|
| Fields | ū, v̄, w̄ (log-pressure), θ̄, v'θ' | `ua`, `va`, `wap`, `ta` (CF/CMIP names) |
| Shape | (time, level, lat), already zonally averaged | (time, level, lat, lon) daily 3-D — or zonal means **plus** 3-D `va` and `ta` |
| Temporal | daily means | daily means |
| Levels | pressure levels bracketing 10, 50, 70 and 100 hPa, with one level beyond each end | same |
| Latitude | global, −90° to +90° | same |
| Derivation | none | w̄ from `wap`, θ̄ from `ta`, v'θ' from 3-D `va` and `ta` (§4.6) |


## Running it

Inputs are read-only CESM history files under `/data/cesm2.1.5_output/histSST`
(daily zonal-mean TEM tape `cam.h6`, 1970–2014).

Set the environment up once, from the repo root:

```bash
python3 -m venv .AIDE-eval_env
.AIDE-eval_env/bin/pip install -r requirements.txt
```

Then:

```bash
cd AIDE-atmosphere/scripts
PY=../../.AIDE-eval_env/bin/python
$PY 02_reference_stats.py && $PY 02b_trends.py && $PY 07_period_split.py \
  && $PY 16_anchors_45yr.py && $PY 17_validate.py 1996 2014 \
  && $PY 18_validation_figures.py
```

About five minutes in total. Order matters — later scripts read earlier JSON. `16`, `17` and
`18` read only existing JSON and take seconds, so iterating on a threshold or scoring another
climate model does not mean re-reading the tape.

`17_validate.py` takes `--climate-model NAME` and the period as arguments, and stamps every
artefact it writes `__<climate model>__<production date>`. Climate-model data enters through
one function, `climate_model_series`, which is what gets replaced when the model is a
rollout rather than a window of CESM's own record.

`01_check_conventions.py`, `01b`, `01c` and `01d` sit outside the pipeline. They establish
the four data conventions in [CLAUDE.md](CLAUDE.md) — that `VTHzm` is already an eddy flux,
that `Wzm` is log-pressure, the `1e35` sentinel, and `MSKtem` — each worth 10–1000% if
assumed wrong, and each silent. Run them standalone against the tape whenever a convention
is in doubt; they feed nothing downstream.

`requirements.txt` pins numpy, scipy, xarray, matplotlib, cftime and netCDF4 to the
versions that produced every committed number, on Python 3.10.12. No LaTeX and no reportlab —
the figures are matplotlib. A change to those pins is a change to the results: the whole
pipeline reproduces all 15 artefacts bit-for-bit under them, and a new environment should be
checked the same way before its numbers are trusted.

## Repository layout

```
README.md              this file
CLAUDE.md              working rules, settled data conventions, protocol constraints
LICENSE
requirements.txt       pinned dependencies
.gitignore
.AIDE-eval_env/        the virtual environment — gitignored, 425 MB, rebuild from the above

AIDE-atmosphere/       everything else lives here
  scripts/             02, 02b, 07, 16, 17, 18 + aide_val_common.py, report_layout.py
                       01, 01b, 01c, 01d — convention evidence, outside the pipeline
  docs/                the protocol
  validation_results/  a scored climate model: validation_result__<stamp>.md and five figures
  stale/               superseded material, kept — see stale/README.md
  output/              JSON results            — generated, gitignored
  logs/                stdout of every run     — generated, gitignored
```

Paths written bare below — `scripts/`, `output/`, `docs/` — are relative to
`AIDE-atmosphere/`.

`output/` and `logs/` are not in version control: they are reproducible from the scripts
and `/data`, and JSON diffs would dominate the history. Never hand-edit a generated
artefact — change the script and re-run, so the result always matches the code that
made it.
