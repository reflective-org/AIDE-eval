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
| The protocol — both tiers, every threshold, the reasoning | [docs/EVALUATION_PROTOCOL.md](docs/EVALUATION_PROTOCOL.md) |
| A scored candidate, with figures | [validation_results/](validation_results/validation_result.md) |
| How to score your own model | [docs/EVALUATION_PROTOCOL.md](docs/EVALUATION_PROTOCOL.md) §3.5, and §4 for the inputs |
| Where the tolerances come from | [docs/EVALUATION_PROTOCOL.md](docs/EVALUATION_PROTOCOL.md) appendix A |
| Why a choice was made | [docs/EVALUATION_PROTOCOL.md](docs/EVALUATION_PROTOCOL.md) appendix B |
| What is knowingly not covered | [docs/EVALUATION_PROTOCOL.md](docs/EVALUATION_PROTOCOL.md) appendix C, §7 |
| Settled data conventions and working rules | [CLAUDE.md](CLAUDE.md) |
| The evidence for those conventions | `scripts/01*`, `logs/01*` |
| Superseded material, kept | [stale/](stale/README.md) |

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
[stale/](stale/README.md) instead: on anchors that did not include the scored period, 28 of
30 tier-1 block verdicts and 4 of 6 tier-2 mean tests pass, the failures being the upwelling
pair and the forced BDC trend behind them. The method carries over; the numbers do not.

## The three results worth knowing

**1. Thresholds depend on rollout length, and two roles of σ get conflated.**
σ as a statement about nature comes from the longest record available; σ as a limit on
detection is σ/√n. The operative rule is `target(n) = max(0.5σ, 1.96σ/√n)`, whose two
branches cross at **n = 15.4 years**. That crossover is why tier 1's 5-year mean check is
advisory only, and why tier 2 needs 35 years to make the variance and the forced trend
testable at all.

**2. A threshold CESM cannot meet on another sample of its own output is not a threshold.**
Both tiers are therefore tried out of sample before being published. Both upwelling
diagnostics fail on 1970–1994 because the BDC acceleration is concentrated in the earlier
period, which makes it a protocol requirement rather than a looser number: **upwelling
targets are period-matched, not absolute** — score a rollout against the CESM years it
covers, or detrend both sides.

**3. The estimator spread is larger than the target.**
The two standard routes to w̄* differ by **10.8%** on identical data, against the 1.1%
mass-flux tolerance at tier 2. A threshold is therefore a statement about a *specific
estimator*, valid only when model and truth pass through the same code path
(`scripts/aide_val_common.py:tem_residual`) on the same grid.

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
cd scripts
PY=../.AIDE-eval_env/bin/python
$PY 02_reference_stats.py && $PY 02b_trends.py && $PY 07_period_split.py \
  && $PY 16_anchors_45yr.py && $PY 17_validate.py 1996 2014 \
  && $PY 18_validation_figures.py
```

About five minutes in total. Order matters — later scripts read earlier JSON. `16`, `17` and
`18` read only existing JSON and take seconds, so iterating on a threshold or scoring another
candidate does not mean re-reading the tape.

`17_validate.py` takes the candidate period as arguments. Candidate data enters through one
function, `candidate_series`, which is what gets replaced when the candidate is a model
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
README.md            this file
CLAUDE.md            working rules, settled data conventions, protocol constraints
requirements.txt     pinned dependencies
.AIDE-eval_env/      the virtual environment — gitignored, 425 MB, rebuild from the above
scripts/             02, 02b, 07, 16, 17, 18 + aide_val_common.py, report_layout.py
                     01, 01b, 01c, 01d — convention evidence, outside the pipeline
docs/                the protocol
validation_results/  a scored candidate: validation_result.md and five figures
stale/               superseded material, kept — see stale/README.md
output/              JSON results            — generated, gitignored
logs/                stdout of every run     — generated, gitignored
```

`output/` and `logs/` are not in version control: they are reproducible from the scripts
and `/data`, and JSON diffs would dominate the history. Never hand-edit a generated
artefact — change the script and re-run, so the result always matches the code that
made it.

## History

The protocol supersedes an earlier 10/20-year target suite: a target table and CSV, the
observability and mechanism diagnostics, a plain-language explainer and a 10-year scorecard,
and separate decision, concepts, tolerance-rule and deferred-work documents. That suite was
removed on 2026-08-19 and **is not retained here** — what it established survives only where
this protocol restates it, principally appendices A–C.

The seasonal-cycle and daily-DJF-distribution plots went with it and are to be rebuilt
against the tier protocol. The diagnostics they visualise are still computed, in
`output/07_period_split.json`.

Three items were kept:

- **`scripts/01*`**, the convention evidence. Re-run on 2026-08-19; they reproduce the
  numbers quoted in [CLAUDE.md](CLAUDE.md) — 15.4% median error for `VTHzm` as-is against
  1086.9% if `Vzm·THzm` is subtracted, a continuity slope of 1.011 at r = 0.9999 confirming
  log-pressure `Wzm`, and 2.9% covariance loss from daily averaging.
- **`stale/`**, everything superseded since: the split-anchor threshold scripts `14` and
  `15`, the screening figure they produced, and the 23 pp technical report. `14` and `15`
  still run and still reproduce their JSON, and §2 of the protocol cites them as the
  out-of-sample evidence the current anchor cannot re-earn.
- **`stale/AIDE_WACCM_validation_targets.pdf`**, the 23 pp technical report and the first PDF
  the project produced. It is **frozen**: the script that built it read JSON from five
  diagnostics that are no longer part of the repo, so it cannot be regenerated here and is
  kept as the only surviving record of the suite the protocol replaced. Its thresholds are
  superseded —
  read them as history, not as current targets, and note that the `max(0.5σ, 1.96σ/√n)`
  columns in it are evaluated at n = 10 and n = 20, not at the tier lengths of 5 and 35.

## Caveats you should not skip

Fixed SST means no ENSO, so the interannual σ anchoring every tolerance is smaller than the
real atmosphere's — the right reference for scoring against *this* CESM run, the wrong
one for claiming realism. The 1970–1995 and 1996–2014 segments are two separate runs
joined at a restart, not one integration, and the tier-2 anchor spans that restart. Full
list in [docs/EVALUATION_PROTOCOL.md](docs/EVALUATION_PROTOCOL.md) appendix C.
