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
| Both tiers tried on CESM 1970–1995 | [docs/AIDE_WACCM_screening_1970-1995.pdf](docs/AIDE_WACCM_screening_1970-1995.pdf) (3 pp) |
| Where the tolerances come from | [docs/EVALUATION_PROTOCOL.md](docs/EVALUATION_PROTOCOL.md) appendix A |
| Why a choice was made | [docs/EVALUATION_PROTOCOL.md](docs/EVALUATION_PROTOCOL.md) appendix B |
| What is knowingly not covered | [docs/EVALUATION_PROTOCOL.md](docs/EVALUATION_PROTOCOL.md) appendix C, §7 |
| Settled data conventions and working rules | [CLAUDE.md](CLAUDE.md) |
| The evidence for those conventions | `scripts/01*`, `logs/01*` |
| The 10/20-year suite the protocol replaced | [docs/AIDE_WACCM_validation_targets.pdf](docs/AIDE_WACCM_validation_targets.pdf) (23 pp, frozen — see History) |

The protocol document is self-contained: every threshold, its derivation, the decisions
behind it and the limitations on it are in that one file, and every number in it is
transcribed from `output/14_evaluation_tiers.json` or `output/15_screen_out_of_sample.json`.

## The two tiers

**Tier 1 — 5-year screening.** Every individual year against a ±3σ band on the CESM
1996–2014 mean, six diagnostics, plus an SSW count. A regression test: fast, quiet on a
healthy model, loud on a broken one. Run against CESM 1970–1995 as five consecutive
screening runs, 28 of 30 block verdicts pass; the two failures are both mass flux, both
before 1980, and both are the forced BDC trend rather than a defect.

**Tier 2 — 35-year validation.** The rollout mean to ±0.5σ and the variance to a
0.66–1.34 interannual σ ratio, anchored on 1980–2014. This is the actual scientific claim.
Scored on CESM 1970–1994, four of six diagnostics pass the mean test and all six pass the
variance test — the two mean failures are the upwelling pair, for the same trend reason.

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
The two standard routes to w̄* differ by **10.8%** on identical data, against a 1.1%
mass-flux tolerance at tier 2. A threshold is therefore a statement about a *specific
estimator*, valid only when model and truth pass through the same code path
(`scripts/aide_val_common.py:tem_residual`) on the same grid.

## Running it

Inputs are read-only CESM history files under `/data/cesm2.1.5_output/histSST`
(daily zonal-mean TEM tape `cam.h6`, 1970–2014).

```bash
cd scripts
PY=/home/ubuntu/atmospheric_scale/paradis_model/paradis_venv/bin/python
$PY 02_reference_stats.py && $PY 02b_trends.py && $PY 07_period_split.py \
  && $PY 14_evaluation_tiers.py && $PY 15_screen_out_of_sample.py
```

About five minutes in total. Order matters — later scripts read earlier JSON. `14` and `15`
read only existing JSON and take seconds, so iterating on a threshold does not mean
re-reading the tape.

`01_check_conventions.py`, `01b`, `01c` and `01d` sit outside the pipeline. They establish
the four data conventions in [CLAUDE.md](CLAUDE.md) — that `VTHzm` is already an eddy flux,
that `Wzm` is log-pressure, the `1e35` sentinel, and `MSKtem` — each worth 10–1000% if
assumed wrong, and each silent. Run them standalone against the tape whenever a convention
is in doubt; they feed nothing downstream.

Requires numpy, xarray, scipy, matplotlib and cftime. No LaTeX or reportlab —
the PDF is built with matplotlib's `PdfPages`.

## Repository layout

```
README.md            this file
CLAUDE.md            working rules, settled data conventions, protocol constraints
scripts/             02, 02b, 07, 14, 15 + aide_val_common.py, report_layout.py
                     01, 01b, 01c, 01d — convention evidence, outside the pipeline
docs/                the protocol, its figure, and the frozen 23 pp report
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

Two items were kept:

- **`scripts/01*`**, the convention evidence. Re-run on 2026-08-19; they reproduce the
  numbers quoted in [CLAUDE.md](CLAUDE.md) — 15.4% median error for `VTHzm` as-is against
  1086.9% if `Vzm·THzm` is subtracted, a continuity slope of 1.011 at r = 0.9999 confirming
  log-pressure `Wzm`, and 2.9% covariance loss from daily averaging.
- **`docs/AIDE_WACCM_validation_targets.pdf`**, the 23 pp technical report and the first PDF
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
