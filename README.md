# AIDE-SAI validation
This repository validates single parts and the fully interactive version of the AIDE-SAI framework. If you want to run the fully interactive version or it's subparts go to:
- [Links to GitHub repos]
- [Links to GitHub repos]

## 1. AIDE-atmosphere validation

The **standard evaluation of stratospheric dynamics**: one scorecard of age of air,
residual circulation, polar vortex and variability, each metric scored as its relative
deviation from CLaMS, WACCM6 or ERA5 over the same years (tolerance ±25 %). Ported from the
jcm-strat relative scorecard; the diagnostics are recomputed here from daily zonal means
and reproduce jcm-strat's values (port check).

- **Standard window 2006–2019**, daily. Other windows are scored and flagged
  "non-standard window (N of 14 standard years)".
- **standard** scoring (jcm-strat phase 31's set, 20 metrics) is the headline score and
  ranking; **extensive** (34 metrics) is the full listing.
- Scored: jcm-strat AL7h 2006–2019 (standard combined A–C 0.69) and seven jcm-strat phase-29
  runs 2006–2010 (0.71–0.99, non-standard window). No AIDE-WACCM emulator rollout yet.
- Open: the w* estimator is interim (AIDE-eval #5); group C reads daily means (AIDE-eval #7).

## Start here

| If you want | Read |
|---|---|
| The protocol: window, inputs, metrics, references, scoring, limitations | [docs/EVALUATION_PROTOCOL.md](AIDE-atmosphere/docs/EVALUATION_PROTOCOL.md) |
| The scores and ranking | [validation_results/SCORING_STANDARD.md](AIDE-atmosphere/validation_results/SCORING_STANDARD.md), [SCORING_EXTENSIVE.md](AIDE-atmosphere/validation_results/SCORING_EXTENSIVE.md) |
| One set's report and figures | `AIDE-atmosphere/validation_results/<set>/` |
| How to evaluate your model | [docs/EVALUATION_PROTOCOL.md](AIDE-atmosphere/docs/EVALUATION_PROTOCOL.md) §3 (inputs) and §9 |
| How the code fits together | [docs/ARCHITECTURE.md](AIDE-atmosphere/docs/ARCHITECTURE.md) |
| Why a choice was made | [docs/DECISIONS.md](AIDE-atmosphere/docs/DECISIONS.md) |
| What is not done yet | [docs/DEFERRED.md](AIDE-atmosphere/docs/DEFERRED.md) |
| What is scored, when, in which PR | [docs/PROGRESS.md](AIDE-atmosphere/docs/PROGRESS.md) |
| Agreement with jcm-strat | [validation_results/port_check.md](AIDE-atmosphere/validation_results/port_check.md) |
| Settled data conventions and working rules | [CLAUDE.md](CLAUDE.md) |

## Evaluate a model

1. Write daily zonal means for 2006-01-01 … 2019-12-31 on pressure levels to netCDF:
   `Uzm`, `Vzm`, `THzm`, `VTHzm` (eddy flux, formed before time averaging), `AOA150`
   (entry age in years, reset below 150 hPa); `T` optional; dims `(time, plev, lat)`.
2. Add `AIDE-atmosphere/sources/<set>.json`:

   ```json
   {"description": "MYMODEL v1, 2006-2019", "reader": "netcdf", "window": [2006, 2019],
    "runs": {"MYMODEL": "/path/to/mymodel_*.nc"}}
   ```

3. From `AIDE-atmosphere/scripts/`:

   ```bash
   PY=../../.AIDE-eval_env/bin/python
   bash copy_era5_refs.sh                    # once: ERA5 copies, checked against the manifest
   $PY 01_fields.py --inspect <set>          # check the fields against the contract
   $PY evaluate.py <set>                     # 01 -> 02 -> 03, summaries rewritten
   ```

4. Read `validation_results/<set>/scorecard_standard.md` and the ranking in
   `validation_results/SCORING_STANDARD.md`.

## Running everything

Set the environment up once, from the repo root:

```bash
python3 -m venv .AIDE-eval_env
.AIDE-eval_env/bin/pip install -r requirements.txt
```

Then, from `AIDE-atmosphere/scripts/` (about 80 s on CPU; wrap in
`/data/JCM_stripped/jcm-strat-phase29/scripts/memcap.sh --` on shared nodes):

```bash
PY=../../.AIDE-eval_env/bin/python
bash copy_era5_refs.sh
$PY evaluate.py                # every set in sources/
$PY checks/port_check.py       # evidence: agreement with jcm-strat
```

- Model files and the CLaMS, WACCM6 and CESM references are read in place under `/data`;
  the ERA5 files are copied into the gitignored `AIDE-atmosphere/data/era5/`.
- `requirements.txt` pins numpy, scipy, xarray, matplotlib, cftime and netCDF4 on Python
  3.10.12. Two runs from an empty `output/` reproduce all 23 committed artefacts bit for bit
  under them.

## Repository layout

```
README.md              this file
CLAUDE.md              working rules, settled data conventions, scientific rules
LICENSE
requirements.txt       pinned dependencies
.gitignore
.AIDE-eval_env/        the virtual environment — gitignored, rebuild from the above

AIDE-atmosphere/
  sources/             one JSON per result set
  scripts/             strat_eval_common.py; 01_fields, 02_scorecard, 03_figures, evaluate.py;
                       checks/port_check.py; copy_era5_refs.sh + era5_refs.sha256
  docs/                EVALUATION_PROTOCOL, ARCHITECTURE, DECISIONS, DEFERRED, PROGRESS
  validation_results/  SCORING_STANDARD.md, SCORING_EXTENSIVE.md, port_check.md, <set>/
  output/              per-run JSON          — generated, gitignored
  data/                ERA5 reference copies — copy_era5_refs.sh, gitignored
```

`validation_results/` is generated and committed; `output/` is generated and not. Never
hand-edit either — change the script and re-run.

## History

Until 2026-11-12 this repository held a two-tier CESM-anchored protocol (5-year screening,
35-year validation) and its results for CESM, four reanalyses and jcm-strat phase 10. It is
archived as git tag `tier-protocol-final`.
