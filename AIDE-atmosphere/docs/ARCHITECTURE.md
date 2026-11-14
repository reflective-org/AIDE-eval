# Architecture

How the stratospheric-dynamics evaluation is built: pipeline, modules, data flow, files.
The method is in [EVALUATION_PROTOCOL.md](EVALUATION_PROTOCOL.md).

## Pipeline

```
sources/<set>.json ──► 01_fields.py ──► output/<set>/<run>.json, references.json
                                              │
                                              ▼
                       02_scorecard.py ──► validation_results/<set>/scorecard_<scoring>.{json,md}
                                              │        └─► validation_results/SCORING_STANDARD.md,
                                              │            SCORING_EXTENSIVE.md (every registered set)
                                              ▼
                       03_figures.py  ──► validation_results/<set>/*.png

evaluate.py <set> [...]  =  01 → 02 → 03 for the sets
checks/port_check.py     =  evidence, outside the chain ──► validation_results/port_check.md
copy_era5_refs.sh        =  setup: ERA5 copies into data/era5/, checked against era5_refs.sha256
```

- **01 measures, 02 scores, 03 draws.** 02 and 03 read only JSON and take seconds; 01 reads
  the model and the references.
- Each stage runs for the named sets (default: every file in `sources/`).
- `output/` is gitignored and rebuilt by 01; `validation_results/` is generated and
  committed. Neither is edited by hand.

## Modules (`AIDE-atmosphere/scripts/`)

| File | Role |
|---|---|
| `strat_eval_common.py` | Everything the scripts share: CESM constants; `REFS` (every reference path); the h6 loader (`load_h6`, `interp_level`); the field contract (`CONTRACT`, `OPTIONAL`, `UNITS`); result-set registry (`set_names`, `load_set`, schema check); readers (`READERS`: `jcm_npz`, `netcdf`) and `load_run`; scoring (`TOL`, `ABS_BOUND`, `score`, `status`, `SCORINGS`, `STANDARD_WINDOW`, `window_flag`); grid helpers ported from jcm-strat `sweep_score.py`; age references (`clams_age`, `waccm_age`); residual circulation (`tem_streamfunction`, `residual_w`, `tropical_wstar`, `upward_mass_flux`) |
| `01_fields.py` | Per run: groups A–D (`group_A` … `group_D`) against the window's references (`references(window)`: A and D per window, B once per call, C per model grid). `--inspect SET` prints the contract fields |
| `02_scorecard.py` | Per set and scoring: `score_all` (rows, group scores), `markdown` (report); `write_summaries` rebuilds both summaries from the scorecard JSON of every registered set, with the ranking in the standard one |
| `03_figures.py` | Per set: score table and heatmap per scoring; w* and age profiles from the standard JSON |
| `evaluate.py` | Runs 01, 02, 03 as processes of the same interpreter; stops on the first failure |
| `checks/port_check.py` | Recomputes the diagnostics of the jcm-strat phase-29 runs and AL7h (both files) through 01's functions, compares with jcm-strat's JSON; runs and paths defined in the script (`P29_RUNS`, `P31_RUNS`) |
| `copy_era5_refs.sh`, `era5_refs.sha256` | Copy the 28 ERA5 files (2006–2019) a metric reads from jcm-strat's cache; verify the manifest |

- Imports: each script puts `scripts/` on `sys.path` and imports `strat_eval_common as SE`;
  `checks/port_check.py` adds its parent directory and loads 01 by file name.
- No package, no new dependency: `requirements.txt` (numpy, scipy, xarray, matplotlib,
  cftime, netCDF4) under Python 3.10.12.

## Data flow

| Input | Read by | Where |
|---|---|---|
| Model daily zonal means | `READERS[reader](paths, years)` via `load_run` | paths in `sources/<set>.json` (read in place) |
| CLaMS v3.1 zonal-mean age | `clams_member`, `clams_age`, `01.clams_monthly_age` | `/data/CLaMS/CLaMS_v3/clams_v3.1_era5_zm_lat.zip` (read from the zip in memory) |
| WACCM6 REF-D1 entry age | `waccm_age` | `/data/CESM2_REFD1_AOA/…_1970-2019_….nc` |
| CESM2.1.5-WACCM6 h6 `Vzm`, `VTHzm`, `THzm` | `load_h6("1996-2014", …)` | `/data/cesm2.1.5_output/histSST/…` |
| ERA5 monthly zonal means | `01.era5_monthly` | `AIDE-atmosphere/data/era5/era5_zm_monthly_<y>.nc` |
| ERA5 CDS daily zonal means | `01.Daily.era5` | `AIDE-atmosphere/data/era5/era5_<y>.npz` |

- Reference years: the set's `window` for A, C, D; 1996–2014 for B.
- In memory a run is a dict of `(day, plev, lat)` arrays with `plev`, `lat` ascending
  (`orient`), `days` as `datetime64[D]`, `AOA150` in years.

## Files

```
AIDE-atmosphere/
  sources/<set>.json                 result-set registry (committed)
  scripts/                           the modules above
  data/era5/                         ERA5 copies (gitignored; copy_era5_refs.sh)
  output/<set>/references.json       references of the set's window (gitignored)
  output/<set>/<run>.json            metrics (model, ref, extras) and profiles per run (gitignored)
  validation_results/
    SCORING_STANDARD.md              every set, standard scoring, ranking
    SCORING_EXTENSIVE.md             every set, extensive scoring
    port_check.md                    agreement with jcm-strat phases 29 and 31
    <set>/scorecard_<scoring>.json   rows, scores, profiles, references, window flag, date
    <set>/scorecard_<scoring>.md     report
    <set>/score_table_<scoring>.png, heatmap_<scoring>.png, wstar_profile.png, age_profile.png
  docs/                              EVALUATION_PROTOCOL, ARCHITECTURE, DECISIONS, DEFERRED, PROGRESS
```

### `sources/<set>.json`

| Key | Type | Meaning |
|---|---|---|
| `description` | string | what the set is |
| `reader` | `"jcm_npz"` \| `"netcdf"` | key of `READERS` |
| `window` | `[Y0, Y1]` | scored years, inclusive; the record must span 1 Jan Y0 – 31 Dec Y1 |
| `runs` | `{name: path \| [path, …]}` | a list is one continuous record; `netcdf` paths may be globs |
| `notes` | `{name: string}` | optional, one line per run |

### Per-run JSON (`output/<set>/<run>.json`)

`name`, `set`, `note`, `reader`, `paths`, `window`, `days` (first, last), `n_days`,
`n_levels`, `n_lats`, `metrics` (one row per metric: `group`, `id`, `label`, `model`, `ref`,
and for B, D4-w, D7 the nearest-level `model_nearest`, `ref_nearest`), `profiles` (tropical
age and w* profiles, w* annual cycle).

## Adding things

- **A model**: one `sources/<set>.json` (protocol §9).
- **A reader**: a function `reader(paths, years)` returning the contract dict, registered in
  `READERS`.
- **A metric**: a row in the relevant `group_*` of 01 with its reference; its scoring in
  `ABS_BOUND` if absolute; its place in `SCORINGS`. Changing the metric set changes what a
  score means — record it in DECISIONS.md.
