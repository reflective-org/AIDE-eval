# Architecture — AIDE-SAI validation

Sibling of `AIDE-atmosphere/`; shares the repo-root environment (`.AIDE-eval_env`,
`requirements.txt`) and nothing else. Paths below are relative to `AIDE-SAI/`.

## Layout

| Path | Content | Git |
|---|---|---|
| `scripts/fetch_tilmes2025_fig2.sh`, `tilmes2025_fig2.sha256` | download Tilmes et al. (2025) Fig. 2 and check it | committed |
| `scripts/digitize_tilmes2025_fig2b.py` | yearly injection rates from Fig. 2b | committed |
| `data/` | downloaded references (figure; REMAP from Task 2) | gitignored |
| `output/` | intermediate JSON/netCDF (from Task 2) | gitignored |
| `validation_results/` | generated reports, tables, figures | generated, committed |
| `docs/` | plan, architecture, decisions, deferred, progress | committed |

## Data flow (current)

```
acp-25-6001-2025-f02-web.png ──fetch (sha256)──▶ data/tilmes2025/
        └─ digitize_tilmes2025_fig2b.py ──▶ validation_results/injection_rates/tilmes2025_fig2b.{csv,png}
                                              │ copied by hand (2030–2034 totals)
                                              ▼
            AIDE-SAI repo: inputs/injection/ccmi2022_sai_30NS_2030-2034.csv ──▶ runs/ccmi2030_{inj,ctrl}
```

Planned (Tasks 2 and 4): CCMI reader and REMAP → `output/` → comparison with the AIDE-SAI
run output → `validation_results/<set>/`.

## Reproduce

From the repo root:

```bash
bash AIDE-SAI/scripts/fetch_tilmes2025_fig2.sh
.AIDE-eval_env/bin/python AIDE-SAI/scripts/digitize_tilmes2025_fig2b.py
```

Both outputs are bit-identical across runs (checked 2026-11-25).
