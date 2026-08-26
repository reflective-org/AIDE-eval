# stale

Superseded material, kept for the record. Nothing here is maintained, and no number
in it should be quoted as current.

| File | What it was | Why it is stale |
|---|---|---|
| `14_evaluation_tiers.py` | Tier thresholds from the split anchors — tier 1 on 1996–2014, tier 2 on 1980–2014 | Superseded by `scripts/16_anchors_45yr.py`, which sets both tiers on the whole 1970–2014 record |
| `15_screen_out_of_sample.py` | Built the screening figure: both tiers scored on CESM 1970–1995 | The out-of-sample test it performs is no longer possible — every CESM year is inside the operational anchor |
| `AIDE_WACCM_screening_1970-1995.pdf` | That figure, 3 pp | Its anchors are the split ones |
| `AIDE_WACCM_validation_targets.pdf` | The original 23 pp technical report, first artefact the project produced | Documents the 10/20-year target suite that the tier protocol replaced |

## What these two still establish

The split anchors existed to satisfy one rule: a threshold CESM cannot meet on a
different sample of its own output is not a threshold. `15` is the evidence that the
±3σ gate and the 0.5σ tolerance behave sensibly on CESM output that did not set them
— 28 of 30 tier-1 block verdicts pass, and 4 of 6 tier-2 mean tests, the failures
being the upwelling pair and the forced BDC trend behind them.

The operational anchor cannot re-earn that evidence, because no CESM years are left
outside it. So the conclusion still stands on these two scripts even though their
numbers are superseded, and `docs/EVALUATION_PROTOCOL.md` §6 says so.

## Running them

Both import from `scripts/`, so they need it on the path:

```bash
PYTHONPATH=../scripts $PY 14_evaluation_tiers.py
PYTHONPATH=../scripts $PY 15_screen_out_of_sample.py
```

`14` reads `output/07_period_split.json` and `output/02b_trends.json`; `15` reads `14`'s
output. Both still run and still reproduce their JSON, and their output paths were
repointed here so that running them cannot write into `docs/` or overwrite anything live.

---

## `22_era5_monthly_series.py` and `23_era5_monthly_tier1.py`

Retired 2026-08-26, superseded by the general scoring path.

They were the first route from a reanalysis to a scorecard: `22` reduced an ERA5
monthly tape to the four diagnostics a monthly archive supports, and `23` scored
them against the anchor and drew the figures. Both worked. The problem was
duplication — `23`'s `score()` reproduced `17_validate.py`'s tier-1 block, its
`score_seasonal()` reproduced 17's shape block, and it overrode 18's footer. Two
scoring paths meant two places to keep a threshold change in step.

What replaced them:

- `20_series_from_zonal_mean.py` takes **any** zonal-mean source — reanalysis,
  another GCM, an emulator rollout — and writes the standard series JSON. Adding a
  source is a few lines of spec, not a new script.
- `17_validate.py --source NAME` scores it. Diagnostics the source cannot supply
  are reported as **not evaluable** rather than omitted, and results go to
  `validation_results/<NAME>/`.

Two ideas from these scripts survive in the replacement and are worth keeping
attached to their origin. The first is *why monthly means are exact* for the four
diagnostics that survive: each is a chain of linear operations, and linear
operators commute with time averaging, so a day-weighted mean of monthly means
equals the mean of the daily series — `20.expand_to_days` recovers that weighting.
The second is that **independence is a property of the source, not of its period**:
`23` had to override 18's footer because a reanalysis whose years lie inside the
anchor is still an independent product. `20` now declares it and `17` carries it.

Neither script runs any more: the tape `22` read is not rebuilt, and `23` imports
`18` under its old signatures. They are kept for the reasoning, not to execute.
