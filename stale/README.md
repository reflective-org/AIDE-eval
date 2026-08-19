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
