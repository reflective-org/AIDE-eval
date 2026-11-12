# Deferred and open

What is knowingly not done, and the open issues that can change a number.

## Open issues

| Issue | Affects | State |
|---|---|---|
| **AIDE-eval #5** — w* estimator | B1, B2, D4-w, D7 | Interim: jcm-strat Ψ* route, log-p interpolated to the target pressures (protocol §7). A different estimator changes the WACCM6 references and every model value of these metrics. |
| **AIDE-eval #7** — group C daily-mean input | C1, C2 | C uses monthly means of daily means; against jcm-strat's 6-hourly monthly means the difference is ≤ 0.022 K (C2), ≤ 0.37 % (jets), one grid row in one jet latitude (protocol §8). Open: daily input, or require 6-hourly. |

## Deferred

- **ERA5 1990–2022 context for D** (jcm-strat `variability_score.py --long`): the scatter of
  D statistics over 5- or 14-year windows in ERA5, to size the sampling error against the
  25 % tolerance. Its yearly files are not copied.
- **ERA5 TEM w* context** (jcm-strat `era5_tem_monthly`, 1996–2001): shown by jcm-strat beside
  B, never scored. Not ported.
- **CDS download script** (jcm-strat `fetch_era5_daily_zm.py`): the ERA5 files are copies of
  jcm-strat's; regenerating them from the CDS is not part of this repository.
- **Windows outside 2006–2019**: need ERA5 files of those years; CLaMS covers 2004–2023,
  the WACCM6 entry age 1970–2019, WACCM6 h6 1970–2014.
- **B period matching**: the WACCM6 h6 reference is 1996–2014 for every window; no WACCM6 TEM
  output after 2014 is available here.
- **Phase 31's absolute-bound group scores** and its imposed / tuned / emergent classes of
  each metric (`make_phase31_report.py`): not ported; standard takes only its metric set.
- **More years for the seven phase-29 runs**: only 2006–2010 exists; scoring them on the
  standard window needs 2011–2019 runs.
- **An AL7h run whose T_eq correction is fitted outside the window**, for an out-of-sample
  standard score.
- **First AIDE-WACCM emulator rollout**: none scored yet; it needs `AOA150` and daily `VTHzm`
  formed before time averaging (protocol §3).
