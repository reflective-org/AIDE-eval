# Decisions

Every decision behind the evaluation, with who took it. **Settled** = decided by Susanne.
**Unconfirmed** = taken during implementation because the brief did not settle it; listed in
the PR for review and open until confirmed or rejected. Numbers 1–34 are PR #6's, 35 onward
those of PR #8 `strat-eval-standard` ([PROGRESS.md](PROGRESS.md)).

---

## Settled

### 2026-11-12 (brief 3)

- The transport scorecard **replaces** the two-tier CESM protocol and is the repository's
  standard evaluation of stratospheric dynamics.
- Tier material archived as annotated tag `tier-protocol-final` (PR #6 head 34d7d6f) and
  deleted from the tree: tier scripts 01, 01b, 01c, 01d, 02, 02b, 07, 16, 17, 18, 19, 20,
  21, `report_layout.py`, `stale/`, the tier results, `era5_requirements.txt`. The scorecard's
  share of `aide_val_common.py` merged with `transport_common.py` into
  `strat_eval_common.py`. The h6 conventions stay documented, citing the tag for the evidence.
- Scripts renumbered `01_fields.py`, `02_scorecard.py`, `03_figures.py`;
  `checks/port_check.py`; `copy_era5_refs.sh` and its manifest stay.
- One JSON per result set in `sources/<set>.json` (reader, window, runs, description); one
  command `evaluate.py <set> [...]` runs 01 → 02 → 03 and rewrites the summaries. No package.
- Results per set in `validation_results/<set>/`; `SCORING_STANDARD.md` and
  `SCORING_EXTENSIVE.md` at the top level.
- Scorings renamed: **standard** (= the former "reduced", jcm-strat phase 31's scored set,
  A2-WACCM for A2) is the headline score and ranking; **extensive** (34 metrics) the full
  listing.
- **Standard window 2006–2019**; reference years follow the window (A, C, D; B stays WACCM6
  h6 1996–2014). Other windows flagged "non-standard window (N of 14 standard years)".
  AL7h scored 2006–2019 from phase 29's `AL7h_0610.npz` and phase 31's `AL7h_1119.npz` as one
  record (join and grids checked first); its T_eq correction was fitted on
  2005-04..2008-01, inside the window. The other seven phase-29 runs on 2006–2010 only;
  QBO30 and QBO10dL76 excluded. The AL7h in-sample / out-of-sample comparison dropped.
- Port check keeps part 1 (phase 29, eight runs, 2006–2010, age references 2005–2009) and
  part 2 (phase 31 AL7h 2006–2010 and 2011–2019) with its own source definitions; part 3
  (age-reference switch) dropped.
- Docs: protocol rewritten as the scorecard protocol (TRANSPORT_SCORECARD.md folded in);
  ARCHITECTURE, DECISIONS, DEFERRED, PROGRESS added; README and CLAUDE.md project context
  rewritten, the working rules kept verbatim.

### 2026-11-11 (brief 2)

- AL7h scored out of sample 2011–2019 (phase 31 `AL7h_1119.npz`, read in place) beside
  2006–2010. *Superseded 2026-11-12 by one 2006–2019 record.*
- Reference years follow the scored window for every set; the phase-29 set moved to
  2006–2010 age references. The WACCM6 w* / Ψ* reference stays h6 1996–2014.
- Two scorings of every set: extensive (34 metrics, unchanged) and reduced (phase 31's scored
  set, A2-WACCM for A2, relative scoring). *Reduced renamed standard 2026-11-12.*

### 2026-11-10 (brief 1)

- A transport scorecard beside the two CESM tiers, tier thresholds unchanged. *Superseded
  2026-11-12: it replaces them.*
- Port the diagnostics, not the cached numbers; score the eight phase-29 runs over 2006–2010
  (QBO30 and QBO10dL76 excluded).
- w*: jcm-strat formula (`tem_streamfunction` + `residual_w`, no w̄), tropical values
  interpolated linearly in log p to the target pressures, \|lat\| > 88° masked — **interim,
  AIDE-eval issue #5 stays open**. Mass flux B2: range of Ψ* within 60°, interpolated to
  100 / 70 hPa. WACCM6 B references recomputed from h6 1996–2014 with the same method.
- Age: entry age only, clock reset below 150 hPa (`aoa150`).
- Field contract: daily `Uzm`, `Vzm`, `THzm`, `VTHzm`, `AOA150` [yr] on `(time, plev, lat)`,
  p in hPa, `Wzm` optional; a `jcm_npz` reader and a generic `netcdf` reader.
- Scoring = jcm-strat `relative_scores.py`: rel = model/ref − 1, err = \|rel\|/0.25; A2,
  C1 latitudes, C2 and D4 dates absolute; A4 against its midpoint 2.75; group scores and
  combined A–C.
- ERA5 files copied into the gitignored `AIDE-atmosphere/data/era5/` with a sha256 manifest;
  CLaMS, WACCM6 and model files read in place.

---

## Unconfirmed — PR #6 (1–34)

Status after this PR in italics.

1. Readers and sources in the common module, in the `reader(spec, years)` form of the former
   script 20, not in 20 itself. *Readers stay in `strat_eval_common`; sources moved to JSON (40).*
2. `netcdf` reader: glob paths; dims `time`, `plev`, `lat`; values beyond 1e20 masked; `T`
   from `THzm` with κ = R/c_p if absent; days labelled by their own date (noleap works);
   checked by a round trip of AL7h (script kept outside the repo).
3. A references at jcm-strat's 2005–2009. *Superseded (settled 2026-11-11).*
4. CESM constants a, g, κ instead of jcm-strat's; ≤ 0.008 % on any metric.
5. WACCM6 h6 interfaces from the model top to 652 hPa, sentinel masked; jcm-strat read all 71
   unmasked; no metric level affected.
6. Mass flux: Ψ* interpolated in log p first, then its range over \|lat\| ≤ 60°.
7. D4-w and D7 w* also interpolated to 70 hPa (one estimator for every w*).
8. A (55 hPa), C (50, 10 hPa), D (10, 50, 55, 100 hPa) at the nearest model level, as
   jcm-strat.
9. D3: jcm-strat's SSW detector, ported, not the tier implementation.
10. C from monthly means of the daily means (one-frame-per-month shift; issue #7).
11. Port-check criterion max(0.5 %, jcm-strat rounding), B/D4-w/D7 at the nearest level; a
    separate evidence script. *Now `checks/port_check.py`.*
12. Artefact names stamped `__<set>__<date>`; scorecard JSON committed beside the report.
    *Superseded by 45 (undated names in per-set folders).*
13. ERA5 copies: manifest `scripts/era5_refs.sha256`, verified on every copy; only the files a
    metric reads.
14. Script 22 interface `[NAME ...]`, `--inspect NAME`, references recomputed per call.
    *Now `01_fields.py [SET ...]`, `--inspect SET`.*
15. jcm runs' `T` taken from the npz for C2 and D6, not derived from `THzm`.
16. Figures: profiles as small multiples; heatmap colours absolute metrics as a fraction of
    their bound (bound = 25 %), 12.5 % steps, bold beyond tolerance.
17. Run descriptions from the jcm-strat phase 22–29 docs. *Kept as `notes` in the set JSON.*
18. A1 reference is the WACCM6 entry age (code), not "CLaMS 1.11" (phase-29 output.md).
19. No package added.
20. Windows on sources, sets in a `SETS` dict. *Superseded by the set JSON (brief 3).*
21. One references file per window. *Now one per set (42).*
22. Extensive shows exactly 34 metrics (A2-WACCM computed, hidden); reduced lists A2 as
    unscored. *Kept for standard.*
23. A2-WACCM in group A; WACCM6 hybrid `lev` taken as pressure; bias stored beside the RMSE,
    not scored.
24. Reduced D = D4-w, D5-N, D5-S. *Kept for standard.*
25. AL7h comparison placement and colours. *Superseded: comparison dropped (brief 3).*
26. Script interfaces and stamps `__<set>__<scoring>__<date>`. *Superseded by 43, 45.*
27. Unscored metrics faded in the heatmap, "(not scored)" in the label.
28. Port check part 1 recomputes A at 2005–2009. *Kept.*
29. Phase 31 D for AL7h 2006–2010 compared with phase 30's `variability.json`. *Kept.*
30. ERA5 copied from `jcm-strat/cache` (phase 31's cache paths resolve there).
31. Interpretation text of TRANSPORT_SCORECARD.md §8–9 written in that update. *Removed with
    that document; see 53.*
32. Summary names `SCORING_EXTENSIVE.md` / `SCORING_REDUCED.md`. *Now `SCORING_STANDARD.md`,
    `SCORING_EXTENSIVE.md` at the top level (brief 3).*
33. Summary content: group scores, fail/pass in every run, links. *Kept; ranking added (46).*
34. Re-stamping 2026-11-11 artefacts as 2026-11-12. *Moot (45).*

---

## Unconfirmed — PR #8 (35–)

35. **Tag.** `tier-protocol-final` on 34d7d6f, annotated "Last state with the two-tier CESM
    evaluation protocol, its scripts (01-21, report_layout.py, stale/) and results …
    Superseded on 2026-11-12 by the stratospheric transport scorecard"; pushed before any
    deletion.
36. **Also deleted:** the `era5_env/` line of `.gitignore` (that environment served only 19
    and 21). `requirements.txt` pins unchanged (scipy is still needed by xarray's `interp`
    in group C); its header comment updated to the scorecard's artefact count.
37. **`strat_eval_common` contents.** From `aide_val_common`: constants, `OUTDIR`, the
    1996–2014 h6 segment only (1970–1995 is read by nothing), `h6_path`, `load_h6` (ilev slice
    now a required argument), `interp_level`. Its docstring and the docs list four h6
    conventions — `Wzm` log-pressure added to the three the brief names, because the contract
    names `Wzm`. Scripts import it as `SE`. The jcm-strat paths used only by the port check
    moved into `checks/port_check.py`.
38. **Setup files** `copy_era5_refs.sh` and `era5_refs.sha256` stay in `scripts/`.
39. **Sources location** `AIDE-atmosphere/sources/`.
40. **Sources schema.** Keys of the brief plus optional `notes` (one line per run, from the
    former run descriptions). `long_name`, `calendar`, `independent`, `time_resolution` dropped
    (read by nothing). `window` as `[Y0, Y1]`.
41. **Record guards.** `load_run` refuses a record that does not span 1 January Y0 to 31
    December Y1; `jcm_npz` refuses files with different grids or a join that is not one day.
42. **Output layout** `output/<set>/references.json` and `output/<set>/<run>.json`.
43. **`evaluate.py`** in `scripts/`; no argument = every set; each stage a subprocess of the
    same interpreter; the ERA5 copy is a separate setup step. 01 and 03 take `[SET ...]`;
    03's `--date` removed.
44. **Commit order** 2.5, 2.4, 2.6, 2.3: the comparison was removed and the layout fixed
    before `evaluate.py` was written against the final interfaces.
45. **Undated result names** (`<set>/scorecard_<scoring>.md` etc.): links stay fixed; the
    production date stays inside each report and JSON, so a rerun on another day changes
    that line only.
46. **Summaries** rebuilt from the JSON on disk for every set registered in `sources/` (sets
    without results listed); `SCORING_STANDARD.md` carries a ranking of every run of every
    set by standard combined score, window and flag shown; the extensive summary points to it.
47. **Flag placement**: a note under each report title, the window row, summary set
    headings, the ranking's window column, figure titles and the score-table footer.
48. **Set names** `jcm-strat-AL7h-2006-2019` (run `AL7h`) and `jcm-strat-phase29-2006-2010`
    (kept for the seven runs).
49. **Port check self-contained**: computes its diagnostics through 01's functions instead of
    reading `output/`; `port_check.md` at the top of `validation_results/`.
50. **AL7h join check beyond the brief**: the age clock and the 10 / 50 hPa means across
    2010-12-31 → 2011-01-01 compared with the daily steps (clock step 0.05 d against a median
    0.5 d); log kept outside the repo (`aide_transport_job/checks/al7h_join_check.log`).
51. **Profiles** drawn from the standard scorecard JSON (identical content in both).
52. **Names kept**: reports and figures still titled "Transport scorecard".
53. **Interpretation** in protocol §10.2 (ranking gaps against the window effect) written in
    this PR; PR #6's §8–9 interpretation not carried over.
54. **Tag citations**: each document mentions `tier-protocol-final` at most once, where it
    refers to the tier protocol (DECISIONS records the decisions about it).
55. **AL7h per-window numbers** (2006–2010 alone, 2011–2019 alone) are no longer result sets;
    they are cited from PR #6 in the protocol, and the port check still checks both files.
56. **Branch in the main checkout.** git refuses one branch in two worktrees, so the
    worktree's HEAD was detached at the same commit before the main checkout switched to
    `strat-eval-standard`.
57. **Worktree logs kept.** The worktree's gitignored `logs/` (18 logs of PR #6's port and
    reader checks, no other copy) moved to `AIDE-atmosphere/logs/strat-transport-scoring/`
    in the main checkout before the worktree was removed; its `output/` was discarded
    (regenerated).
58. **Repo-relative paths in result JSON** (`strat_eval_common.portable`, commit 910765a):
    the references block carried the checkout's absolute ERA5 path, so the committed JSON
    differed between checkouts; found while reproducing from the main checkout.
59. **`AIDE-atmosphere/era5_env/` in the main checkout** (120 MB, the tier-only ERA5 download
    environment) is left in place and now shows as untracked, since its ignore line is
    removed (36). Delete it, or restore the ignore line.
60. **ERA5 move verified** with the manifest check of `copy_era5_refs.sh` run directly
    (`sha256sum -c`), not by re-running the copy.
