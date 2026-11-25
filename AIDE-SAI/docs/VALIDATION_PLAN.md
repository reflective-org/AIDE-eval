# AIDE-SAI validation plan

Status (2026-11-25): **step 1 approved and running** (CESM, 2030s rate, 3-year runs; corrected-rate run added 2026-11-25). No
comparison has been made yet. JCM and the mid-century rate are deferred ([DEFERRED.md](DEFERRED.md)).

Goal: compare the AIDE-SAI aerosol from an SO2 injection with the CCMI-2022 SAI aerosol
distribution that CESM2-WACCM6 produced for the same injection. Report only (no pass/fail).

## 1. Reference

| Item | Value | Source |
|---|---|---|
| Source run | "SSP2-4.5 SAI": fully coupled CESM2-WACCM6, interactive MAM4 aerosol, 2 members | Tilmes et al. 2025, Sect. 2.2.2 |
| Injection 2030–2034 | SO2 at 30°N and 30°S, 180°, ~75 hPa (18.1 km); 15°N/15°S ≤ 0.2 Tg/yr per site until 2034 | Tilmes et al. 2025, Sect. 2.2.2 and Fig. 2b |
| Rate | feedback-controlled; total 3.81, 4.72, 5.60, 6.34, 6.70 Tg SO2/yr in 2030–2034, read at the integer-year vertex (± 0.3 read error; the two members overlap in 2030–2032) | [`validation_results/injection_rates/tilmes2025_fig2b.csv`](../validation_results/injection_rates/tilmes2025_fig2b.csv) |
| Prescribed product | REMAP-CCMI-2022-sai: monthly zonal-mean sulfate per mode, SAD, effective radius, AOD | doi:10.3929/ethz-b-000714654 |
| On /data | `CCMI_SAI_SST2025_001` (senD2-sai, WACCM6 reading REMAP, fixed 2020–2029 SSTs), monthly h0: `SAD_SULFC`, `RAD_SULFC`, `AEROD_v`, `REFF_AERO`, `SAD_AERO`, `T`, `TROP_P` | `/data/cesm2.1.5_output/CCMI_SAI_SST2025_001/month_1` |
| Not on /data | prescribed mass mixing ratios (`VOLC_MMR*` not written); senD2-fix baseline | — |

References:
- Tilmes, S., et al.: Stratospheric Aerosol Intervention experiment for the Chemistry–Climate
  Model Initiative, ACP, 25, 6001–6023, doi:10.5194/acp-25-6001-2025, 2025.
- Jörimann, A., et al.: Middle atmosphere chemical and dynamical effects in the CCMI-2022
  stratospheric aerosol injection scenario, ACP, 26, 11207–11234, doi:10.5194/acp-26-11207-2026, 2026.

## 2. Step 1 setup (CESM, 2030s)

| Item | Value |
|---|---|
| AIDE-SAI | `ccmi2030_inj`, `ccmi2030_ctrl`: `main` 891f782 (PR #9: `INJ_SO2_SCHED_FILE`); control resumed under 68e6a66 (help text only). `ccmi2030_inj_tick`: c3103a6 (PR #12, on `main` 4ac8d61; input file only) |
| Preset | `cesm_numerical` (f09, 24 levels, 1.2–143 hPa), current defaults |
| Injection | annual steps; `--inj-lat 30 --inj-mirror 1` (half at each site), lon 180, point source, 75 hPa (model level 73.8 hPa) |
| Runs (26 280 h = 3 × 365 d, all-level zonal-mean frames) | `ccmi2030_inj_tick`: corrected rates 3.81, 4.72, 5.60 Tg SO2/yr (`ccmi2022_sai_30NS_2030-2034.csv`), **the comparison run**. `ccmi2030_inj`: first reading 4.27, 5.15, 6.00 (`*_yearmean.csv`, ~10 % high), rate-sensitivity case. `ccmi2030_ctrl`: zero injection |
| Compared with | CCMI 2030–2032, year by year |

## 3. Comparability limits (stated with every result)

- **Spin-up**: CCMI 2030 carries aerosol from earlier injection (small, ≤ 0.4 Tg/yr, about
  2020.5–2022.7; then rising from 2024 to 2.7 Tg/yr in 2029); AIDE-SAI starts from the
  background state in run year 0. Proposed: run year 0 is spin-up and years 1–2 are compared,
  with the caveat that the CCMI burden includes the pre-2030 injection (DECISIONS U5).
- **Meteorology**: AIDE-SAI reads CESM meteorology (prescribed, one-way); CCMI 2030s are under
  SSP2-4.5 with SAI feedback on dynamics. Same model family, different years and climate.
- **Domain**: AIDE-SAI covers 1.2–143 hPa (`cesm_numerical`) with open faces and no wet removal; compare in that
  band and above the tropopause only.
- **AOD**: AIDE-SAI AOD is band-only, sulfate-only and uncalibrated (AIDE-SAI `MANIFEST.md`).
- **Reference**: one model (WACCM6); REMAP zeroes aerosol below each model's tropopause.
- **Rate**: digitised from a figure (± 0.3 Tg/yr); hemispheric split taken as equal (the figure has one
  red line per member for 30°N/30°S together).

## 4. Tasks

### Task 1 — Scaffold, docs, injection rates (PR 1, this PR)
1. `AIDE-SAI/{scripts,docs,validation_results}`; gitignored `AIDE-SAI/data/` and `AIDE-SAI/output/` → verify: AIDE-atmosphere chain unchanged.
2. `fetch_tilmes2025_fig2.sh` + sha256 manifest; `digitize_tilmes2025_fig2b.py` → CSV and check figure → verify: 2030–2034 totals equal the AIDE-SAI schedule file (PR #12); twice the red curve agrees within 0.04; two runs bit-identical.
3. Root `README.md` section 2 and `CLAUDE.md` layout; docs (this plan, ARCHITECTURE, DECISIONS, DEFERRED, PROGRESS).
4. Fix the moved WACCM6 entry-age path in AIDE-atmosphere → verify: all 23 artefacts unchanged except date and path.

### Task 2 — CCMI reference (PR 2)
1. Reader for the on-disk SAI fields → monthly zonal means on the AIDE-SAI levels → verify: units, h0 stamp is the interval end, fill handling.
2. Download REMAP with sha256 manifest (`fetch_remap.sh`) → verify: manifest check.
3. Ingestion check: REMAP SAD vs WACCM `SAD_SULFC` → verify: agreement within interpolation error.
4. Reference figures for 2030–2032: lat–p sections (sulfate MMR, SAD, r_eff), AOD vs latitude, band burden vs time.

### Task 3 — AIDE-SAI runs (AIDE-SAI repo; running)
1. Yearly schedule (AIDE-SAI PR #9, merged) → verified: constant-rate and one-row runs bit-identical to `main`.
2. Corrected rates (AIDE-SAI PR #12) → verify: CI test.
3. `ccmi2030_inj_tick`, `ccmi2030_inj`, `ccmi2030_ctrl` → verify: budget closes; rate steps at h 8760 and h 17 520 (`ccmi2030_inj`: h 8760 seen).

### Task 4 — CESM evaluation (PR 3, report only)
`ccmi2030_inj_tick` per year (1, 2; year 0 shown, marked spin-up), band 1.2–143 hPa, above the tropopause;
`ccmi2030_inj` shown as the +10 % rate case:
- Sulfate burden (Tg S), injection-normalised; burden vs time.
- Zonal-mean sulfate MMR, SAD, effective radius: lat–p (model, reference, difference; pattern correlation, relative deviation).
- Profiles at the equator, ±30°; altitude and latitude of the sulfate peak.
- Burden and r_eff response between the two injection runs (rate sensitivity).
- AOD vs latitude.
- `EVALUATION_PROTOCOL.md` written with this PR.

## 5. Open items
- REMAP: CC BY 4.0, 9 files, 14.9 GB, MD5 given; no file on the WACCM/CESM grid (remapped to ECHAM6 and MIROC32) and no variable description. Which file to use (Task 2).
- Whether to use SAI-minus-control on the AIDE-SAI side only, or AIDE-SAI absolute fields (Task 4; decided at review of Task 2's figures).
