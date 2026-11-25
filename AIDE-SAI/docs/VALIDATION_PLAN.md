# AIDE-SAI validation plan

Status (2026-11-25): **step 1 approved and running** (CESM, 2030s rate, 3-year runs). No
comparison has been made yet. JCM and the mid-century rate are deferred ([DEFERRED.md](DEFERRED.md)).

Goal: compare the AIDE-SAI aerosol from an SO2 injection with the CCMI-2022 SAI aerosol
distribution that CESM2-WACCM6 produced for the same injection. Report only (no pass/fail).

## 1. Reference

| Item | Value | Source |
|---|---|---|
| Source run | "SSP2-4.5 SAI": fully coupled CESM2-WACCM6, interactive MAM4 aerosol, 2 members | Tilmes et al. 2025, Sect. 2.2.2 |
| Injection 2030–2034 | SO2 at 30°N and 30°S, 180°, ~75 hPa (18.1 km, Fig. 2a); 15°N/15°S start ~2033 and are < 0.2 Tg/yr per site until 2034 | Tilmes et al. 2025, Fig. 2 |
| Rate | feedback-controlled; total 4.27, 5.15, 6.00, 6.39, 6.87 Tg SO2/yr in 2030–2034 (digitised, ± 0.3 read error, member envelope ± 0.4) | [`validation_results/injection_rates/tilmes2025_fig2b.csv`](../validation_results/injection_rates/tilmes2025_fig2b.csv) |
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
| AIDE-SAI | `main` 891f782 (PR #9: `INJ_SO2_SCHED_FILE`); control resumed under 68e6a66 (help text only) |
| Preset | `cesm_numerical` (f09, 24 levels, 1.2–143 hPa), current defaults |
| Injection | `inputs/injection/ccmi2022_sai_30NS_2030-2034.csv`: 4.27, 5.15, 6.00 Tg SO2/yr in run years 0, 1, 2 (annual steps); `--inj-lat 30 --inj-mirror 1` (half at each site), lon 180, point source, 75 hPa (model level 73.8 hPa) |
| Runs | `runs/ccmi2030_inj` and zero-injection `runs/ccmi2030_ctrl`, 26 280 h (3 × 365 d), all-level zonal-mean frames |
| Compared with | CCMI 2030–2032, year by year |

## 3. Comparability limits (stated with every result)

- **Spin-up**: CCMI 2030 carries aerosol from injection since 2025; AIDE-SAI starts from the
  background state in run year 0. Run year 0 is spin-up; years 1–2 are compared with the
  caveat that the CCMI burden includes 2025–2029 injection.
- **Meteorology**: AIDE-SAI reads CESM meteorology (prescribed, one-way); CCMI 2030s are under
  SSP2-4.5 with SAI feedback on dynamics. Same model family, different years and climate.
- **Domain**: AIDE-SAI covers 1–150 hPa with open faces and no wet removal; compare in that
  band and above the tropopause only.
- **AOD**: AIDE-SAI AOD is band-only, sulfate-only and uncalibrated (AIDE-SAI `MANIFEST.md`).
- **Reference**: one model (WACCM6); REMAP zeroes aerosol below each model's tropopause.
- **Rate**: digitised from a figure (± 0.3 Tg/yr); hemispheric split taken as equal.

## 4. Tasks

### Task 1 — Scaffold, docs, injection rates (PR 1, this PR)
1. `AIDE-SAI/{scripts,docs,validation_results}`; gitignored `AIDE-SAI/data/` and `AIDE-SAI/output/` → verify: AIDE-atmosphere chain unchanged.
2. `fetch_tilmes2025_fig2.sh` + sha256 manifest; `digitize_tilmes2025_fig2b.py` → CSV and check figure → verify: reproduces the AIDE-SAI schedule file, two runs bit-identical.
3. Root `README.md` section 2 and `CLAUDE.md` layout; docs (this plan, ARCHITECTURE, DECISIONS, DEFERRED, PROGRESS).

### Task 2 — CCMI reference (PR 2)
1. Reader for the on-disk SAI fields → monthly zonal means on the AIDE-SAI levels → verify: units, h0 stamp is the interval end, fill handling.
2. Download REMAP with sha256 manifest (`fetch_remap.sh`) → verify: manifest check.
3. Ingestion check: REMAP SAD vs WACCM `SAD_SULFC` → verify: agreement within interpolation error.
4. Reference figures for 2030–2032: lat–p sections (sulfate MMR, SAD, r_eff), AOD vs latitude, band burden vs time.

### Task 3 — AIDE-SAI runs (AIDE-SAI repo; running)
1. Yearly schedule (AIDE-SAI PR #9, merged) → verified: constant-rate and one-row runs bit-identical to `main`.
2. `ccmi2030_inj`, `ccmi2030_ctrl` → verify: budget closes; rate steps at h 8760 (seen) and h 17 520.

### Task 4 — CESM evaluation (PR 3, report only)
Per year (1, 2; year 0 shown, marked spin-up), band 1–143 hPa, above the tropopause:
- Sulfate burden (Tg S), injection-normalised; burden vs time.
- Zonal-mean sulfate MMR, SAD, effective radius: lat–p (model, reference, difference; pattern correlation, relative deviation).
- Profiles at the equator, ±30°; altitude and latitude of the sulfate peak.
- AOD vs latitude.
- `EVALUATION_PROTOCOL.md` written with this PR.

## 5. Open items
- REMAP file size and licence (Task 2).
- Whether to use SAI-minus-control on the AIDE-SAI side only, or AIDE-SAI absolute fields (Task 4; decided at review of Task 2's figures).
