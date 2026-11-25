# Decisions — AIDE-SAI validation

**Settled** = decided by Susanne. **Unconfirmed** = taken during implementation; listed in
the PR for review and open until confirmed or rejected.

## Settled

### 2026-11-24
- Separate section in a sibling `AIDE-SAI/` directory on `main` (not AIDE-atmosphere).
- Reference: CCMI-2022 SAI (CESM2-WACCM6 aerosol prescribed into senD2-sai). Compare with the
  SAI fields as they are; no SAI-minus-control on the CCMI side (senD2-fix is not on /data).
- Mass mixing ratios from REMAP-CCMI-2022-sai, downloaded into gitignored `AIDE-SAI/data/`
  with a checksum manifest.
- Report only: deviations and figures, no pass/fail.
- New AIDE-SAI runs designed for the comparison are allowed.
- Commits and PRs keep the Claude attribution lines (overrides `CLAUDE.md` rule 4 for this work).

### 2026-11-24 (scope narrowed)
- CESM only; JCM later.
- 2030s rate first; 3-year runs.
- Window 2030–2034 of Fig. 2b; yearly ramp (annual steps), not a constant rate.
- Rates read off Tilmes et al. (2025) Fig. 2b (no rates from the authors in time).
- Injection at 75 hPa (Fig. 2a, 30°N/S at 18.1 km); current `cesm_numerical` preset.
- Ramp implemented in the AIDE-SAI repo as a schedule file (`INJ_SO2_SCHED_FILE`, PR #9, merged).

### 2026-11-25 (review of the Fig. 2b reading)
- Fig. 2b is read at the integer-year vertex (one vertex per year), not averaged over [Y, Y+1);
  the first reading was 0.3–0.5 Tg/yr (~10 %) high.
- Corrected schedule in AIDE-SAI (PR #12) and a new run `ccmi2030_inj_tick` with it; the
  running `ccmi2030_inj` (first reading) is kept as a +10 % rate-sensitivity case.
- The moved WACCM6 entry-age path is fixed in this PR.
- The resume seam in `ccmi2030_ctrl` (step 2165) is recorded in PROGRESS only; the control is kept.

## Unconfirmed — PR 1 (this PR)

1. **Red/green curves of Fig. 2b are per site.** Rate at each of 30°N and 30°S, not the pair
   sum: while green is zero, twice red equals black within 0.04 Tg/yr (2030: 1.92 × 2 vs 3.81).
   The schedule uses the black total with `INJ_MIRROR=1` (equal halves), consistent with this.
2. **Equal hemispheric split.** The figure has one red line per member for 30°N and 30°S
   together, so a N/S asymmetry cannot be read from it; equal halves are assumed.
3. **Year value = member mean at the tick**, linearly interpolated between the two pixel
   columns bracketing x(Y); lo/hi = member spread there (zero where the members overlap).
4. **Figure not committed**; fetched from the journal with a sha256 manifest (CC BY 4.0 would
   allow committing it).
5. **Run year 0 treated as spin-up** in the comparison (AIDE-SAI starts without the pre-2030
   aerosol that CCMI 2030 carries: ≤ 0.4 Tg/yr about 2020.5–2022.7, rising from 2024).
6. **Digitiser reads the PNG with `matplotlib.pyplot.imread`** (no new pin; Pillow comes in
   through matplotlib).
