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
- Commits and PRs keep the Claude attribution lines.

### 2026-11-24 (scope narrowed)
- CESM only; JCM later.
- 2030s rate first; 3-year runs.
- Window 2030–2034 of Fig. 2b; yearly ramp (annual steps), not a constant rate.
- Rates read off Tilmes et al. (2025) Fig. 2b (no rates from the authors in time).
- Injection at 75 hPa (Fig. 2a, 30°N/S at 18.1 km); current `cesm_numerical` preset.
- Ramp implemented in the AIDE-SAI repo as a schedule file (`INJ_SO2_SCHED_FILE`, PR #9, merged).

## Unconfirmed — PR 1 (this PR)

1. **Red/green curves of Fig. 2b are per site.** Rate at each of 30°N and 30°S, not the pair
   sum: while green is zero, red = black/2 (2030: 2.13 vs 4.27). The schedule uses the black
   total with `INJ_MIRROR=1` (equal halves), consistent with this.
2. **Equal hemispheric split.** The red member envelope in 2030–2032 is ± 0.2 Tg/yr per site;
   whether that is N/S asymmetry or member spread cannot be told from the figure.
3. **Year value = mean over the pixel columns of that calendar year** of the member-mean curve.
4. **Figure not committed**; fetched from the journal with a sha256 manifest (CC BY 4.0 would
   allow committing it).
5. **Digitiser reads the PNG with `matplotlib.pyplot.imread`** (no new pin; Pillow comes in
   through matplotlib).
6. **Run year 0 treated as spin-up** in the comparison (AIDE-SAI starts without the 2025–2029
   aerosol that CCMI 2030 carries).
