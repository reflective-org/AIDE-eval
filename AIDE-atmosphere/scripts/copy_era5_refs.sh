#!/usr/bin/env bash
# Copy the ERA5 reference files the transport scorecard reads into AIDE-atmosphere/data/era5/
# (gitignored) and check them against the committed manifest era5_refs.sha256.
#
#   bash scripts/copy_era5_refs.sh                    copy, then verify
#   bash scripts/copy_era5_refs.sh --write-manifest   copy, then (re)write the manifest
#
# Only the files a metric reads (traced in jcm-strat scripts/teq_score.py and
# scripts/variability_score.py), for the scored windows 2006-2010 (phase-29 set) and
# 2011-2019 (AL7h, phase 31):
#   era5_zm_monthly_<y>.nc  monthly zonal-mean uzm, tzm, 25 levels 1-1000 hPa, 0.25 deg   C1, C2  (teq_update.ERA5_ZM)
#   era5_<y>.npz            CDS daily zonal means, 17 levels 1-250 hPa  D1-D7  (--era5-dir .../full)
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA="$HERE/../data"; DEST="$DATA/era5"; MANIFEST="$HERE/era5_refs.sha256"
SRC_ZM=/data/JCM_stripped/jcm-strat/cache/era5_ref
SRC_DAILY=/data/JCM_stripped/jcm-strat/cache/variability_cds/full
# jcm-strat-phase31/cache/{era5_ref,variability_cds} are symlinks to these two directories
mkdir -p "$DEST"
for y in $(seq 2006 2019); do
  cp --preserve=timestamps "$SRC_ZM/era5_zm_monthly_$y.nc" "$DEST/"
  cp --preserve=timestamps "$SRC_DAILY/era5_$y.npz" "$DEST/"
done
cd "$DATA"
if [ "${1:-}" = "--write-manifest" ]; then
  sha256sum era5/era5_zm_monthly_20*.nc era5/era5_20*.npz > "$MANIFEST"
  echo "wrote $MANIFEST"
else
  sha256sum -c "$MANIFEST"
fi
