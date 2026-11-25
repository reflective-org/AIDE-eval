#!/usr/bin/env bash
# Download Tilmes et al. (2025), ACP 25, 6001-6023, Fig. 2 (CC BY 4.0) into the
# gitignored AIDE-SAI/data/tilmes2025/ and check it against tilmes2025_fig2.sha256.
# The figure is the source of the yearly SAI injection rates (digitize_tilmes2025_fig2b.py).
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
dest="$here/../data/tilmes2025"
mkdir -p "$dest"
url=https://acp.copernicus.org/articles/25/6001/2025/acp-25-6001-2025-f02-web.png
[ -f "$dest/acp-25-6001-2025-f02-web.png" ] || curl -sSfL -o "$dest/acp-25-6001-2025-f02-web.png" "$url"
cd "$dest" && sha256sum -c "$here/tilmes2025_fig2.sha256"
