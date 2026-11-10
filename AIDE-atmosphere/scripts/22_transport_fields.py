"""
22 - Read a source's transport fields and compute the transport diagnostics.

The first stage of the transport scorecard (docs/TRANSPORT_SCORECARD.md). A
source is read through transport_common.READERS into the field contract of
section 1 of that document, then reduced to the metrics of groups A-D, each
beside its reference.

  22_transport_fields.py --inspect NAME      shapes, units and ranges of the
                                             contract fields, nothing written
"""

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import transport_common as T

YEARS = range(2006, 2011)


def inspect(name):
    """Print what the reader hands the diagnostics, so a new source can be checked
    against the contract before anything is scored."""
    spec = T.SOURCES[name]
    F = T.load_source(name, YEARS)
    p, lat, days = F["plev"], F["lat"], F["days"]
    print(f"{name}: {spec['long_name']}  reader {spec['reader']}  {spec['path']}")
    print(f"  days  {days.size}  {days[0]} .. {days[-1]}")
    print(f"  plev  {p.size} levels  {p.min():.4g} .. {p.max():.4g} hPa, ascending "
          f"{bool(np.all(np.diff(p) > 0))}; inside 5-150 hPa: {int(((p >= 5) & (p <= 150)).sum())}")
    print(f"  lat   {lat.size}  {lat.min():.3f} .. {lat.max():.3f}, ascending "
          f"{bool(np.all(np.diff(lat) > 0))}")
    for k in T.CONTRACT + T.OPTIONAL:
        if k not in F:
            print(f"  {k:7s} absent")
            continue
        x = F[k]
        print(f"  {k:7s} {str(x.shape):18s} {str(x.dtype):8s} [{T.UNITS[k]:5s}] "
              f"min {np.nanmin(x):11.4g}  mean {np.nanmean(x, dtype='f8'):11.4g}  "
              f"max {np.nanmax(x):11.4g}  NaN {int(np.isnan(x).sum())}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1].strip())
    ap.add_argument("--inspect", metavar="NAME", help="print the contract fields of a source")
    a = ap.parse_args()
    if a.inspect:
        inspect(a.inspect)


if __name__ == "__main__":
    main()
