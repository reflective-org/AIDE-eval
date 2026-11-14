"""
evaluate - Score result sets end to end: 01_fields -> 02_scorecard -> 03_figures.

  evaluate.py SET [SET ...]     the named sets (sources/<SET>.json)
  evaluate.py                   every set in sources/

Each stage runs as its own process with this interpreter, in order, and stops the chain
on the first failure. 02 rewrites validation_results/SCORING_STANDARD.md and
SCORING_EXTENSIVE.md over every registered set, so evaluating one set updates the
cross-set summaries. The ERA5 copies must be in place first (copy_era5_refs.sh).
"""

import argparse
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import strat_eval_common as SE

STAGES = ("01_fields.py", "02_scorecard.py", "03_figures.py")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1].strip())
    ap.add_argument("sets", nargs="*", help="result sets to evaluate (default: every set)")
    sets = ap.parse_args().sets or SE.set_names()
    for s in sets:
        SE.load_set(s)                     # fail on an unknown or malformed set before any work
    for stage in STAGES:
        print(f"== {stage} {' '.join(sets)}", flush=True)
        subprocess.run([sys.executable, os.path.join(HERE, stage), *sets], check=True)


if __name__ == "__main__":
    main()
