#!/usr/bin/env python3
"""Compare two staging directories for event count differences.
Validates that extra events in staging_150 are due to a lower energy threshold.

Usage:
    python tmp/check_staging_diff.py
"""
import pathlib
from ROOT import RDataFrame as RDF

OLD = pathlib.Path("/global/cfs/projectdirs/cuore/scratch/xcorat/data/MC/staging")
NEW = pathlib.Path("/global/cfs/projectdirs/cuore/scratch/xcorat/data/MC/staging_150")
TREENAME = "uvtree"

def count_and_min_u0(fpath: pathlib.Path):
    rdf = RDF(TREENAME, str(fpath))
    n   = rdf.Count().GetValue()
    mn  = rdf.Min("u0").GetValue() if n > 0 else float("nan")
    return n, mn

def main():
    old_files = {f.name: f for f in OLD.glob("*.root")}
    new_files = {f.name: f for f in NEW.glob("*.root")}

    only_old  = sorted(old_files.keys() - new_files.keys())
    only_new  = sorted(new_files.keys() - old_files.keys())
    both      = sorted(old_files.keys() & new_files.keys())

    if only_old:
        print(f"\n=== Only in OLD ({len(only_old)}) ===")
        for n in only_old:
            print(f"  {n}")

    if only_new:
        print(f"\n=== Only in NEW ({len(only_new)}) ===")
        for n in only_new:
            print(f"  {n}")

    print(f"\n=== In both ({len(both)}) ===")
    print(f"{'file':<50}  {'old_n':>10}  {'new_n':>10}  {'delta':>8}  {'old_min_u0':>12}  {'new_min_u0':>12}  valid_threshold")

    for name in both:
        old_n, old_min = count_and_min_u0(old_files[name])
        new_n, new_min = count_and_min_u0(new_files[name])
        delta = new_n - old_n

        if delta != 0:
            # Extra events should have lower u0 than old minimum
            valid = (new_min < old_min) if delta > 0 else "N/A (fewer)"
            print(f"  {name:<48}  {old_n:>10}  {new_n:>10}  {delta:>+8}  {old_min:>12.4f}  {new_min:>12.4f}  {valid}")
        else:
            print(f"  {name:<48}  {old_n:>10}  {new_n:>10}  {'(same)':>8}")

if __name__ == "__main__":
    main()