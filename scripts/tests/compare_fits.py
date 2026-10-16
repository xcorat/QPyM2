"""
Compare posterior modes for all components across co-shifts between two fit sets:
  - data/fits/model_gamma10_v150_u1740       (nominal)
  - data_syst/fits/model_gamma10_v150_u1740  (syst)

Run from QPyM2/ root:
  python scripts/tests/compare_fits.py
"""
import sys
import os
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

from qpym2.fit.results import read_posterior_stats

SHIFTS = [0.0, 0.3, 0.5, 0.8, 1.0, 1.2, 1.5]
MODEL_PARS_BASE = dict(co_cut='gamma10_v150', ndbd_cutoff='u1740', smooth=True)

DIR_A = 'data/fits'
DIR_B = 'data_syst/fits'


def get_modes(parent_dir, shift):
    mp = {**MODEL_PARS_BASE, 'co_shift': shift}
    try:
        df = read_posterior_stats(parent_dir, mp, hm=None)
        return df['mode']
    except Exception as e:
        print(f"  [WARN] {parent_dir}  shift={shift}: {e}")
        return None


def main():
    rows = []
    for shift in SHIFTS:
        modes_a = get_modes(DIR_A, shift)
        modes_b = get_modes(DIR_B, shift)
        if modes_a is None or modes_b is None:
            continue
        common = modes_a.index.intersection(modes_b.index)
        for var in common:
            a, b = modes_a[var], modes_b[var]
            diff = b - a
            rel = (diff / a * 100) if a != 0 else np.nan
            rows.append(dict(shift=shift, var=var,
                             mode_data=a, mode_syst=b,
                             diff=diff, rel_diff_pct=rel))

    df = pd.DataFrame(rows).set_index(['shift', 'var'])
    pd.set_option('display.float_format', '{:.4g}'.format)
    pd.set_option('display.max_rows', 300)

    for var in df.index.get_level_values('var').unique():
        sub = df.xs(var, level='var')[['mode_data', 'mode_syst', 'diff', 'rel_diff_pct']]
        print(f"\n=== {var} ===")
        print(sub.to_string())


if __name__ == '__main__':
    os.chdir(os.path.join(os.path.dirname(__file__), '..', '..'))
    main()