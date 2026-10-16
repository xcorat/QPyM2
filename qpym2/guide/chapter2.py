"""Utilities for computing Chapter 2 data summary values.

Provides functions to:
 - compute per-dataset and total exposures (from reduced ROOT tree or a datasets list file)
 - compute overall and per-dataset cut-flow counts and efficiencies using ROOT RDataFrame
 - read PCA efficiencies and compute exposure-weighted PSA efficiency
 - export computed values to YAML and generate a markdown snippet

Designed to be imported by a short notebook that calls these functions and writes
values to `tech_note/md/generated_chapter2_values.md` and `tech_note/values/chapter2_values.yaml`.
"""

from pathlib import Path
import yaml
import pandas as pd
import numpy as np

try:
    from ROOT import RDataFrame
except Exception:
    RDataFrame = None

from qpym2.cfg import get_config
from qpym2.utils.syst_rwp import read_ds_exposures
from qpym2.utils.syst_rwp import read_pca_eff_TMP as read_pca_eff
from qpym2.io import hdf
from qpym2.utils.logging import log


def _find_column(columns, candidates):
    for c in candidates:
        if c in columns:
            return c
    return None


def get_cfg_paths(environment='tech_note', analysis='tech_note'):
    cfg = get_config(environment=environment, analysis=analysis)
    return cfg


def compute_exposures_from_file(fname_ds):
    """Read a simple datasets list file and sum exposures.

    File format expected by read_ds_exposures(): whitespace delimited two columns: ds exposure
    """
    if not Path(fname_ds).exists():
        raise FileNotFoundError(f"Datasets file not found: {fname_ds}")

    df = read_ds_exposures(fname_ds)
    df = df.rename(columns={'exposure': 'exposure_kg_yr'})
    df['exposure_kg_yr'] = df['exposure_kg_yr'].astype(float)
    df.index.name = 'ds'
    return df


def compute_exposures_from_reduced(reduced_data_path, tree_name='outTree'):
    """Compute per-dataset exposures from a reduced ROOT file.

    This uses ROOT.RDataFrame.AsNumpy to read columns. It looks for columns
    that indicate dataset id (candidates: 'ds','Dataset','dataset') and exposure
    (candidates: 'Exposure','exposure'). Returns a DataFrame indexed by dataset with
    exposure summed (kg·yr).
    """
    if RDataFrame is None:
        raise RuntimeError('ROOT is not available in this Python environment')

    rdf = RDataFrame(tree_name, str(reduced_data_path))
    columns = [str(c) for c in rdf.GetColumnNames()]

    ds_col = _find_column(columns, ['ds', 'DS', 'Dataset', 'dataset'])
    exp_col = _find_column(columns, ['Exposure', 'exposure', 'Exp'])

    if ds_col is None:
        raise ValueError(f'No dataset column found in reduced file. Available columns: {columns[:20]}')

    if exp_col is None:
        raise ValueError('No Exposure column found in reduced file. Cannot compute exposure from reduced tree.')

    arrays = rdf.AsNumpy([ds_col, exp_col])
    arr_ds = arrays[ds_col]
    arr_exp = arrays[exp_col]

    df = pd.DataFrame({'ds': arr_ds, 'exposure_kg_yr': arr_exp})
    df = df.groupby('ds', as_index=True).sum()
    df.index = df.index.astype(int)
    return df


_NDBD_CUTOFF_U = {'u1775': 1775, 'u1740': 1740}


def compute_cutflow(reduced_data_path, tree_name='outTree', model_roi=None, cut_flow=None,
                    ds_filter=None, evar='Energy', esum_var='TotalEnergy', ndbd_cutoff=None):
    """Compute sequential cut-flow counts.

    Args:
        reduced_data_path: path to reduced ROOT file
        tree_name: ROOT tree name
        model_roi: dict with keys esum_min, esum_max, threshold
        cut_flow: list of tuples (name, expression), if not provided default flow is used
        ds_filter: optional dataset id to filter to a single ds (int)
        evar: column name for single-crystal energy (default 'Energy').
              Used to define 'e1' alias and derive 'e2 = esum - e1'.
        esum_var: column name for total (sum) energy (default 'TotalEnergy').
                  Used to define 'esum' alias.
        ndbd_cutoff: optional key for 0νββ signal region cut, e.g. 'u1740' or 'u1775'.
                     If None, the 0νββ step is omitted from the default cut flow.

    Returns:
        pd.DataFrame with Step, Cut, Events, Absolute Eff (%), Relative Eff (%)
    """
    if RDataFrame is None:
        raise RuntimeError('ROOT is not available in this Python environment')

    rdf_base = RDataFrame(tree_name, str(reduced_data_path))
    cols = [str(c) for c in rdf_base.GetColumnNames()]

    # Define 'esum' alias if not already present
    if 'esum' not in cols:
        if esum_var in cols:
            rdf_base = rdf_base.Define('esum', esum_var)
        elif 'TotalEnergy' in cols:
            rdf_base = rdf_base.Define('esum', 'TotalEnergy')

    # Define 'e1' and 'e2' aliases for threshold cuts
    # e1 = evar (first crystal energy), e2 = esum - e1
    if 'e1' not in cols:
        if evar in cols:
            rdf_base = rdf_base.Define('e1', evar)
        elif 'Energy' in cols:
            rdf_base = rdf_base.Define('e1', 'Energy')
    if 'e2' not in cols:
        rdf_base = rdf_base.Define('e2', f'{esum_var} - {evar}')

    # Define 'u' (rotated sum coordinate) for signal region cuts
    if 'u' not in cols:
        rdf_base = rdf_base.Define('u', f'{esum_var} / sqrt(2.0)')

    rdf = rdf_base
    if ds_filter is not None:
        ds_col = None
        for cand in ['ds', 'DS', 'Dataset', 'dataset']:
            if cand in cols:
                ds_col = cand
                break
        if ds_col is None:
            raise ValueError('Dataset column not found to apply ds_filter')
        rdf = rdf.Filter(f"{ds_col}=={int(ds_filter)}")

    # default cut flow
    if cut_flow is None:
        multiplicity_cut = 'Multiplicity==2' if 'Multiplicity' in cols else 'Multiplicity==2'
        if model_roi is not None:
            esum_min, esum_max, threshold = model_roi['esum_min'], model_roi['esum_max'], model_roi['threshold']
            hist_cuts = f'esum>={esum_min} && esum<={esum_max} && e1>{threshold} && e2>{threshold}'
        else:
            hist_cuts = 'esum>=2400 && esum<=2570'

        cut_flow = [
            ("Total events", "true"),
            ("Multiplicity == 2", multiplicity_cut),
            ("ROI + Threshold", hist_cuts),
        ]
        if ndbd_cutoff is not None:
            u_val = _NDBD_CUTOFF_U.get(ndbd_cutoff, int(ndbd_cutoff.lstrip('u')))
            cut_flow.append((f"0νββ region (u>{u_val})", f"u>{u_val}"))


    results = []
    rdf_current = rdf
    for step_name, cut_expr in cut_flow:
        rdf_current = rdf_current.Filter(cut_expr)
        count = rdf_current.Count().GetValue()
        results.append({'Step': step_name, 'Cut': cut_expr, 'Events': int(count)})

    df = pd.DataFrame(results)
    total = df.iloc[0]['Events']
    df['Absolute Eff (%)'] = 100 * df['Events'] / total if total > 0 else 0
    df['Relative Eff (%)'] = 100 * df['Events'] / df['Events'].shift(1)
    df.loc[0, 'Relative Eff (%)'] = 100.0
    return df


def compute_pca_weighted_efficiency(exposures_df=None):
    """Compute exposure-weighted PCA (PSA) efficiency using read_pca_eff().

    If exposures_df is provided (indexed by ds) use it to compute weighted mean; otherwise
    return the unweighted mean and the per-ds table.
    """
    pca = read_pca_eff()
    # pca index is ds numbers as ints
    pca = pca.rename_axis('ds')

    # choose a sensible column for PSA efficiency: e_mcmc (model-based)
    eff_col = 'e_mcmc'
    if eff_col not in pca.columns:
        eff_col = 'e_x2' if 'e_x2' in pca.columns else pca.columns[0]

    pca_eff = pca[[eff_col]].rename(columns={eff_col: 'psa_eff_percent'})

    if exposures_df is None:
        mean_eff = pca_eff['psa_eff_percent'].mean()
        return pca_eff, float(mean_eff)

    # Ensure ds types align
    idx = pca_eff.index.astype(int)
    exposures_idx = exposures_df.copy()
    exposures_idx.index = exposures_idx.index.astype(int)

    # Align exposures with pca table
    common = pca_eff.join(exposures_idx, how='inner')
    if 'exposure_kg_yr' not in common.columns:
        raise ValueError('exposure_kg_yr must be present in exposures_df')

    weights = common['exposure_kg_yr']
    effs = common['psa_eff_percent']

    weighted_mean = np.average(effs, weights=weights)
    return pca_eff, float(weighted_mean)


def compute_signal_efficiency_from_hdf(input_h5):
    """Compute signal reconstruction efficiency from fit input HDF5 file.

    Expects `read_signal()` to return a table with an 'integral' column giving the number
    of MC events in the mchist after selections. The number of generated events must be known:
    Some inputs include a 'generated' or 'nevt' column; if not found, we return None for efficiency.
    """
    sig = hdf.read_signal(input_h5, 'ndbd')
    integral = float(sig['integral'].iloc[0])

    # Try to detect generated number
    gen_col = None
    for c in ['generated', 'nevt', 'ngenerated', 'n_generated']:
        if c in sig.columns:
            gen_col = c
            break

    if gen_col is None:
        # try reading attributes in the file (not implemented) - return None
        return {'integral': integral, 'generated': None, 'efficiency': None}

    generated = int(sig[gen_col].iloc[0])
    eff = integral / generated if generated > 0 else None
    return {'integral': integral, 'generated': generated, 'efficiency': eff}


def export_values_yaml(md_outpath, yaml_outpath, values_dict):
    """Write a markdown snippet and YAML representation to disk.

    md_outpath: Path to markdown file to create
    yaml_outpath: Path to yaml file to create
    values_dict: dict with scalar and tabular data (tables expected as list of dicts)
    """
    md_lines = []
    md_lines.append("# Generated Chapter 2 values\n")

    # Simple table representation for cutflow and exposures
    if 'cutflow' in values_dict:
        md_lines.append("## Cut flow\n")
        md_lines.append("| Step | Events | Absolute Eff (%) | Relative Eff (%) |")
        md_lines.append("|------|--------:|-----------------:|-----------------:|")
        for r in values_dict['cutflow']:
            md_lines.append(f"| {r['Step']} | {r['Events']:,} | {r['Absolute Eff (%)']:.2f} | {r['Relative Eff (%)']:.2f} |")
        md_lines.append('\n')

    if 'exposures' in values_dict:
        md_lines.append('## Dataset exposures (kg·yr)\n')
        md_lines.append('| Dataset | Exposure (kg·yr) |')
        md_lines.append('|--------:|-----------------:|')
        for ds, exp in values_dict['exposures'].items():
            md_lines.append(f"| {ds} | {exp:.3f} |")
        md_lines.append('\n')

    Path(md_outpath).parent.mkdir(parents=True, exist_ok=True)
    with open(md_outpath, 'w') as f:
        f.write('\n'.join(md_lines))

    # YAML
    Path(yaml_outpath).parent.mkdir(parents=True, exist_ok=True)
    with open(yaml_outpath, 'w') as f:
        yaml.safe_dump(values_dict, f)

    return md_outpath, yaml_outpath
