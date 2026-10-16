import numpy as np
import pandas as pd

from qpym2.fit.create_m2uv import get_shift, get_shiftname, get_systname
from qpym2.cfg.model_utils import build_model_dir, get_shift_suffix


def make_priors(comp_cfg, mctable, mksamples, sigma_scale=1.0, use_sample_std=False):
    """Generate prior distributions from JAGS MC systematic samples.

    For each background component, computes the expected number of events
    across all MC chain samples and derives the mean and std used as
    Gaussian prior parameters in the PyMC model.

    Prior sigma = sqrt(mean_nevs) * sigma_scale  (Poisson-based width),
    or std(nevs_samples) when use_sample_std=True (syst mode).

    Args:
        comp_cfg: List of component dicts, each with 'filter_regex' or
            'filter_regex_inverse' selecting rows from mctable.
        mctable: DataFrame indexed by source name, with a 'mchist' column
            containing the MC histogram arrays.
        mksamples: DataFrame of JAGS chain samples (columns = sample index,
            rows = source names). Typically from ``read_mkchain(...).T``.
        sigma_scale: Multiplicative scale applied to sqrt(N) prior width
            (default 1.0). Use 5.0 for '5sig' mode, 10.0 for '10sig'.
            Ignored when use_sample_std=True.
        use_sample_std: If True, use the empirical std of nevs_samples as the
            prior sigma (syst mode). If False, use sqrt(mean)*sigma_scale.

    Returns:
        comp_cfg with 'mean' and 'std' keys added to each component dict.
    """
    print(f"\nMake priors - MC table: {mctable.shape}, "
          f"MK samples: {mksamples.shape}, Components: {len(comp_cfg)}, "
          f"sigma_scale={sigma_scale}  use_sample_std={use_sample_std}")

    def _get_component_table(table, comp):
        if 'filter_regex' in comp:
            return table.filter(regex=comp['filter_regex'], axis=0)
        elif 'filter_regex_inverse' in comp:
            filtered = table.filter(regex=comp['filter_regex_inverse'], axis=0)
            return table.drop(filtered.index)

    tb = mctable.copy()
    tbs = pd.concat([tb, mksamples], axis=1, join='inner')

    component_tables = [_get_component_table(tbs, comp) for comp in comp_cfg]

    def _get_comp_nevs(ct, ismpl):
        pdf = ct['mchist'] * ct[ismpl].values
        return np.sum(pdf.sum())

    nevs_samples = np.array([[_get_comp_nevs(ct, ismpl)
                              for ct in component_tables]
                             for ismpl in mksamples.columns])

    print(f"Prior nevs_samples shape: {nevs_samples.shape}")
    means = np.mean(nevs_samples, axis=0)
    if use_sample_std:
        stds = nevs_samples.std(axis=0) * sigma_scale
        print(f"Prior stds (empirical): {stds}")
    else:
        stds = np.sqrt(means) * sigma_scale
    for i, comp in enumerate(comp_cfg):
        comp['std'] = stds[i]
        comp['mean'] = means[i]

    return comp_cfg


def model_cfg(hist_model, name='model', co_cut=None, ndbd_cutoff=None,
              co_shift=0, ndbd_shift=0,
              skip_roi=False, blinded=False, smooth=True, 
              common_cuts=[], mc_cuts=None,
              sample_mkchain=0,
              co_cuts_dict=None, ndbd_cuts_dict=None):
    """Create model configuration.
    
    Args:
        co_cuts_dict: Dictionary of co_cut -> cut_string (from config)
        ndbd_cuts_dict: Dictionary of ndbd_cutoff -> cut_string (from config)
        Other parameters: see inline usage
    """
    mname = name
    ndbd_cuts = []
    data_cuts = []

    if co_cut and co_cut not in (None, 'None', 'null'):
        if co_cuts_dict is None:
            # Fallback: lazy load from config if not provided
            from qpym2.cfg.hydra_adapter import get_config
            _cfg = get_config()
            co_cuts_dict = dict(_cfg.cuts.co_cuts)
        mname = mname + '_' + co_cut
        common_cuts.append(co_cuts_dict[co_cut])
    if ndbd_cutoff and ndbd_cutoff not in (None, 'None', 'null'):
        if ndbd_cuts_dict is None:
            # Fallback: lazy load from config if not provided
            from qpym2.cfg.hydra_adapter import get_config
            _cfg = get_config()
            ndbd_cuts_dict = dict(_cfg.cuts.ndbd_cutoff)
        mname = mname + '_' + ndbd_cutoff
        ndbd_cuts.append(ndbd_cuts_dict[ndbd_cutoff])
    if co_shift:
        mname = mname + '_' + get_shiftname(get_shift(co_shift))
    if ndbd_shift:
        # mname = mname + '_' + get_shiftname(get_shift(co_shift))
        mname = mname + '_' + get_systname( {
            'ndbd_shift': ndbd_shift
        })
    if not smooth:
        mname = mname + '_nosmooth'

    # dir_name: directory without shifts; shift_suffix: filename suffix
    dir_name = build_model_dir(
        base_name=name, co_cut=co_cut if co_cut not in (None, 'None', 'null') else None,
        ndbd_cutoff=ndbd_cutoff if ndbd_cutoff not in (None, 'None', 'null') else None,
        smooth=smooth)
    shift_suffix = get_shift_suffix(co_shift, ndbd_shift)

    if skip_roi:
        hist_model.skip_roi = skip_roi
    
    return {
        'name': mname,
        'dir_name': dir_name,
        'shift_suffix': shift_suffix,
        'cuts': {
            'common': '&&'.join(common_cuts) if common_cuts else None,
            'mc_only': '&&'.join(mc_cuts) if mc_cuts else None,
            'ndbd_only': '&&'.join(ndbd_cuts) if ndbd_cuts else None,
            'data_only': '&&'.join(data_cuts) if data_cuts else None,
        },
        'h2uv_model': hist_model,
        'fit_pars': {
            'co_shift': co_shift
        },
        'syst_pars': {
            'ndbd_shift': ndbd_shift
        },
        'skip_roi': skip_roi,
        'blinded': blinded,
        'smooth': smooth,
        'sample_mkchain': 0,
    }
