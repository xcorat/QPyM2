"""Systematic effects pipeline: read → combine → smear → rate → save → plot.

Refactored from the per-systematic notebooks
(``06-lss.ipynb``, ``06-bmmc_ndbd.ipynb``, ``06-bias_corr_syst.ipynb``,
``06-syst-rate-new.ipynb``) into reusable components.

Typical usage::

    from qpym2.syst import (
        SystDelta, compute_delta, combine_deltas,
        smear_posterior, to_rate, rate_stats,
        save_syst_delta, load_syst_delta,
        save_combined_result, load_combined_result,
    )
    from qpym2.syst import plots as syst_plots
"""

from qpym2.syst.core import SystDelta, compute_deltas, combine_deltas,  compute_1sigma_edge_delta
from qpym2.syst.smear import smear_posterior
from qpym2.syst.rate import F130, N_TEO2, to_rate, rate_stats
from qpym2.syst.io import (
    save_syst_delta, load_syst_delta,
    save_syst_result, load_syst_result, post_stats,
    save_combined_result, load_combined_result,
    cs_tag,
)
from qpym2.syst import plots

__all__ = [
    'SystDelta', 'compute_deltas', 'combine_deltas', 'compute_1sigma_edge_delta',
    'smear_posterior',
    'F130', 'N_TEO2', 'to_rate', 'rate_stats',
    'save_syst_delta', 'load_syst_delta',
    'save_syst_result', 'load_syst_result', 'post_stats',
    'save_combined_result', 'load_combined_result',
    'cs_tag', 'plots',
]
