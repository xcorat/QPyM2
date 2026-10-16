"""
Integration tests for prior generation and fit setup in run_fits.

These tests require actual data files and a configured environment.
Run with: pytest -m integration qpym2/tests/test_run_fits.py
"""

import pytest

from qpym2.cfg.hydra_adapter import get_config
from qpym2.cfg.config_builders import (
    get_jags_path, get_data_names, get_fit_dir, get_h2uv_model, get_comp_cfg
)
from qpym2.io.bkg_model import read_mkchain
from qpym2.fit.params import model_cfg
from qpym2.fit import fit as fit_module
from qpym2.io.fit import read_input
from qpym2.fit.params import make_priors


@pytest.mark.integration
def test_prior_generation_and_fit(environment='tech_note'):
    """Verify prior generation from MC samples and subsequent fit setup.

    Checks:
    - mkchain samples produce valid priors via make_priors()
    - fit.setup_fit() succeeds with generated priors
    - fit_comps contains a non-None 'prior' key after the run
    """
    cfg = get_config(environment=environment, analysis=environment)

    jags_fpath = get_jags_path(cfg)
    data_name, _ = get_data_names(cfg)
    fit_dir = get_fit_dir(cfg)
    h2uv_model = get_h2uv_model(cfg)
    comp_cfg_list = get_comp_cfg(cfg)

    mcfg = model_cfg(h2uv_model, co_cut='gamma10_v150', ndbd_cutoff='u1740', co_shift=0.8)
    fitdir = f"{fit_dir}/{mcfg['dir_name']}"
    hm = mcfg['h2uv_model']
    shift_suffix = mcfg.get('shift_suffix', '')

    mctable, signal_df, data, edges = read_input(f"{fitdir}/input{shift_suffix}.h5", dataname=data_name)
    mksamples = read_mkchain(jags_fpath).sample(n=100).T

    comp_cfg_new = (
        make_priors(comp_cfg_list, mctable, mksamples)
        if mksamples is not None
        else comp_cfg_list
    )

    m2model, fit_comps = fit_module.setup_fit(
        mctable, data, signal_df, hm,
        smooth=hm.smooth,
        comp_cfg=comp_cfg_new,
        edges=edges,
    )
    trace = fit_module.run_fit(m2model, single_core=False, chains=4, draws=5000)

    assert 'prior' in fit_comps, "fit_comps must contain 'prior'"
    assert fit_comps['prior'] is not None
    trace.posterior.attrs['priors'] = fit_comps['prior']
