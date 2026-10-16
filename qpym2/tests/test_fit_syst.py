"""Unit tests for qpym2.fit.syst.

Tests are split into:
  - Pure-numpy tests (no ROOT, run in CI)
  - Integration tests marked with @pytest.mark.integration (require ROOT RDFs)

Run pure tests:
    pytest qpym2/tests/test_fit_syst.py

Run all including integration:
    pytest -m integration qpym2/tests/test_fit_syst.py
"""

import types

import numpy as np
import pandas as pd
import pytest

from qpym2.fit.syst import (
    _DEFAULT_LSS_SYST_SAMPLES,
    _apply_sscl_transform,
    _empty_hist,
    _hist_edges,
    _hist_from_np,
    get_lss_syst_defs,
    make_lss_mctable,
    make_lss_mctable_np,
)
from qpym2.utils.hists import HistType


# ---------------------------------------------------------------------------
# Shared fixtures / helpers
# ---------------------------------------------------------------------------

def _make_hm(umin=2540, umax=2700, vmin=0, vmax=160,
             binsize=2.0, threshold=-50.0):
    """Minimal histogram-model namespace."""
    hm = types.SimpleNamespace(
        umin=umin, umax=umax,
        vmin=vmin, vmax=vmax,
        binsize=binsize,
        threshold=threshold,
        skip_roi=None,
        hist_type=HistType.h2uv,
    )
    return hm


def _make_lss_params(n_ds=3, n_ch=10, ds_start=3601, rng=None):
    """Minimal lss_params dict for sscl tests."""
    if rng is None:
        rng = np.random.default_rng(0)
    sigma0_2d = rng.uniform(0.01, 0.1, size=(n_ds, n_ch))
    delta_p1_1r = rng.uniform(0.001, 0.01, size=(n_ds,))
    delta_p1_1l = -delta_p1_1r
    # q_pars not needed for sscl
    return {
        'sigma0_2d': sigma0_2d,
        'delta_p1_1r': delta_p1_1r,
        'delta_p1_1l': delta_p1_1l,
        'q_pars_arr': [],
        'q_pars_sample_1r': [],
        'q_pars_sample_1l': [],
    }


class _MockRDF:
    """Fake RDF: wraps pre-built numpy arrays, exposes AsNumpy."""

    def __init__(self, u0, v0, ch1, ch2, ds_arr):
        self._data = {
            'u0': u0.astype(np.float32),
            'v0': v0.astype(np.float32),
            'ch1': ch1.astype(np.int32),
            'ch2': ch2.astype(np.int32),
            'ds':  ds_arr.astype(np.int32),
        }

    def AsNumpy(self, columns):
        return {c: self._data[c] for c in columns}


# ---------------------------------------------------------------------------
# get_lss_syst_defs
# ---------------------------------------------------------------------------

class TestGetLssSystDefs:
    def test_length(self):
        datasets = [3601, 3602, 3603]
        result = get_lss_syst_defs(_DEFAULT_LSS_SYST_SAMPLES, datasets)
        assert len(result) == len(_DEFAULT_LSS_SYST_SAMPLES) * len(datasets)

    def test_structure(self):
        datasets = [3601]
        result = get_lss_syst_defs(_DEFAULT_LSS_SYST_SAMPLES, datasets)
        for sname, scfg, ds in result:
            assert sname in _DEFAULT_LSS_SYST_SAMPLES
            assert scfg is _DEFAULT_LSS_SYST_SAMPLES[sname]
            assert ds == 3601

    def test_ordering_sname_outer(self):
        """Outer loop is over snames, inner over datasets."""
        datasets = [3601, 3602]
        result = get_lss_syst_defs(_DEFAULT_LSS_SYST_SAMPLES, datasets)
        snames = [r[0] for r in result]
        # Each sname should appear in a block of len(datasets)
        for i in range(0, len(result), len(datasets)):
            block = snames[i:i + len(datasets)]
            assert len(set(block)) == 1


# ---------------------------------------------------------------------------
# _empty_hist / _hist_edges
# ---------------------------------------------------------------------------

class TestHistHelpers:
    def test_empty_hist_shape(self):
        hm = _make_hm(umin=0, umax=10, vmin=0, vmax=4, binsize=2.0)
        h = _empty_hist(hm)
        assert h.shape == (5, 2)
        assert h.sum() == 0

    def test_hist_edges_count(self):
        hm = _make_hm(umin=0, umax=10, vmin=0, vmax=4, binsize=2.0)
        u_edges, v_edges = _hist_edges(hm)
        assert len(u_edges) == 6   # nbins + 1
        assert len(v_edges) == 3
        assert u_edges[0] == pytest.approx(0.0)
        assert u_edges[-1] == pytest.approx(10.0)


# ---------------------------------------------------------------------------
# _hist_from_np
# ---------------------------------------------------------------------------

class TestHistFromNp:
    def _make_uv(self, hm, n=200, rng=None):
        if rng is None:
            rng = np.random.default_rng(42)
        u = rng.uniform(hm.umin + 1, hm.umax - 1, n)
        v = rng.uniform(1.0, hm.vmax - 1, n)
        return u, v

    def test_returns_empty_for_none(self):
        hm = _make_hm()
        h = _hist_from_np(None, hm)
        assert h.sum() == 0

    def test_shape_matches_hm(self):
        hm = _make_hm(umin=2540, umax=2700, vmin=0, vmax=160, binsize=4.0)
        u, v = self._make_uv(hm)
        h = _hist_from_np({'u': u, 'v': v}, hm)
        nbins_u = round((hm.umax - hm.umin) / hm.binsize)
        nbins_v = round((hm.vmax - hm.vmin) / hm.binsize)
        assert h.shape == (nbins_u, nbins_v)

    def test_events_outside_window_excluded(self):
        hm = _make_hm(umin=100, umax=200, vmin=0, vmax=50, binsize=10.0,
                      threshold=-999)
        # All u outside [umin, umax]
        u = np.full(50, 300.0)
        v = np.full(50, 10.0)
        h = _hist_from_np({'u': u, 'v': v}, hm)
        assert h.sum() == 0

    def test_skip_roi_removes_events(self):
        hm = _make_hm(umin=0, umax=400, vmin=0, vmax=400, binsize=10.0,
                      threshold=-999)
        hm.skip_roi = (100 * np.sqrt(2), 200 * np.sqrt(2))  # ROI in esum units
        # Place all events at u ≈ 150 (inside ROI)
        u = np.full(100, 150.0)
        v = np.zeros(100)
        h = _hist_from_np({'u': u, 'v': v}, hm)
        assert h.sum() == 0

    def test_events_counted_correctly(self):
        hm = _make_hm(umin=0, umax=100, vmin=0, vmax=100, binsize=10.0,
                      threshold=-999)
        u = np.array([5.0, 15.0, 55.0])
        v = np.array([5.0, 5.0, 5.0])
        h = _hist_from_np({'u': u, 'v': v}, hm)
        assert h.sum() == pytest.approx(3)


# ---------------------------------------------------------------------------
# _apply_sscl_transform  (no ROOT needed — uses _MockRDF)
# ---------------------------------------------------------------------------

class TestApplySsclTransform:
    def _make_evtrees(self, n_events=300, n_ch=10, ds_start=3601,
                      n_ds=3, rng=None):
        if rng is None:
            rng = np.random.default_rng(7)
        u0 = rng.uniform(2540, 2700, n_events)
        v0 = rng.uniform(0.1, 50, n_events)
        ch1 = rng.integers(1, n_ch + 1, n_events)
        ch2 = rng.integers(1, n_ch + 1, n_events)
        ds_arr = rng.integers(ds_start, ds_start + n_ds, n_events)
        return {
            'comp_a': _MockRDF(u0, v0, ch1, ch2, ds_arr),
            'comp_b': _MockRDF(u0 * 0.9, v0, ch1, ch2, ds_arr),
            'none_comp': None,
        }

    def test_none_component_passes_through(self):
        evtrees = {'c': None}
        lss = _make_lss_params()
        result = _apply_sscl_transform(
            evtrees, {'type': 'sscl', 'spt': 1}, ds=3601,
            lss_params=lss, ds_start=3601, rng=np.random.default_rng(0),
        )
        assert result['c'] is None

    def test_output_has_same_keys(self):
        evtrees = self._make_evtrees()
        lss = _make_lss_params()
        result = _apply_sscl_transform(
            evtrees, {'type': 'sscl', 'spt': 1}, ds=3601,
            lss_params=lss, ds_start=3601, rng=np.random.default_rng(0),
        )
        assert set(result.keys()) == set(evtrees.keys())

    def test_output_arrays_same_length(self):
        n = 200
        evtrees = self._make_evtrees(n_events=n)
        lss = _make_lss_params()
        result = _apply_sscl_transform(
            evtrees, {'type': 'sscl', 'spt': 1}, ds=3601,
            lss_params=lss, ds_start=3601, rng=np.random.default_rng(0),
        )
        for key, val in result.items():
            if val is None:
                continue
            assert len(val['u']) == n
            assert len(val['v']) == n

    def test_non_target_ds_unchanged(self):
        """Events from ds != target_ds must keep u0, v0."""
        rng = np.random.default_rng(99)
        n = 100
        u0 = rng.uniform(2540, 2700, n)
        v0 = rng.uniform(0.1, 50, n)
        ch1 = rng.integers(1, 11, n)
        ch2 = rng.integers(1, 11, n)
        # All events belong to ds=3602, target is 3601
        ds_arr = np.full(n, 3602, dtype=np.int32)
        evtrees = {'c': _MockRDF(u0, v0, ch1, ch2, ds_arr)}
        lss = _make_lss_params()
        result = _apply_sscl_transform(
            evtrees, {'type': 'sscl', 'spt': 1}, ds=3601,
            lss_params=lss, ds_start=3601, rng=np.random.default_rng(0),
        )
        np.testing.assert_allclose(result['c']['u'], u0)
        np.testing.assert_allclose(result['c']['v'], v0)

    def test_negative_spt_uses_1l_key(self):
        """spt < 0 should use delta_p1_1l (negative shift)."""
        rng = np.random.default_rng(5)
        n = 200
        u0 = np.full(n, 2600.0)
        v0 = np.full(n, 10.0)
        ch1 = np.ones(n, dtype=np.int32)
        ch2 = np.ones(n, dtype=np.int32)
        ds_arr = np.full(n, 3601, dtype=np.int32)
        evtrees = {'c': _MockRDF(u0, v0, ch1, ch2, ds_arr)}

        lss = _make_lss_params(rng=rng)
        # Force 1r and 1l to have opposite signs
        lss['delta_p1_1r'] = np.array([+0.05, 0.0, 0.0])
        lss['delta_p1_1l'] = np.array([-0.05, 0.0, 0.0])

        res_r = _apply_sscl_transform(
            evtrees, {'type': 'sscl', 'spt': 1}, ds=3601,
            lss_params=lss, ds_start=3601, rng=np.random.default_rng(1),
        )
        res_l = _apply_sscl_transform(
            evtrees, {'type': 'sscl', 'spt': -1}, ds=3601,
            lss_params=lss, ds_start=3601, rng=np.random.default_rng(1),
        )
        # Mean shifts should differ in sign
        du_r = res_r['c']['u'].mean() - u0.mean()
        du_l = res_l['c']['u'].mean() - u0.mean()
        assert du_r * du_l < 0 or (abs(du_r) > 0 and abs(du_l) > 0)


# ---------------------------------------------------------------------------
# make_lss_mctable  (sscl path only, no ROOT)
# ---------------------------------------------------------------------------

class TestMakeLssMctable:
    def _make_fit_cfg(self):
        hm = _make_hm(umin=2540, umax=2700, vmin=0, vmax=160,
                      binsize=4.0, threshold=-999.0)
        return {'h2uv_model': hm, 'cuts': {}}

    def _make_mcfit_table(self, comps=('bkg_a', 'bkg_b')):
        return pd.DataFrame(
            {'integral': [1.0] * len(comps)},
            index=list(comps),
        )

    def _make_evtrees(self, comps, n=300, n_ch=10, ds_target=3601,
                      ds_start=3601, rng=None):
        if rng is None:
            rng = np.random.default_rng(12)
        evtrees = {}
        for c in comps:
            u0 = rng.uniform(2545, 2695, n)
            v0 = rng.uniform(1.0, 50.0, n)
            ch1 = rng.integers(1, n_ch + 1, n)
            ch2 = rng.integers(1, n_ch + 1, n)
            ds_arr = np.full(n, ds_target, dtype=np.int32)
            evtrees[c] = _MockRDF(u0, v0, ch1, ch2, ds_arr)
        # ndbd component
        u0 = rng.uniform(2545, 2695, n)
        v0 = rng.uniform(1.0, 50.0, n)
        evtrees['ndbd'] = _MockRDF(u0, v0,
                                   rng.integers(1, n_ch + 1, n),
                                   rng.integers(1, n_ch + 1, n),
                                   np.full(n, ds_target, dtype=np.int32))
        return evtrees

    def test_returns_tuple_of_df_and_signal(self):
        comps = ['bkg_a', 'bkg_b']
        fit_cfg = self._make_fit_cfg()
        mcfit_table = self._make_mcfit_table(comps)
        evtrees = self._make_evtrees(comps)
        lss = _make_lss_params()
        scfg = {'type': 'sscl', 'spt': 1}

        mctable, signal_df = make_lss_mctable(
            evtrees, mcfit_table, fit_cfg, scfg,
            ds=3601, lss_params=lss, ds_start=3601,
            rng=np.random.default_rng(0),
        )
        assert isinstance(mctable, pd.DataFrame)
        assert isinstance(signal_df, pd.DataFrame)

    def test_mctable_has_mchist_and_integral(self):
        comps = ['bkg_a']
        fit_cfg = self._make_fit_cfg()
        mcfit_table = self._make_mcfit_table(comps)
        evtrees = self._make_evtrees(comps)
        lss = _make_lss_params()
        scfg = {'type': 'sscl', 'spt': 1}

        mctable, _ = make_lss_mctable(
            evtrees, mcfit_table, fit_cfg, scfg,
            ds=3601, lss_params=lss, ds_start=3601,
            rng=np.random.default_rng(0),
        )
        assert 'mchist' in mctable.columns
        assert 'integral' in mctable.columns
        assert mctable.index.tolist() == comps

    def test_signal_df_index_is_ndbd(self):
        comps = ['bkg_a']
        fit_cfg = self._make_fit_cfg()
        mcfit_table = self._make_mcfit_table(comps)
        evtrees = self._make_evtrees(comps)
        lss = _make_lss_params()

        _, signal_df = make_lss_mctable(
            evtrees, mcfit_table, fit_cfg, {'type': 'sscl', 'spt': 1},
            ds=3601, lss_params=lss, ds_start=3601,
            rng=np.random.default_rng(0),
        )
        assert list(signal_df.index) == ['ndbd']

    def test_unknown_syst_type_raises(self):
        comps = ['bkg_a']
        fit_cfg = self._make_fit_cfg()
        mcfit_table = self._make_mcfit_table(comps)
        evtrees = self._make_evtrees(comps)
        lss = _make_lss_params()

        with pytest.raises(ValueError, match="Unknown syst type"):
            make_lss_mctable(
                evtrees, mcfit_table, fit_cfg, {'type': 'bad'},
                ds=3601, lss_params=lss, ds_start=3601,
            )

    def test_integral_matches_hist_sum(self):
        comps = ['bkg_a', 'bkg_b']
        fit_cfg = self._make_fit_cfg()
        mcfit_table = self._make_mcfit_table(comps)
        evtrees = self._make_evtrees(comps)
        lss = _make_lss_params()

        mctable, signal_df = make_lss_mctable(
            evtrees, mcfit_table, fit_cfg, {'type': 'sscl', 'spt': 1},
            ds=3601, lss_params=lss, ds_start=3601,
            rng=np.random.default_rng(0),
        )
        for comp in comps:
            assert mctable.at[comp, 'integral'] == pytest.approx(
                mctable.at[comp, 'mchist'].sum()
            )
        assert signal_df.at['ndbd', 'integral'] == pytest.approx(
            signal_df.at['ndbd', 'mchist'].sum()
        )


# ---------------------------------------------------------------------------
# make_lss_mctable_np  (production path used by pool workers)
# ---------------------------------------------------------------------------

class TestMakeLssMctableNp:
    def _make_fit_cfg(self, cuts=None):
        hm = _make_hm(umin=2540, umax=2700, vmin=0, vmax=160,
                      binsize=4.0, threshold=-999.0)
        return {'h2uv_model': hm, 'cuts': cuts or {}}

    def _make_mcfit_table(self, comps=('bkg_a', 'bkg_b')):
        return pd.DataFrame(
            {'integral': [1.0] * len(comps)},
            index=list(comps),
        )

    def _make_evdata(self, comps, n=300, n_ch=10, ds_target=3601,
                     ds_start=3601, rng=None):
        """Pre-materialized numpy arrays, matching the evdata contract."""
        if rng is None:
            rng = np.random.default_rng(20)
        evdata = {}
        for c in [*comps, 'ndbd']:
            evdata[c] = {
                'u0':  rng.uniform(2545, 2695, n).astype(np.float32),
                'v0':  rng.uniform(1.0, 50.0, n).astype(np.float32),
                'ch1': rng.integers(1, n_ch + 1, n).astype(np.int32),
                'ch2': rng.integers(1, n_ch + 1, n).astype(np.int32),
                'ds':  np.full(n, ds_target, dtype=np.int32),
            }
        return evdata

    def _make_lss_params_full(self, n_ds=3, n_ch=10, rng=None):
        """lss_params dict with both sscl and dq fields populated."""
        if rng is None:
            rng = np.random.default_rng(30)
        sigma0_2d = rng.uniform(0.01, 0.1, (n_ds, n_ch))
        delta_p1_1r = rng.uniform(0.001, 0.01, (n_ds,))
        delta_p1_1l = -delta_p1_1r
        q_pars_arr      = rng.uniform(-1e-3, 1e-3, (n_ds, 3))
        q_pars_sample_1r = q_pars_arr + rng.uniform(1e-5, 1e-4, (n_ds, 3))
        q_pars_sample_1l = q_pars_arr - rng.uniform(1e-5, 1e-4, (n_ds, 3))
        return {
            'sigma0_2d':        sigma0_2d,
            'delta_p1_1r':      delta_p1_1r,
            'delta_p1_1l':      delta_p1_1l,
            'q_pars_arr':       q_pars_arr,
            'q_pars_sample_1r': q_pars_sample_1r,
            'q_pars_sample_1l': q_pars_sample_1l,
        }

    def test_returns_tuple_of_df_and_signal(self):
        comps = ['bkg_a', 'bkg_b']
        fit_cfg = self._make_fit_cfg()
        mcfit_table = self._make_mcfit_table(comps)
        evdata = self._make_evdata(comps)
        lss = self._make_lss_params_full()

        mctable, signal_df = make_lss_mctable_np(
            evdata, mcfit_table, fit_cfg, {'type': 'sscl', 'spt': 1},
            ds=3601, lss_params=lss, ds_start=3601,
            rng=np.random.default_rng(0),
        )
        assert isinstance(mctable, pd.DataFrame)
        assert isinstance(signal_df, pd.DataFrame)

    def test_mctable_has_mchist_and_integral(self):
        comps = ['bkg_a']
        fit_cfg = self._make_fit_cfg()
        mcfit_table = self._make_mcfit_table(comps)
        evdata = self._make_evdata(comps)
        lss = self._make_lss_params_full()

        mctable, _ = make_lss_mctable_np(
            evdata, mcfit_table, fit_cfg, {'type': 'sscl', 'spt': 1},
            ds=3601, lss_params=lss, ds_start=3601,
            rng=np.random.default_rng(0),
        )
        assert 'mchist' in mctable.columns
        assert 'integral' in mctable.columns
        assert mctable.index.tolist() == comps

    def test_signal_df_index_is_ndbd(self):
        comps = ['bkg_a']
        fit_cfg = self._make_fit_cfg()
        mcfit_table = self._make_mcfit_table(comps)
        evdata = self._make_evdata(comps)
        lss = self._make_lss_params_full()

        _, signal_df = make_lss_mctable_np(
            evdata, mcfit_table, fit_cfg, {'type': 'sscl', 'spt': 1},
            ds=3601, lss_params=lss, ds_start=3601,
            rng=np.random.default_rng(0),
        )
        assert list(signal_df.index) == ['ndbd']

    def test_unknown_syst_type_raises(self):
        comps = ['bkg_a']
        fit_cfg = self._make_fit_cfg()
        mcfit_table = self._make_mcfit_table(comps)
        evdata = self._make_evdata(comps)
        lss = self._make_lss_params_full()

        with pytest.raises(ValueError, match="Unknown syst type"):
            make_lss_mctable_np(
                evdata, mcfit_table, fit_cfg, {'type': 'bad'},
                ds=3601, lss_params=lss, ds_start=3601,
            )

    def test_integral_matches_hist_sum_sscl(self):
        comps = ['bkg_a', 'bkg_b']
        fit_cfg = self._make_fit_cfg()
        mcfit_table = self._make_mcfit_table(comps)
        evdata = self._make_evdata(comps)
        lss = self._make_lss_params_full()

        mctable, signal_df = make_lss_mctable_np(
            evdata, mcfit_table, fit_cfg, {'type': 'sscl', 'spt': 1},
            ds=3601, lss_params=lss, ds_start=3601,
            rng=np.random.default_rng(0),
        )
        for comp in comps:
            assert mctable.at[comp, 'integral'] == pytest.approx(
                mctable.at[comp, 'mchist'].sum()
            )
        assert signal_df.at['ndbd', 'integral'] == pytest.approx(
            signal_df.at['ndbd', 'mchist'].sum()
        )

    def test_integral_matches_hist_sum_dq(self):
        comps = ['bkg_a']
        fit_cfg = self._make_fit_cfg()
        mcfit_table = self._make_mcfit_table(comps)
        evdata = self._make_evdata(comps)
        lss = self._make_lss_params_full()

        mctable, signal_df = make_lss_mctable_np(
            evdata, mcfit_table, fit_cfg, {'type': 'dq', 'ptn': '1r'},
            ds=3601, lss_params=lss, ds_start=3601,
        )
        assert mctable.at['bkg_a', 'integral'] == pytest.approx(
            mctable.at['bkg_a', 'mchist'].sum()
        )
        assert signal_df.at['ndbd', 'integral'] == pytest.approx(
            signal_df.at['ndbd', 'mchist'].sum()
        )

    def test_common_filter_applied(self):
        """Events failing the common cut must be absent from the histogram.

        The old make_lss_mctable sscl path skipped the common cut; the new
        make_lss_mctable_np must apply it.
        """
        comps = ['bkg_a']
        # Wide hm so all test events pass ecuts; common cut removes everything
        hm = _make_hm(umin=2400, umax=2800, vmin=0, vmax=200,
                      binsize=10.0, threshold=-999.0)
        fit_cfg_cut = {'h2uv_model': hm, 'cuts': {'common': 'u > 10000'}}
        fit_cfg_nocut = {'h2uv_model': hm, 'cuts': {}}

        mcfit_table = self._make_mcfit_table(comps)
        evdata = self._make_evdata(comps, n=200)
        lss = self._make_lss_params_full()
        scfg = {'type': 'sscl', 'spt': 1}

        mctable_cut, _ = make_lss_mctable_np(
            evdata, mcfit_table, fit_cfg_cut, scfg,
            ds=3601, lss_params=lss, ds_start=3601,
            rng=np.random.default_rng(0),
        )
        mctable_nocut, _ = make_lss_mctable_np(
            evdata, mcfit_table, fit_cfg_nocut, scfg,
            ds=3601, lss_params=lss, ds_start=3601,
            rng=np.random.default_rng(0),
        )
        assert mctable_nocut.at['bkg_a', 'integral'] > 0, "events should exist without cut"
        assert mctable_cut.at['bkg_a', 'integral'] == pytest.approx(0.0)
