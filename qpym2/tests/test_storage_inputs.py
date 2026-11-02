"""
Pure unit tests for storage input generation routines (qpym2/storage/inputs.py).

Does not require external NERSC datasets or PyROOT installations.
"""

from pathlib import Path
import numpy as np
import pytest

from qpym2.storage.inputs import (
    build_shift_defs,
    build_source_hash_params,
    generate_source_template,
    generate_observed_data,
    create_storage_inputs_for_point,
)
from qpym2.storage.hasher import compute_input_uuid
from qpym2.storage.backend import read_template, read_observed_data
from qpym2.utils.hists import HistType
from qpym2.cfg.hydra_adapter_new import HydraConfigManager


def test_build_shift_defs():
    """Verify branch definition generation for RDF."""
    # Zero shift
    defs_zero = build_shift_defs(0.0)
    assert defs_zero == [('u', 'u0'), ('v', 'v0')]

    # ushiftp shift
    defs_shift = build_shift_defs(0.8, shift_type='ushiftp')
    assert defs_shift[0][0] == 'u'
    assert str(0.8 / 1000.0) in defs_shift[0][1]
    assert defs_shift[1][0] == 'v'
    assert str(0.8 / 1000.0) in defs_shift[1][1]


def test_source_uuid_isolation():
    """Verify co_shift changes 60Co UUID but leaves other sources invariant."""
    staging_id = "v1"
    cuts = {'common': 'v>150'}

    # 1. 60Co with different shifts
    p_co_00 = build_source_hash_params('60Co', staging_id, cuts=cuts, transforms=None)
    p_co_08 = build_source_hash_params('60Co', staging_id, cuts=cuts, transforms={'shift_type': 'ushiftp', 'shift_val': 0.8})
    p_co_10 = build_source_hash_params('60Co', staging_id, cuts=cuts, transforms={'shift_type': 'ushiftp', 'shift_val': 1.0})

    uuid_co_00 = compute_input_uuid(p_co_00)
    uuid_co_08 = compute_input_uuid(p_co_08)
    uuid_co_10 = compute_input_uuid(p_co_10)

    assert uuid_co_00 != uuid_co_08
    assert uuid_co_08 != uuid_co_10

    # 2. Invariant sources (208Tl, 238U, ndbd) have no shift
    p_tl = build_source_hash_params('208Tl', staging_id, cuts=cuts, transforms=None)
    p_u = build_source_hash_params('238U', staging_id, cuts=cuts, transforms=None)
    p_ndbd = build_source_hash_params('ndbd', staging_id, cuts=cuts, transforms=None)

    uuid_tl = compute_input_uuid(p_tl)
    uuid_u = compute_input_uuid(p_u)
    uuid_ndbd = compute_input_uuid(p_ndbd)

    # Identical repeated calls
    assert compute_input_uuid(p_tl) == uuid_tl
    assert compute_input_uuid(p_u) == uuid_u
    assert compute_input_uuid(p_ndbd) == uuid_ndbd


def test_generate_source_template_roundtrip_and_cache(tmp_path):
    """Test generating, persisting, and cache-detecting a source template."""
    input_dir = tmp_path / "fit_data" / "roi_2400_2570" / "thr200" / "h2uv_b2.5" / "input"
    input_dir.mkdir(parents=True, exist_ok=True)

    dummy_hist = np.array([[10.0, 20.0], [30.0, 40.0]])
    dummy_edges = [np.array([0.0, 1.0, 2.0]), np.array([0.0, 1.0, 2.0])]

    class DummyHM:
        hist_type = HistType.h2uv
        umin = 1697.0
        umax = 1817.0
        threshold = 200.0
        skip_roi = None

    def mock_read_hist(fpath, hm, treename, defs, filters, rtype):
        return (dummy_hist, *dummy_edges)

    cuts = {'co_cut': 'v>150'}
    transforms = {'shift_type': 'ushiftp', 'shift_val': 0.8}

    # First call: writes template
    path1, uuid1, cached1 = generate_source_template(
        source_name='60Co',
        staging_fpath=tmp_path / 'staging' / '60Co.root',
        treename='uvtree',
        hm=DummyHM(),
        cuts_dict=cuts,
        transforms_dict=transforms,
        input_dir=input_dir,
        read_hist_fn=mock_read_hist,
    )
    assert not cached1
    assert path1.exists()
    assert f"60Co_{uuid1}.h5" == path1.name

    # Verify content
    read_h, read_e, meta = read_template(path1)
    np.testing.assert_array_equal(read_h, dummy_hist)
    assert meta['component'] == '60Co'
    assert meta['transforms']['shift_val'] == 0.8

    # Second call: must be detected as CACHE HIT
    path2, uuid2, cached2 = generate_source_template(
        source_name='60Co',
        staging_fpath=tmp_path / 'staging' / '60Co.root',
        treename='uvtree',
        hm=DummyHM(),
        cuts_dict=cuts,
        transforms_dict=transforms,
        input_dir=input_dir,
        read_hist_fn=mock_read_hist,
    )
    assert cached2
    assert path1 == path2
    assert uuid1 == uuid2


def test_generate_observed_data_roundtrip_and_cache(tmp_path):
    """Test generating and caching observed data."""
    input_dir = tmp_path / "input"
    input_dir.mkdir(parents=True, exist_ok=True)

    dummy_events = np.array([[1700.0, 50.0], [1750.0, 60.0]])

    class DummyHM:
        hist_type = HistType.h2uv
        umin = 1697.0
        umax = 1817.0
        threshold = 200.0
        skip_roi = None

    def mock_read_evlist(fpath, treename, defs, filters, rtype):
        return dummy_events

    cuts = {'common': 'v>150'}

    # 1. Unblinded data
    p1, c1 = generate_observed_data(
        data_fpath=tmp_path / 'data.root',
        treename='uvtree',
        hm=DummyHM(),
        cuts_dict=cuts,
        blinded=False,
        input_dir=input_dir,
        read_evlist_fn=mock_read_evlist,
    )
    assert not c1
    assert p1.name == 'data.h5'

    # Cache hit
    p2, c2 = generate_observed_data(
        data_fpath=tmp_path / 'data.root',
        treename='uvtree',
        hm=DummyHM(),
        cuts_dict=cuts,
        blinded=False,
        input_dir=input_dir,
        read_evlist_fn=mock_read_evlist,
    )
    assert c2
    assert p1 == p2

    # Read back
    data_arr, _, meta = read_observed_data(p1)
    np.testing.assert_array_equal(data_arr, dummy_events)
    assert not meta['blinded']


def test_create_storage_inputs_dry_run():
    """Verify dry_run previews paths and UUIDs using real config."""
    cfg = HydraConfigManager.get_config(environment='local', analysis='base')

    manifest = create_storage_inputs_for_point(
        cfg=cfg,
        model_key='tech_note_200',
        co_cut='gamma10_v150',
        ndbd_cutoff='u1740',
        co_shift=0.8,
        blinded=False,
        dry_run=True,
    )

    assert manifest['coords']['roi'] == (2400.0, 2570.0)
    assert manifest['coords']['threshold'] == 200.0
    assert 'roi_2400_2570/thr200/h2uv_b2.5/input' in manifest['input_dir']
    assert 'templates' in manifest
    assert len(manifest['templates']) > 0
    assert 'data' in manifest
    assert manifest['data']['path'].endswith('data.h5')
