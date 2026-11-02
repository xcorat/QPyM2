"""
Pure unit tests for qpym2.studies module (config, points, assembly, study).

Does not require external NERSC datasets, live ROOT data, or PyROOT.
"""

import pytest
from omegaconf import OmegaConf

from qpym2.studies.config import (
    resolve_study_dir,
    find_study_config_file,
    load_study_config,
)
from qpym2.studies.point import StudyPoint
from qpym2.studies.assembly import expand_study_points
from qpym2.studies.study import Study


def test_resolve_study_dir(tmp_path):
    """Verify study directory resolution and error handling."""
    study_folder = tmp_path / "my_custom_study"
    study_folder.mkdir()

    # Direct path
    resolved = resolve_study_dir(study_folder)
    assert resolved == study_folder.resolve()

    # Non-existent path raises FileNotFoundError
    with pytest.raises(FileNotFoundError):
        resolve_study_dir(tmp_path / "does_not_exist")


def test_find_study_config_file(tmp_path):
    """Verify primary study configuration file discovery."""
    study_dir = tmp_path / "test_study"
    study_dir.mkdir()

    with pytest.raises(FileNotFoundError):
        find_study_config_file(study_dir)

    cfg_file = study_dir / "study.yaml"
    cfg_file.write_text("study:\n  name: test\n")
    assert find_study_config_file(study_dir) == cfg_file


def test_load_study_config_and_overrides(tmp_path):
    """Verify loading study config with defaults inheritance and OmegaConf merging."""
    study_dir = tmp_path / "shift_study"
    study_dir.mkdir()

    study_yaml = study_dir / "study.yaml"
    study_yaml.write_text(
        """
study:
  name: shift_study
  description: "Test study description"

defaults:
  environment: local
  staging: default
  analysis: base
  models: tech_note
  model_key: tech_note_200

reference:
  co_cut: gamma10_v150
  co_shift: 0.0

sweep:
  co_shift: [0.0, 0.8, 1.2]
"""
    )

    cfg, resolved_dir = load_study_config(study_dir)
    assert resolved_dir == study_dir.resolve()
    assert cfg.study.name == "shift_study"
    assert cfg.study.description == "Test study description"
    assert cfg.reference.co_shift == 0.0
    assert list(cfg.sweep.co_shift) == [0.0, 0.8, 1.2]

    # Verify base config was loaded and retained
    assert cfg.infrastructure.data_root is not None
    assert cfg.staging.version is not None


def test_study_config_auxiliary_yaml_overrides(tmp_path):
    """Verify auxiliary YAML files in study directory override base settings."""
    study_dir = tmp_path / "aux_study"
    study_dir.mkdir()

    (study_dir / "study.yaml").write_text(
        """
study:
  name: aux_study
defaults:
  environment: local
"""
    )

    # Auxiliary override file
    (study_dir / "cuts_override.yaml").write_text(
        """
reference:
  custom_cut_key: "custom_value"
"""
    )

    cfg, _ = load_study_config(study_dir)
    assert cfg.reference.custom_cut_key == "custom_value"


def test_expand_study_points_nominal():
    """Verify single baseline point expansion when sweep is empty."""
    ref = {
        "co_cut": "gamma10_v150",
        "ndbd_cutoff": "u1740",
        "co_shift": 0.5,
        "blinded": False,
    }

    points = expand_study_points(reference=ref, sweep=None)
    assert len(points) == 1
    p = points[0]
    assert p.point_id == "nominal"
    assert p.task_id == 0
    assert p.co_cut == "gamma10_v150"
    assert p.co_shift == 0.5
    assert p.blinded is False


def test_expand_study_points_cartesian():
    """Verify multi-axis Cartesian sweep expansion."""
    ref = {
        "co_cut": "gamma10_v150",
        "ndbd_cutoff": "u1740",
        "co_shift": 0.0,
        "blinded": False,
    }
    sweep = {
        "co_shift": [0.0, 0.8],
        "co_cuts": ["gamma10_v150", "v150"],
    }

    points = expand_study_points(reference=ref, sweep=sweep)
    assert len(points) == 4

    p_ids = [p.point_id for p in points]
    assert len(set(p_ids)) == 4

    # Verify task_id assignment
    for i, p in enumerate(points):
        assert p.task_id == i

    # Verify parameter combinations
    combos = [(p.co_shift, p.co_cut) for p in points]
    expected = [
        (0.0, "gamma10_v150"),
        (0.0, "v150"),
        (0.8, "gamma10_v150"),
        (0.8, "v150"),
    ]
    assert combos == expected


def test_study_point_check_inputs(monkeypatch):
    """Verify StudyPoint.check_inputs calculates readiness metrics properly."""
    point = StudyPoint(
        point_id="p0_shift0.8",
        task_id=0,
        co_shift=0.8,
    )

    mock_manifest = {
        "input_dir": "/mock/input",
        "param_dir": "/mock/param",
        "data": {"path": "/mock/input/data.h5", "cached": True},
        "templates": {
            "sourceA": {"path": "/mock/input/sourceA_123.h5", "cached": True},
            "sourceB": {"path": "/mock/input/sourceB_456.h5", "cached": False},
        },
    }

    monkeypatch.setattr(
        "qpym2.studies.point.create_storage_inputs_for_point",
        lambda **kwargs: mock_manifest,
    )

    dummy_cfg = OmegaConf.create({})
    report = point.check_inputs(dummy_cfg)

    assert report["point_id"] == "p0_shift0.8"
    assert report["total_templates"] == 2
    assert report["cached_templates"] == 1
    assert report["missing_templates"] == 1
    assert report["missing_template_names"] == ["sourceB"]
    assert report["missing_template_paths"] == ["/mock/input/sourceB_456.h5"]
    assert report["data_cached"] is True
    assert report["is_ready"] is False
    assert report["status"] == "MISSING_INPUTS"


def test_study_orchestrator(tmp_path, monkeypatch):
    """Verify Study check_inputs, missing inputs tracking, and manifest serialization."""
    study_dir = tmp_path / "study_run"
    study_dir.mkdir()

    points = [
        StudyPoint(point_id="p0", task_id=0, co_shift=0.0),
        StudyPoint(point_id="p1", task_id=1, co_shift=0.8),
    ]

    mock_manifests = {
        0.0: {
            "input_dir": "/mock/input",
            "param_dir": "/mock/param",
            "data": {"path": "/mock/input/data.h5", "cached": True},
            "templates": {
                "s1": {"path": "/mock/input/s1.h5", "cached": True},
            },
        },
        0.8: {
            "input_dir": "/mock/input",
            "param_dir": "/mock/param",
            "data": {"path": "/mock/input/data.h5", "cached": False},
            "templates": {
                "s1": {"path": "/mock/input/s1.h5", "cached": False},
            },
        },
    }

    monkeypatch.setattr(
        "qpym2.studies.point.create_storage_inputs_for_point",
        lambda **kwargs: mock_manifests[kwargs["co_shift"]],
    )

    study = Study(
        name="test_orchestrator",
        cfg=OmegaConf.create({}),
        points=points,
        study_dir=study_dir,
        description="Test orchestrator unit test",
    )

    assert len(study) == 2

    # Test check_inputs DataFrame
    df = study.check_inputs()
    assert len(df) == 2
    assert bool(df.loc[0, "is_ready"]) is True
    assert bool(df.loc[1, "is_ready"]) is False

    # Test missing inputs dictionary
    missing = study.get_missing_inputs()
    assert "p0" not in missing
    assert "p1" in missing
    assert "/mock/input/s1.h5" in missing["p1"]
    assert "/mock/input/data.h5" in missing["p1"]

    # Test manifest export
    m_dict = study.to_manifest()
    assert m_dict["study_name"] == "test_orchestrator"
    assert m_dict["n_points"] == 2

    # Test save manifest file
    m_file = study.save_manifest()
    assert m_file.is_file()
    assert m_file == study_dir / "test_orchestrator.manifest"
