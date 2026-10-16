"""Unit tests for qpym2.io.staging_transforms — pure-numpy transform functions.

Tests are pure numpy (no ROOT), run with:
    pytest qpym2/tests/test_staging_transforms.py
"""

import numpy as np
import pytest

from qpym2.io.staging_transforms import (
    apply_residual_sampling_ds,
    apply_sigma_scaling_ds,
)

_S2 = np.sqrt(2)


def _make_sigma_params(n_ds=3, n_ch=10, rng=None):
    if rng is None:
        rng = np.random.default_rng(0)
    sigma0_2d = rng.uniform(0.01, 0.1, size=(n_ds, n_ch))
    delta_p1 = rng.uniform(0.001, 0.01, size=(n_ds,))
    return sigma0_2d, delta_p1


# ---------------------------------------------------------------------------
# apply_sigma_scaling_ds
# ---------------------------------------------------------------------------

class TestApplySigmaScalingDs:
    def test_non_target_ds_unchanged(self):
        """Events from ds != target_ds must be returned unchanged."""
        rng = np.random.default_rng(1)
        n = 100
        u0 = rng.uniform(2540, 2700, n)
        v0 = rng.uniform(0.1, 50, n)
        ch1 = rng.integers(1, 11, n)
        ch2 = rng.integers(1, 11, n)
        ds_arr = np.full(n, 3602, dtype=np.int32)
        sigma0_2d, delta_p1 = _make_sigma_params()

        u_out, v_out = apply_sigma_scaling_ds(
            u0, v0, ch1, ch2, ds_arr,
            sigma0_2d, delta_p1, target_ds=3601, ds_start=3601,
            rng=np.random.default_rng(0),
        )
        np.testing.assert_array_equal(u_out, u0)
        np.testing.assert_array_equal(v_out, v0)

    def test_target_ds_events_transformed(self):
        """Target-ds events should be stochastically shifted (distributional change)."""
        rng = np.random.default_rng(2)
        n = 500
        u0 = np.full(n, 2600.0)
        v0 = np.full(n, 10.0)
        ch1 = np.ones(n, dtype=np.int32)
        ch2 = np.ones(n, dtype=np.int32)
        ds_arr = np.full(n, 3601, dtype=np.int32)
        # Large delta_p1 so the smearing effect is detectable
        sigma0_2d = np.full((3, 10), 0.05)
        delta_p1 = np.array([0.1, 0.0, 0.0])

        u_out, v_out = apply_sigma_scaling_ds(
            u0, v0, ch1, ch2, ds_arr,
            sigma0_2d, delta_p1, target_ds=3601, ds_start=3601,
            rng=np.random.default_rng(42),
        )
        assert not np.all(u_out == u0), "target-ds events should be stochastically shifted"

    def test_empty_mask_returns_copy(self):
        """No events from target_ds → return independent copies of u0, v0."""
        rng = np.random.default_rng(3)
        n = 50
        u0 = rng.uniform(2540, 2700, n)
        v0 = rng.uniform(0.1, 50, n)
        ch1 = rng.integers(1, 11, n)
        ch2 = rng.integers(1, 11, n)
        ds_arr = np.full(n, 3603, dtype=np.int32)
        sigma0_2d, delta_p1 = _make_sigma_params()

        u_out, v_out = apply_sigma_scaling_ds(
            u0, v0, ch1, ch2, ds_arr,
            sigma0_2d, delta_p1, target_ds=3601, ds_start=3601,
        )
        np.testing.assert_array_equal(u_out, u0)
        np.testing.assert_array_equal(v_out, v0)
        assert u_out is not u0
        assert v_out is not v0

    def test_output_length_preserved(self):
        """Output arrays must have the same length as input."""
        rng = np.random.default_rng(4)
        n = 300
        u0 = rng.uniform(2540, 2700, n)
        v0 = rng.uniform(0.1, 50, n)
        ch1 = rng.integers(1, 11, n)
        ch2 = rng.integers(1, 11, n)
        ds_arr = rng.integers(3601, 3604, n)
        sigma0_2d, delta_p1 = _make_sigma_params()

        u_out, v_out = apply_sigma_scaling_ds(
            u0, v0, ch1, ch2, ds_arr,
            sigma0_2d, delta_p1, target_ds=3601, ds_start=3601,
        )
        assert len(u_out) == n
        assert len(v_out) == n

    def test_non_target_unaffected_in_mixed_array(self):
        """Mixed ds_arr: non-target events kept exactly, target events shifted."""
        rng = np.random.default_rng(5)
        n = 200
        u0 = rng.uniform(2540, 2700, n)
        v0 = rng.uniform(0.1, 50, n)
        ch1 = rng.integers(1, 11, n)
        ch2 = rng.integers(1, 11, n)
        # First half target ds, second half other ds
        ds_arr = np.where(np.arange(n) < n // 2, 3601, 3602).astype(np.int32)
        sigma0_2d = np.full((3, 10), 0.05)
        delta_p1 = np.array([0.1, 0.0, 0.0])

        u_out, v_out = apply_sigma_scaling_ds(
            u0, v0, ch1, ch2, ds_arr,
            sigma0_2d, delta_p1, target_ds=3601, ds_start=3601,
            rng=np.random.default_rng(0),
        )
        np.testing.assert_array_equal(u_out[n // 2:], u0[n // 2:])
        np.testing.assert_array_equal(v_out[n // 2:], v0[n // 2:])
        assert not np.all(u_out[:n // 2] == u0[:n // 2])


# ---------------------------------------------------------------------------
# apply_residual_sampling_ds
# ---------------------------------------------------------------------------

class TestApplyResidualSamplingDs:
    def _make_qpars(self, n_ds=3, n_par=3, rng=None):
        if rng is None:
            rng = np.random.default_rng(10)
        qpars = rng.uniform(-0.001, 0.001, (n_ds, n_par))
        qpars_s = qpars + rng.uniform(1e-4, 5e-4, (n_ds, n_par))
        return qpars, qpars_s

    def test_non_target_ds_unchanged(self):
        rng = np.random.default_rng(11)
        n = 100
        u0 = rng.uniform(2540, 2700, n)
        v0 = rng.uniform(0.1, 50, n)
        ds_arr = np.full(n, 3602, dtype=np.int32)
        qp, qps = self._make_qpars()

        u_out, v_out = apply_residual_sampling_ds(
            u0, v0, ds_arr, qp, qps, target_ds=3601, ds_start=3601,
        )
        np.testing.assert_array_equal(u_out, u0)
        np.testing.assert_array_equal(v_out, v0)

    def test_empty_mask_returns_copy(self):
        rng = np.random.default_rng(12)
        n = 50
        u0 = rng.uniform(2540, 2700, n)
        v0 = rng.uniform(0.1, 50, n)
        ds_arr = np.full(n, 3603, dtype=np.int32)
        qp, qps = self._make_qpars()

        u_out, v_out = apply_residual_sampling_ds(
            u0, v0, ds_arr, qp, qps, target_ds=3601, ds_start=3601,
        )
        np.testing.assert_array_equal(u_out, u0)
        assert u_out is not u0

    def test_polynomial_applied_correctly(self):
        """Spot-check: u' = u + a*√2 + b*u + c/√2*(u²+v²)."""
        n = 10
        u0 = np.full(n, 2600.0)
        v0 = np.full(n, 10.0)
        ds_arr = np.full(n, 3601, dtype=np.int32)

        # Controlled shift: a=0.01, b=0.002, c=0 → c term vanishes, easy to check
        qp  = np.array([[0.00, 0.000, 0.0]])
        qps = np.array([[0.01, 0.002, 0.0]])
        a, b, c = 0.01, 0.002, 0.0

        u_exp = u0 + a * _S2 + b * u0 + (c / _S2) * (u0**2 + v0**2)
        v_exp = v0           + b * v0  +  c * _S2  *  u0 * v0

        u_out, v_out = apply_residual_sampling_ds(
            u0, v0, ds_arr, qp, qps, target_ds=3601, ds_start=3601,
        )
        np.testing.assert_allclose(u_out, u_exp, rtol=1e-10)
        np.testing.assert_allclose(v_out, v_exp, rtol=1e-10)

    def test_non_target_unaffected_in_mixed_array(self):
        """Only target-ds events are shifted; others kept exact."""
        rng = np.random.default_rng(13)
        n = 200
        u0 = rng.uniform(2540, 2700, n)
        v0 = rng.uniform(0.1, 50, n)
        ds_arr = np.where(np.arange(n) < n // 2, 3601, 3602).astype(np.int32)

        qp  = np.zeros((2, 3))
        qps = np.array([[0.01, 0.002, 0.0],   # ds index 0 → large shift
                        [0.00, 0.000, 0.0]])

        u_out, v_out = apply_residual_sampling_ds(
            u0, v0, ds_arr, qp, qps, target_ds=3601, ds_start=3601,
        )
        np.testing.assert_array_equal(u_out[n // 2:], u0[n // 2:])
        assert not np.allclose(u_out[:n // 2], u0[:n // 2])
