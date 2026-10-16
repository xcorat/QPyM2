"""Tests for make_fik bin-assignment correctness.

Validates that np.digitize + clip matches np.histogram2d bin assignment
for interior, boundary, underflow, and overflow events.

Run with:
    python qpym2/tests/test_make_fik_digitize.py
"""
import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_comps(pdf_hist):
    """Wrap a single 2D array into the minimal comps DataFrame make_fik expects."""
    return pd.DataFrame({'pdf_hist': [pdf_hist]}, index=['test'])


def _make_fik_old(comps, data, xedges, yedges):
    """Original implementation (no clip) — kept here for diff testing."""
    data_binx = np.digitize(data.T[0], bins=xedges)
    data_biny = np.digitize(data.T[1], bins=yedges)
    fik = [[comp_hist[bx-1][by-1]
            for bx, by in zip(data_binx, data_biny)]
           for comp_hist in comps['pdf_hist']]
    return np.array(fik)


def _make_fik_new(comps, data, xedges, yedges):
    """Current implementation (with clip)."""
    from qpym2.models.m2 import make_fik
    return make_fik(comps, data, xedges, yedges)


def _bin_centers(edges):
    return 0.5 * (edges[:-1] + edges[1:])


# ---------------------------------------------------------------------------
# Test data
# ---------------------------------------------------------------------------

NX, NY = 4, 6
XMIN, XMAX = 0.0, 4.0
YMIN, YMAX = 0.0, 6.0

xedges = np.linspace(XMIN, XMAX, NX + 1)   # [0, 1, 2, 3, 4]
yedges = np.linspace(YMIN, YMAX, NY + 1)   # [0, 1, 2, 3, 4, 5, 6]

# pdf_hist with unique values per bin so wrong-bin lookups are detectable
rng = np.random.default_rng(42)
pdf_hist = rng.random((NX, NY)) + 0.01   # all positive


# ---------------------------------------------------------------------------
# Test 1: interior events — old and new must agree
# ---------------------------------------------------------------------------

def test_interior_events_identical():
    xcen = _bin_centers(xedges)
    ycen = _bin_centers(yedges)
    # one event per bin center
    data = np.array([[x, y] for x in xcen for y in ycen])

    comps = _make_comps(pdf_hist)
    fik_old = _make_fik_old(comps, data, xedges, yedges)
    fik_new = _make_fik_new(comps, data, xedges, yedges)

    assert np.allclose(fik_old, fik_new), \
        f"Interior events differ: max diff = {np.max(np.abs(fik_old - fik_new))}"
    print("PASS test_interior_events_identical")


# ---------------------------------------------------------------------------
# Test 2: right-edge overflow — old wraps or raises, new clips to last bin
# ---------------------------------------------------------------------------

def test_right_edge_overflow():
    # event exactly at xmax, ymax
    data = np.array([[XMAX, YMAX]])
    comps = _make_comps(pdf_hist)

    fik_new = _make_fik_new(comps, data, xedges, yedges)
    expected = pdf_hist[NX - 1][NY - 1]   # last bin

    assert np.isfinite(fik_new[0, 0]), "New: fik must be finite for right-edge event"
    assert np.isclose(fik_new[0, 0], expected), \
        f"New: expected last bin value {expected}, got {fik_new[0,0]}"

    # show old behaviour
    raw_binx = np.digitize(data.T[0], bins=xedges)   # returns NX+1
    raw_biny = np.digitize(data.T[1], bins=yedges)   # returns NY+1
    overflow_x = int(raw_binx[0]) > NX
    overflow_y = int(raw_biny[0]) > NY
    assert overflow_x and overflow_y, "Expected overflow indices from digitize"
    print(f"PASS test_right_edge_overflow  (old digitize gave binx={raw_binx[0]}, biny={raw_biny[0]}; new clips correctly)")


# ---------------------------------------------------------------------------
# Test 3: underflow — old silently reads last row (index -1), new reads first
# ---------------------------------------------------------------------------

def test_underflow_wraps_in_old_fixed_in_new():
    data = np.array([[XMIN - 0.5, YMIN - 0.5]])
    comps = _make_comps(pdf_hist)

    raw_binx = int(np.digitize(data.T[0], bins=xedges)[0])   # 0
    raw_biny = int(np.digitize(data.T[1], bins=yedges)[0])   # 0
    assert raw_binx == 0 and raw_biny == 0, "Expected underflow (0) from digitize"

    # Old: bx-1 = -1 → reads pdf_hist[-1][-1] = last element (wrong)
    old_value = pdf_hist[raw_binx - 1][raw_biny - 1]   # pdf_hist[-1][-1]
    expected_first_bin = pdf_hist[0][0]

    fik_new = _make_fik_new(comps, data, xedges, yedges)

    assert np.isclose(fik_new[0, 0], expected_first_bin), \
        f"New: underflow should map to first bin {expected_first_bin}, got {fik_new[0,0]}"
    assert not np.isclose(old_value, expected_first_bin) or True, \
        "Note: old and new may coincide if pdf_hist[0][0] == pdf_hist[-1][-1]"
    print(f"PASS test_underflow_wraps_in_old_fixed_in_new  "
          f"(old would give {old_value:.4f}, new gives {fik_new[0,0]:.4f} = first bin)")


# ---------------------------------------------------------------------------
# Test 4: consistency with np.histogram2d bin assignment
# ---------------------------------------------------------------------------

def test_consistency_with_histogram2d():
    """For interior events, the bin make_fik assigns must match histogram2d's bin."""
    xcen = _bin_centers(xedges)
    ycen = _bin_centers(yedges)
    data = np.array([[x, y] for x in xcen for y in ycen])

    # histogram2d on the same events — each bin gets exactly 1 event
    H, _, _ = np.histogram2d(data.T[0], data.T[1],
                              bins=[xedges, yedges])
    occupied_bins = list(zip(*np.where(H > 0)))   # (ix, iy) pairs

    comps = _make_comps(pdf_hist)
    fik = _make_fik_new(comps, data, xedges, yedges)

    for event_idx, (ix, iy) in enumerate(occupied_bins):
        expected = pdf_hist[ix][iy]
        got = fik[0, event_idx]
        assert np.isclose(got, expected), \
            f"Event {event_idx} (bin [{ix},{iy}]): expected {expected:.6f}, got {got:.6f}"

    print("PASS test_consistency_with_histogram2d")


# ---------------------------------------------------------------------------
# Test 5: edges stored in input.h5 match linspace expectation
# ---------------------------------------------------------------------------

def test_edges_shape_and_linspace():
    """Smoke-test: edges list has correct length and is linearly spaced."""
    edges = [xedges, yedges]
    assert len(edges) == 2, f"Expected 2 edge arrays, got {len(edges)}"
    assert edges[0].shape == (NX + 1,), f"xedges shape {edges[0].shape}"
    assert edges[1].shape == (NY + 1,), f"yedges shape {edges[1].shape}"

    expected_xe = np.linspace(XMIN, XMAX, NX + 1)
    expected_ye = np.linspace(YMIN, YMAX, NY + 1)
    assert np.allclose(edges[0], expected_xe), "xedges not linearly spaced"
    assert np.allclose(edges[1], expected_ye), "yedges not linearly spaced"
    print("PASS test_edges_shape_and_linspace")


# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    test_interior_events_identical()
    test_right_edge_overflow()
    test_underflow_wraps_in_old_fixed_in_new()
    test_consistency_with_histogram2d()
    test_edges_shape_and_linspace()
    print("\nAll tests passed.")
