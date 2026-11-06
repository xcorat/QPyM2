# notebooks/ — interactive notebooks and workbenches

Guidance for interactive development and data exploration in QPyM2.

## Layout

- `notebooks/marimo/` — Modern interactive Marimo workbenches (versioned in git):
  - `staging_inventory.py` — Staging tier and input dataset inventory, config inspection, ROOT tree entries, and JAGS table verification.
  - `input_creation.py` — Storage input generation and validation workbench for v0.6: preview deterministic UUIDs, check cache hit/miss status, trigger source-level template generation, inspect embedded HDF5 metadata, and plot 2D distributions and 1D projections.
  - `fit_run_results.py` — Interactive fit result exploration and posterior trace inspection.
  - `bias_study_workbench.py` — Signal injection bias study inspection, calibration curve fitting, parameter stability checks, raw vs bias-corrected posterior comparison, and systematic marginalization.
  - `bmmc_study_workbench.py` — Background Model Monte Carlo (BMMC) systematics investigation & diagnostics workbench: inspect stat-only unperturbed JAGS chain draws, fit likelihoods, weight collapse & effective sample size ($N_{\text{eff}}$), read-only bias calibration resolution, composite posterior reconstruction ($\Delta_{\text{mode}}$, $\Delta_{\text{left}}$, $\Delta_{\text{right}}$), sample-size convergence, and background component sensitivity correlations.
  - `lss_study_workbench.py` — Lineshape Systematics (LSS: `lss_sscl` and `lss_dq`) interactive workbench: inspect per-dataset MCMC fits, bias-corrected mode shifts ($\Delta_{ds}$), signed quadrature sum envelopes ($\delta_{\text{right}}$, $\delta_{\text{left}}$), diagnostic per-dataset barplots, dominance & Pareto breakdown, posterior overlays, and Format v2 artifact inspection.
- `notebooks/old/` — Legacy Jupyter (`.ipynb`) scratch notebooks (gitignored, local filesystem only).

## Marimo Guidelines

- Notebooks are written as pure Python scripts driven by the reactive `marimo` runtime (`app = marimo.App(...)`).
- Follow Marimo DAG rules:
  - **No cycle dependencies**: Cell dataflow must form a directed acyclic graph.
  - **Separation of UI definition and value access**: Never access a UI element's `.value` in the cell that created it; read `.value` in downstream cells.
  - **Single owning cell per public variable**: Local helpers inside cells should start with `_` to avoid unintentional global bindings.
- Test notebooks before committing:
  ```bash
  ruff check notebooks/marimo/
  ruff format --check notebooks/marimo/
  marimo export html notebooks/marimo/<notebook>.py -o /tmp/test.html
  ```
