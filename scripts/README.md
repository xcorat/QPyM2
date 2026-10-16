# scripts/

Utilities for generating SLURM job parameter grids and submitting batch jobs on Perlmutter.

## Files

| File | Purpose |
|------|---------|
| `gen_params.py` | Unified parameter-grid generator — writes TSV files consumed by SLURM array jobs |
| `submit.sh` | Regenerates the param grid then submits the corresponding SLURM array |
| `run_job.sh` | Generic runner called inside each SLURM task; reads its command from the TSV by `SLURM_ARRAY_TASK_ID` |

## Directories

| Directory | Purpose |
|-----------|---------|
| `params/` | Generated TSV parameter files (one per run type) |
| `slurm/` | SLURM job scripts (`staging.slurm`, `input.slurm`, `fit.slurm`, `toys.slurm`, …) |
| `tests/` | Validation and smoke-test scripts (`compare_staging.py`, `validate_inputs.py`) |

---

## Workflow

### 1. Generate a param grid and submit in one step

```bash
# From the project root (QPyM2/)
bash scripts/submit.sh <run_type> [gen_params flags] [-- sbatch flags]
```

`run_type` is one of: `staging`, `input`, `fit`, `toys`, `sens`, `bias`.

Examples:

```bash
bash scripts/submit.sh staging                                 # default env
bash scripts/submit.sh staging --environment tech_note         # explicit env
bash scripts/submit.sh staging --test                          # 2-task debug queue smoke test
bash scripts/submit.sh fit --chains 4 --draws 10000
bash scripts/submit.sh toys --n-fits 200
bash scripts/submit.sh sens --n-fits 100
bash scripts/submit.sh bias --n-fits 100 -- -q debug         # pass extra sbatch flags
```

### 2. Generate a param grid only

```bash
python scripts/gen_params.py <subcommand> [flags]
```

Subcommands mirror the run types.  The output TSV is written to `scripts/params/<run_type>_params.tsv`.

```bash
python scripts/gen_params.py staging
python scripts/gen_params.py staging --environment tech_note --analysis tech_note_200
python scripts/gen_params.py fit -i                            # interactive grid selection
python scripts/gen_params.py fit --co-shift 0.8 1.0 1.2
python scripts/gen_params.py toys --n-fits 200
python scripts/gen_params.py sens --n-fits 100
python scripts/gen_params.py bias --n-fits 100
python scripts/gen_params.py input --co-shift 0.9 1.0 1.1
```

### 3. Compare staging outputs

```bash
python scripts/tests/compare_staging.py <folder_ref> <folder_new> [options]

# Show up to 20 differing events for TeO2 files:
python scripts/tests/compare_staging.py \
    /path/to/staging_ref \
    /path/to/staging_new \
    --component TeO2 --show-diff 20
```

---

## Running Fits

### Test a single fit locally

```bash
python qpym2/run/run_single_fit.py \
    --environment m2v5-1 --analysis tech_note \
    --co-cut None --ndbd-cutoff None --co-shift 0.0 \
    --chains 8 --draws 20000 --prior-width 2syst
```

### Generate the fit parameter grid

```bash
python scripts/gen_params.py fit \
    --environment m2v5-1 --analysis tech_note \
    --chains 8 --draws 20000 --prior-width 2syst
```

Output: `scripts/params/fit_params.tsv` (63 tasks: all combinations of co-cut × ndbd-cutoff × co-shift).

### Submit to debug queue (2 jobs, 30 min)

```bash
sbatch --array=1-2%2 --qos=debug --time=00:30:00 scripts/slurm/fit.slurm
```

### Submit full fit grid

```bash
sbatch --array=1-63%20 scripts/slurm/fit.slurm
```

---

## Running Sensitivity and Bias Tests

Sens/bias jobs use `run_toy_batch.py`, which **loads the already-completed fit** from disk and
re-uses the `comps` DataFrame (which already carries the prior mu/sigma from that fit).
**There is no `--prior-width` argument** — the prior width is inherited from the parent fit automatically.

### Correct environment/analysis flags

Always match `--environment` and `--analysis` to the fit grid:

```bash
# Generate TSVs (42 tasks: [None, gamma10_v150] × NDBD_CUTS × CO_SHIFTS)
# Note: v150 co-cut is excluded — only None and gamma10_v150 are needed
python scripts/gen_params.py sens --environment m2v5-1 --analysis tech_note --n-fits 200 --co-cut None gamma10_v150
python scripts/gen_params.py bias --environment m2v5-1 --analysis tech_note --n-fits 200 --co-cut None gamma10_v150
```

### Submit to debug queue (2 tasks, 30 min)

```bash
bash scripts/submit.sh sens --environment m2v5-1 --analysis tech_note --n-fits 200 --co-cut None gamma10_v150 -- --qos=debug --time=00:30:00
bash scripts/submit.sh bias --environment m2v5-1 --analysis tech_note --n-fits 200 --co-cut None gamma10_v150 -- --qos=debug --time=00:30:00
```

### Submit full sens/bias grid

```bash
bash scripts/submit.sh sens --environment m2v5-1 --analysis tech_note --n-fits 200 --co-cut None gamma10_v150
bash scripts/submit.sh bias --environment m2v5-1 --analysis tech_note --n-fits 200 --co-cut None gamma10_v150
```

---

## SLURM scripts (`slurm/`)

| Script | Description |
|--------|-------------|
| `staging.slurm` | Single-node ROOT → staging ROOT |
| `staging_parallel_all.slurm` | Parallel staging across all components |
| `staging_retry.slurm` | Re-run failed staging tasks |
| `input.slurm` | Staging ROOT → fit input HDF5 |
| `fit.slurm` | Main Bayesian fits |
| `toys.slurm` | Toy batch fits |
| `sens.slurm` | Sensitivity batch fits |
| `bias.slurm` | Bias batch fits |

Each SLURM script sets a `PARAMS` environment variable pointing to its TSV and delegates execution to `run_job.sh`.
