# scripts/ — parameter grids & SLURM submission

See `scripts/README.md` for the full command reference; this is the quick agent
orientation. Root `AGENTS.md` has the pipeline overview.

## Flow

```
gen_params.py <type>  ─►  scripts/params/<type>_params.tsv
                              │
        submit.sh / submit_packed.sh / submit_syst.sh
                              │
                        scripts/slurm/*.slurm  (sets PARAMS)
                              │
                        scripts/run_job.sh     (reads row by $SLURM_ARRAY_TASK_ID)
```

## Key files

- `gen_params.py` — argparse subparsers:
  `staging | input | fit | toys | sens | bias | syst_lss | syst_mc`.
  Writes `scripts/params/<type>_params.tsv`.
- `gen_syst_params.py` — systematics grid (`syst_lss` / `syst_mc`).
- `submit.sh` — regenerate TSV then `sbatch` for array types
  `fit | toys | sens | bias`; pass `--no-regen` to reuse an existing TSV.
- `submit_packed.sh` — packed (single-node, GNU `parallel`) mode for
  `staging | input | fit | toys`.
- `submit_syst.sh` — `lss | mc` systematics.
- `run_job.sh` — task runner; executes row `SLURM_ARRAY_TASK_ID` from `$PARAMS`,
  activates `/global/common/software/cuore/sachi/qpym2`, and sets
  `PYTENSOR_FLAGS` per task.
- `collect_results.py` — merge per-task toy CSVs.
- `params/` — generated TSVs (committed; regenerate with `gen_params.py`).
- `slurm/` — one `.slurm` per stage. Note `syst_mc_old.slurm`,
  `staging_retry.slurm` are special-purpose.
- `tests/` — `compare_staging.py`, `compare_fits.py`, `validate_inputs.py`.

## Gotchas

- `run_job.sh` reads the **last** TSV column as the command. Keep row format
  consistent when editing generators.
- `run_threshold250_bias.py` and `run_params_seq.py` are one-off/legacy — check
  before reusing.
- Never run `submit*.sh` from an agent without explicit user request; they
  enqueue real Perlmutter jobs.
