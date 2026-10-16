"""Single source of truth for parameter grid axes.

``ALL_*`` constants are the complete, index-stable lists of every possible
value for each axis.  Run scripts that need a subset should slice directly
into these (e.g. ``ALL_CO_CUTS[3:4]``), so the slice indices remain stable
and the full universe of options is always visible here.

The non-prefixed names define the **default batch-job grid** used by
``scripts/gen_job_params.py`` and similar launchers.
"""

# ---------------------------------------------------------------------------
# Complete ordered lists — index-stable; slice into these in run scripts
# ---------------------------------------------------------------------------
ALL_CO_CUTS   = [None, 'gamma10', 'v150', 'gamma10_v150']
ALL_NDBD_CUTS = [None, 'u1775', 'u1740']
ALL_CO_SHIFTS = [0.0, 0.3, 0.5, 0.8, 1.0, 1.2, 1.5]
ALL_FIT_TYPES = ['toys', 'sens', 'bias']

# ---------------------------------------------------------------------------
# Default grid for batch job submission (gen_job_params.py etc.)
# ---------------------------------------------------------------------------
CO_CUTS   = [None, 'gamma10_v150', 'v150']
NDBD_CUTS = [None, 'u1740', 'u1775']
CO_SHIFTS = ALL_CO_SHIFTS
FIT_TYPES = ['toys']
