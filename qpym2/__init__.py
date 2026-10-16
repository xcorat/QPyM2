"""
QPyM2: CUORE Multiplicity 2 Analysis Framework

This package provides tools and utilities for analyzing CUORE Multiplicity 2 data.

Main modules:
    - qpym2.cfg: Configuration management
    - qpym2.fit: Fitting and statistical analysis
    - qpym2.io: Data input/output utilities  
    - qpym2.plotting: Visualization tools
    - qpym2.utils: General utilities
    - qpym2.guide: Interactive guide utilities for notebooks

For interactive analysis, see the notebooks under notebooks/ (start with
notebooks/tech_note_setup/). Plotting is library-only in qpym2.plotting and is
driven from those notebooks.
"""

__version__ = "0.5.0"
__author__ = "QPyM2 Team"

# Note: We don't import guide utilities at the top level to avoid
# heavy dependencies for users who don't need the interactive features.
# Use: from qpym2.guide import ... for guide utilities.