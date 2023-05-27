# CUORE Python M2 Analysis

This repo is arranged as follows,

data : input and output for fitting,
  - testname1
  - groupname2
      - input.h5 : staged data and MC -  stored as dataframes
      - fit.h5 : fit output and component histograms. (Everything needed to replicate a fit?)
      - fit.nc : arviz inference data

qpym2 : python package that contains most of the code and sub packages.
  - io : directory for all the low level IO code implemented as functions.
  - run : directory for all the high level run scripts.
  - cfg : configurable settings for analysis.

doc : Documentation for the M2 analysis software and references. Built using Sphinx.

tests : tests and other related code, playground

### Requirements 

- root, numpy, pandas, pymc, h5py, pytables, sphinx

### Setting up documentation

Requirements:
  * sphinx, autodoc2, myst-parser, pip-tools

The documentation and reference manual is built using [Sphinx] and [autodoc2]. 
The documentation is build automatically at [readthedocs.io](https://qpym2.readthedocs.io) when the `ref_docs` branch is updated. To build and setup (for a new project only,)

1. Install sphinx and autodoc2 if not already installed.
2. To generate the initial files, run `sphinx-quickstart` in the `doc` directory.
3. Edit the `doc/conf.py` and add:
    + extensions `myst-parser` and `autodoc2`
    + autodoc and sphinx parameters (including the path to the python package)
4. Create the `docs/requirements.in` file and add the required packages, and
    run `pip-compile` to generate the `docs/requirements.txt` file.



 
