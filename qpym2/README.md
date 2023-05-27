qpym2: CUORE Python M2 analysis python package
================================================


### Requirements 

- root, numpy, pandas, pymc, h5py, pytables, sphinx

io : directory for all the low level IO code implemented as functions.
run : directory for all the high level run scripts.
cfg : configurable settings for analysis.

### Setting up documentation

Requirements:
    + sphinx, autodoc2, myst-parser, pip-tools

The documentation and reference manual is built using [Sphinx] and [autodoc2]. 
The documentation is build automatically at [readthedocs.io](qpym2.readthedocs.io) when the `ref_docs` branch is updated. To build and setup (for a new project only,)

1. Install sphinx and autodoc2 if not already installed.
2. To generate the initial files, run `sphinx-quickstart` in the `docs` directory.
3. Edit the `docs/conf.py` and add:
    + extensions `myst-parser` and `autodoc2`
    + autodoc and sphinx parameters (including the path to the python package)
4. Create the `docs/requirements.in` file and add the required packages, and
    run `pip-compile` to generate the `docs/requirements.txt` file.

