import numpy as np

import site
site.addsitedir('/global/homes/x/xcorat/Software//QPyM2/')

from qpym2.utils import debug

_FNQ_NAME = 'fit_function_Q'
_FNS_NAME = 'fit_function_sigma'
_COVQ_NAME = 'covariance_matrix_Q'
_COVS_NAME = 'covariance_matrix_sigma'

def sample_scaling(rf, fn_name, cov_name, skip_last_mean=False, update=False):
    fn = rf.Get(fn_name)
    cov = rf.Get(cov_name)
    means = np.array([ fn.GetParameter(i) for i in range(fn.GetNpar()) ])
    debug(f'read means ({fn_name}): {means}')
    means = means[:-1] if skip_last_mean else means
    cv = np.array([ [ cov(i, j) for j in range(cov.GetNrows()) ] for i in range(cov.GetNcols()) ])

    sampled = np.random.multivariate_normal(means, cv)

    if skip_last_mean:
        sampled = np.append(sampled, means[-1])
    if update:
        debug('Updating the file...')
        fn.SetParameters(*sampled)
        rf.Write()
    debug(f'sampled means ({fn_name}): {sampled}')
    return sampled

if __name__ == '__main__':
    """
    
    save scripts:
    1. copy the config data folder to multiple folders for sampling
        ```
        $ pwd
        /global/homes/x/xcorat/Software/cuore-nersc-modern_v4.0.0/cuoremc/ares/data
        for i in {1..10}; do cp -r Nature2021 -p m2/s${i}; done;
        ```
    2. run the sampling script (HINT: update)
        ```
        $ pwd
        /global/homes/x/xcorat/Software/QPyM2
        do python /global/homes/x/xcorat/Software/QPyM2/qpym2/run/sample_lsscaling_ares.py ../cuore-nersc-modern_v4.0.0/cuoremc/ares/data/m2/s${i}/residual_and_width_vs_energy_ds36*.root | tee >> ../cuore-nersc-modern_v4.0.0/cuoremc/ares/data/m2/s${i}/run_sample_lsscaling.log; done;
        ```  

    """
    from ROOT import TFile
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('file', help='path to the fit file/s.', nargs='+')
    parser.add_argument('--update', action='store_true', help='update the fit file')
    args = parser.parse_args()

    debug(args)

    for fname in args.file:
        debug('Opening file: ', fname)
        f = TFile(fname, 'UPDATE')
        q = sample_scaling(f, _FNQ_NAME, _COVQ_NAME, skip_last_mean=True, update=args.update)
        sig = sample_scaling(f, _FNS_NAME, _COVS_NAME, update=args.update)
        print(q, sig)
        f.Close()
