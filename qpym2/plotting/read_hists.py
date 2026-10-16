from qpym2.fit.create_m2uv import get_mctable, create_evlist
from qpym2.models.m2 import make_bkgcomps_table
from qpym2.utils.hists import get_sum_hist, get_empty_hist

from qpym2.utils import hists

def mch2sum(staging_cfg, fit_cfg, jags_fpath, comp_cfg, signal_nevs=False):
    """ Read the staging files and create the sum histogram accordingly 
    
        TODO: this need to be refactored along with the `get_mctable`
    """

    staging_comps = get_mctable(fit_cfg, staging_cfg, jags_fpath)
    data = create_evlist(fit_cfg, staging_cfg)

    smooth = fit_cfg['h2uv_model'].smooth
    fit_comps = make_bkgcomps_table(comp_cfg, staging_comps[:-1], smooth)
    h2uv_model = fit_cfg['h2uv_model']

    if signal_nevs:
        signal_df = staging_comps.loc['ndbd'].copy()
        sig_hist = signal_df['mchist']

        if smooth:
            sig_hist = hists.smooth_nph2(sig_hist, smooth=smooth)

        signal_df['pdf_hist'] = sig_hist/signal_df['integral'] 
        if type(signal_nevs) == int:
            signal_df.integral = signal_nevs
        else:
            signal_df.integral = staging_cfg['ndbd_injected']

        fit_comps.loc['ndbd'] = signal_df
    
    else:
        fit_comps = fit_comps[fit_comps.index != 'ndbd']

    # from qpym2.utils.hists import get_sum_hist, get_empty_hist
    h2_fit = get_sum_hist(fit_comps, fit_comps.integral)
    # _, xedges, yedges = get_empty_hist(hm_noroi)
    th2d_fit = get_empty_hist(h2uv_model, return_numpy=False)

    for i in range(h2_fit.shape[0]):
        for j in range(h2_fit.shape[1]):
            th2d_fit.SetBinContent(i+1, j+1, h2_fit[i, j])

    th2d_data = th2d_fit.Clone('th2d_fit')
    th2d_data.Reset('RSEM')
    for ev in data:
        th2d_data.Fill(ev[0], ev[1])

    return th2d_fit, th2d_data