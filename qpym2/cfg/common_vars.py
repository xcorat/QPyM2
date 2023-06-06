""" All the configuration and settings for the analysis are
    defined here. This absolutely needs to be refactored, 

    + Most of the variables need to be moved to a config file (not python)
    + Different settings-sets per configuration can be moved to a single
        file that the main script import according to the analysis,

    etc.. etc...

    This is a mess, but it works for now. [<- says copilot]
"""
import math

from qpym2.utils import TempConfig
from qpym2.hists import HistType

# datapath = '/global/project/projectdirs/cuore/scratch/xcorat/data'
datapath = '/global/cfs/projectdirs/cuore/scratch/xcorat/data'
ares_outpath = f'{datapath}/BM_Nature2021/out_ares'
stagedir = f'{datapath}/MC/staging_nat21_roi310_unblinded'
fitdir = '/global/homes/x/xcorat/Software//QPyM2/data/nat21_roi310_unblinded'

# ---------------- Cuts ---------------
cuts = "PCA && Multiplicity == 2"
datacuts = 'Validation && cutsTree.Included && cutsTree.FullMultiplet'
mccuts = 'AllFilters && MultipletValidation'
staging_config = {
    'outpath': stagedir,
    'treename': 'uvtree',
    # 'fname_assembler': lambda fname, shift_name: f'{fname}_{shift_name}',
    # TODO: this is named 'out_cols' 
    'out_cols' : ['u0', 'v0', 'ch1', 'ch2' ], #, 'Dataset', 'PCA', 'MultipletValidation'],
    # 'out_cols_data' : ['u0', 'v0', 'ch1', 'ch2' ], #, 'Dataset', 'PCA', 'PCA8', 'Validation', 'Validation8']
}

read_config = {
    'inpath': ares_outpath,
    'treename': 'outTree',
    'evar': 'Energy',
    'esum_var': 'TotalEnergy',
    'ch1': 'Channel',
    'ch2': 'Multiplet[1]',
    'extra_cuts': mccuts,
}

read_config_0v = {
    # 'inpath_full': f'{datapath}/MC/ares/nat21_highstats/0vbb/0vbb_1M_*.root',
    'inpath':f'{datapath}/MC/ares/nat21_highstats/0vbb',
    'fname': '0vbb_1M_*.root',
    'treename': 'outTree',
    'evar': 'Energy',
    'esum_var': 'TotalEnergy',
    'ch1': 'Channel',
    'ch2': 'Multiplet[1]',
    'extra_cuts': mccuts,
}

# ---------------- End staging data -----------
read_config_data = {
    'inpath': '/global/cfs/projectdirs/cuore/scratch/xcorat/data//BM_Nature2021',
    'fname': 'Reduced_Nature21.root',
    'treename': 'outTree',
    'evar': 'Energy',
    'esum_var': 'TotalEnergy',
    'ch1': 'Channel',
    'ch2': 'Multiplet[1]',
    'out_fpath': f'{stagedir}/data_unblinded.root',
    'extra_cuts': datacuts,
}


# --------------- Staging model ---------------
hcfg = {
    'esum_max': 2700,
    'esum_min': 2390,
    'threshold': 340
}
# convert to make_hist_cut(read_cfg, hcfg)
hist_cuts = f"{read_config['esum_var']}>{hcfg['esum_min']} && {read_config['esum_var']}<{hcfg['esum_max']}"
hist_cuts += f" && {read_config['evar']}>{read_config['esum_var']}/2 && {read_config['esum_var']}-{read_config['evar']} > {hcfg['threshold']}"

cuts = cuts + " && " + hist_cuts 
# --------------- End staging model -----------

# --------------- Fit model ---------------
esum_min, esum_max = 2400, 2570
s2 = math.sqrt(2)
threshold = 350
vmin, vmax = 0, (esum_max-2*threshold)/s2
hist_model_fit = {
    'esum_max': esum_max,
    'esum_min': esum_min,
    'threshold': 350
}

smooths = [None, (1, 'k5b')]
binsizes = [1]
smooth = smooths[1]
binsize = binsizes[0]

ROI_MIN, ROI_MAX = 2515, 2540
## IMPORTANT
# skip_roi = None # (ROI_MIN, ROI_MAX)
# TODO: take the hist_type from hist enum
hm = TempConfig(name='hmu', threshold=threshold, 
                umin=esum_min/s2, umax=esum_max/s2,
                vmin=vmin, vmax=vmax,
                binsize=binsize/s2, skip_roi=None, 
                smooth=smooth, hist_type=HistType.h2uv)
hm_noroi = TempConfig(name='hmu', threshold=threshold, 
                umin=esum_min/s2, umax=esum_max/s2,
                vmin=vmin, vmax=vmax,
                binsize=binsize/s2, skip_roi=(ROI_MIN, ROI_MAX),
                smooth=smooth, hist_type=HistType.h2uv)

# branches u, v, e1, e2, esum, are assumed to be present.
# NOTE: is this e1>e2 implented elsewhere? mcuv creation?
ecuts = f"u>{esum_min}/sqrt(2) && u<{esum_max}/sqrt(2) && e1>e2"
ecuts += '&& e2 > 350'

from enum import Enum

fit_cuts_co = {
    '_v150': f'v>150',
    '_cosum_10': f'abs(e1-{1334})>{10} && abs(e2-{1173})>{10}',
}

fit_cuts_ndbd = {
    '_0v1775': f'u>1775',
    '_0v1740': f'u>1740'
}

shift_types = [ 'ushiftp' ]
shifts = [0, 0.3, 0.5, 1, 1.5]

fit_cuts_common = f'{ecuts}&&{fit_cuts_co["_v150"]}&&{fit_cuts_co["_cosum_10"]}'


fit_cuts = TempConfig(
    common = fit_cuts_common,
    mc_only = None,
    ndbd_only = fit_cuts_ndbd['_0v1740'],
    data_only = None,
)

# fit_input_cfg = TempConfig(
#     fit_dir = fitdir,

fit_pars = TempConfig(
    name = 'noco2_no0v1740_v150_smooth',
    fit_dir = fitdir,
    cuts = fit_cuts,
    shift = TempConfig(type=shift_types[0], val=shifts[2], filter_str='co60')
)
    
# cut_rmCo2_10 = f'abs(e1-{1334})>{10} && abs(e2-{1173})>{10}'
# cut_0v1775 = f'u>1775'
# cut_0v1740 = f'u>1740'
# cut_v150 = f'v>150'

# TODO: 6/1/23 - we don't need the PCA and MultipletValidation cuts on staging data.
# run_configs = \
# {
#     'noco2_no0v1775_v150_smooth': 
#     {
#         'cuts': {
#             'run': f'{ecuts}&&{cut_350}&&{cut_rmCo2_10}&&{cut_v150}',
#             '0v': f'{cut_0v1775} && PCA && MultipletValidation', 
#             'mc': f'PCA && MultipletValidation',
#             'data': f'PCA && Validation'
#         },
#         'title': f'All floors, -Co peaks (10keV), -0v < 2509, #DeltaE<150*#sqrt 2',
#         'smooth': smooth,
#     },
#     'noco2_no0v1740_v150_smooth': 
#     {
#         'cuts': {
#             'run': f'{ecuts}&&{cut_350}&&{cut_rmCo2_10}&&{cut_v150}',
#             '0v': f'{cut_0v1740} && PCA && MultipletValidation', 
#             'mc': f'PCA && MultipletValidation',
#             'data': f'PCA && Validation'
#         },
#         'title': f'All floors, -Co peaks (10keV), -0v < 2440, #DeltaE<150*#sqrt 2',
#         'smooth': smooth,
#     },
#     'noco2_full0v_v150_smooth': 
#     {
#         'cuts': {
#             'run': f'{ecuts}&&{cut_350}&&{cut_rmCo2_10}&&{cut_v150}',
#             '0v': f'PCA && MultipletValidation', 
#             'mc': f'PCA && MultipletValidation',
#             'data': f'PCA && Validation'
#         },
#         'title': f'All floors, -Co peaks (10keV), #DeltaE<150*#sqrt 2',
#         'smooth': smooth,
#     },
#     'no0v1740_v150_smooth': 
#     {
#         'cuts': {
#             'run': f'{ecuts}&&{cut_350}&&{cut_v150}',
#             '0v': f'{cut_0v1740} && PCA && MultipletValidation', 
#             'mc': f'PCA && MultipletValidation',
#             'data': f'PCA && Validation'
#         },
#         'title': f'All floors, -0v < 2440, #DeltaE<150*#sqrt 2',
#         'smooth': smooth,
#     },
#     'no0v1740_smooth': 
#     {
#         'cuts': {
#             'run': f'{ecuts}&&{cut_350}',
#             '0v': f'{cut_0v1740} && PCA && MultipletValidation', 
#             'mc': f'PCA && MultipletValidation',
#             'data': f'PCA && Validation'
#         },
#         'title': f'All floors, -0v < 2440',
#         'smooth': smooth,
#     },
# }

# cut_roi = f'!(u>{ROI_MIN}/sqrt(2) && u<{ROI_MAX}/sqrt(2))'
# run_configs_noroi = \
# {
#     'noco2_no0v1775_v150_smooth_noroi': 
#     {
#         'cuts': {
#             'run': f'{ecuts}&&{cut_350}&&{cut_rmCo2_10}&&{cut_v150}&&{cut_roi}',
#             '0v': f'{cut_0v1775} && PCA && MultipletValidation', 
#             'mc': f'PCA && MultipletValidation',
#             'data': f'PCA && Validation'
#         },
#         'title': f'All floors, -Co peaks (10keV), -0v < 2509, #DeltaE<150*#sqrt 2',
#         'smooth': smooth,
#     },
#     'noco2_no0v1740_v150_smooth_noroi': 
#     {
#         'cuts': {
#             'run': f'{ecuts}&&{cut_350}&&{cut_rmCo2_10}&&{cut_v150}&&{cut_roi}',
#             '0v': f'{cut_0v1740} && PCA && MultipletValidation', 
#             'mc': f'PCA && MultipletValidation',
#             'data': f'PCA && Validation'
#         },
#         'title': f'All floors, -Co peaks (10keV), -0v < 2440, #DeltaE<150*#sqrt 2',
#         'smooth': smooth,
#     },
#     'noco2_full0v_v150_smooth_noroi': 
#     {
#         'cuts': {
#             'run': f'{ecuts}&&{cut_350}&&{cut_rmCo2_10}&&{cut_v150}&&{cut_roi}',
#             '0v': f'PCA && MultipletValidation', 
#             'mc': f'PCA && MultipletValidation',
#             'data': f'PCA && Validation'
#         },
#         'title': f'All floors, -Co peaks (10keV), #DeltaE<150*#sqrt 2',
#         'smooth': smooth,
#     },
#     'no0v1740_v150_smooth_noroi': 
#     {
#         'cuts': {
#             'run': f'{ecuts}&&{cut_350}&&{cut_v150}&&{cut_roi}',
#             '0v': f'{cut_0v1740} && PCA && MultipletValidation', 
#             'mc': f'PCA && MultipletValidation',
#             'data': f'PCA && Validation'
#         },
#         'title': f'All floors, -0v < 2440, #DeltaE<150*#sqrt 2',
#         'smooth': smooth,
#     },
#     'no0v1740_smooth_noroi': 
#     {
#         'cuts': {
#             'run': f'{ecuts}&&{cut_350}&&{cut_roi}',
#             '0v': f'{cut_0v1740} && PCA && MultipletValidation', 
#             'mc': f'PCA && MultipletValidation',
#             'data': f'PCA && Validation'
#         },
#         'title': f'-0v < 2440',
#         'smooth': smooth,
#     }
# }


# shifts = {
#     'noshift': { 'fname_postfix' : 'noshift', 'shift': 0 },
#     'ushiftp03': { 'fname_postfix' : 'ushiftp03', 'shift': 0.3, 'transform': ( '(esum + 0.0003*e1)/sqrt(2)', '(2*e1 + 0.0003*e1 - esum)/sqrt(2)') },
#     'ushiftp05': { 'fname_postfix' : 'ushiftp05', 'shift': 0.5, 'transform': ( '(esum + 0.0005*e1)/sqrt(2)', '(2*e1 + 0.0005*e1 - esum)/sqrt(2)') },
#     'ushiftp1': { 'fname_postfix' : 'ushiftp1',  'shift': 1,  'transform': ( '(esum + 0.001 *e1)/sqrt(2)', '(2*e1 + 0.001*e1 - esum)/sqrt(2)') },
#     'ushiftp15': { 'fname_postfix' : 'ushiftp15', 'shift': 1.5, 'transform': ( '(esum + 0.0015*e1)/sqrt(2)', '(2*e1 + 0.0015*e1 - esum)/sqrt(2)') },
# }
cfg_m2nat21_may23 = TempConfig(
    ares_config = read_config,
    ares_config_0v = read_config_0v,
    data_config = read_config_data,
    staging_cuts = cuts,
    outdir =  fitdir,
    datapath = datapath,
    staging_config = staging_config,
    jags_fpath = f'{datapath}/BM_Nature2021/out_jags/Out.h5',
    hm = hm,
    fit_pars = fit_pars,
    # --
    # run_configs = run_configs,
    # shifts = shifts,
    # smooth = smooth
)


# global_vars_noroi = TempConfig(
#     outdir =  fitdir,
#     datapath = datapath,
#     staging_config = staging_config,
#     jags_fpath = '/global/cfs/projectdirs/cuore/syncData/CUORE/simulation/Nature2021/out_jags/Out.h5',
#     run_configs = run_configs_noroi,
#     shifts = shifts,
#     hm = hm_noroi,
#     smooth = smooth
# )