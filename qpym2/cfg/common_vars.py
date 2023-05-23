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

# datapath = '/global/project/projectdirs/cuore/scratch/xcorat/data'
datapath = '/global/cfs/projectdirs/cuore/scratch/xcorat/data'
ares_outpath = f'{datapath}/BM_Nature2021/out_ares'
stagedir = f'{datapath}/MC/staging_nat21_roi310_unblinded'
fitdir = '/global/homes/x/xcorat/Software//QPyM2/pdata/nat21_roi310_unblinded'

staging_config = {
    'outpath': stagedir,
    'outtree_name': 'uvtree',
    # 'fname_assembler': lambda fname, shift_name: f'{fname}_{shift_name}',
    # TODO: this is named 'out_cols' 
    'out_cols_mc' : ['u', 'v', 'Ch1', 'Ch2', 'Dataset', 'PCA', 'MultipletValidation'],
    'out_cols_data' : ['u', 'v', 'Ch1', 'Ch2', 'Dataset', 'PCA', 'PCA8', 'Validation', 'Validation8']
}

read_config = {
    'inpath': ares_outpath,
    'treename': 'outTree',
    'evar': 'Energy',
    'esum_var': 'TotalEnergy',
    'ch1': 'Channel',
    'ch2': 'Multiplet[1]'
}

read_config_0v = {
    'inpath':'',
    'inpath_full': f'{datapath}/MC/ares/nat21_highstats/0vbb/0vbb_1M_*.root',
    'treename': 'outTree',
    'evar': 'Energy',
    'esum_var': 'TotalEnergy',
    'ch1': 'Channel',
    'ch2': 'Multiplet[1]'
}

# --------------- Staging model ---------------
hcfg = {
    'esum_max': 2700,
    'esum_min': 2390,
    'threshold': 340
}
hist_cuts = f"{read_config['esum_var']}>{hcfg['esum_min']} && {read_config['esum_var']}<{hcfg['esum_max']}"
hist_cuts += f" && {read_config['evar']}>{read_config['esum_var']}/2 && {read_config['esum_var']}-{read_config['evar']} > {hcfg['threshold']}"

cuts = "AllFilters && Multiplicity == 2 && " + hist_cuts 
 # --------------- End staging model -----------


s2 = math.sqrt(2)
threshold = 350
esummin, esummax = 2400, 2570
vmin, vmax = 0, (esummax-2*threshold)/s2
smooths = [None, (1, 'k5b')]
smooth = smooths[1]

ROI_MIN, ROI_MAX = 2515, 2540
## IMPORTANT
skip_roi = None # (ROI_MIN, ROI_MAX)

# NOTE: is this e1>e2 implented elsewhere? mcuv creation?
ecuts = f"u>{esummin}/sqrt(2) && u<{esummax}/sqrt(2) && e1>e2"
cut_350 = 'e2 > 350'
cut_rmCo2_10 = f'abs(e1-{1334})>{10} && abs(e2-{1173})>{10}'
cut_0v1775 = f'u>1775'
cut_0v1740 = f'u>1740'
cut_v150 = f'v>150'


run_configs = \
{
    'noco2_no0v1775_v150_smooth': 
    {
        'cuts': {
            'run': f'{ecuts}&&{cut_350}&&{cut_rmCo2_10}&&{cut_v150}',
            '0v': f'{cut_0v1775} && PCA && MultipletValidation', 
            'mc': f'PCA && MultipletValidation',
            'data': f'PCA && Validation'
        },
        'title': f'All floors, -Co peaks (10keV), -0v < 2509, #DeltaE<150*#sqrt 2',
        'smooth': smooth,
    },
    'noco2_no0v1740_v150_smooth': 
    {
        'cuts': {
            'run': f'{ecuts}&&{cut_350}&&{cut_rmCo2_10}&&{cut_v150}',
            '0v': f'{cut_0v1740} && PCA && MultipletValidation', 
            'mc': f'PCA && MultipletValidation',
            'data': f'PCA && Validation'
        },
        'title': f'All floors, -Co peaks (10keV), -0v < 2440, #DeltaE<150*#sqrt 2',
        'smooth': smooth,
    },
    'noco2_full0v_v150_smooth': 
    {
        'cuts': {
            'run': f'{ecuts}&&{cut_350}&&{cut_rmCo2_10}&&{cut_v150}',
            '0v': f'PCA && MultipletValidation', 
            'mc': f'PCA && MultipletValidation',
            'data': f'PCA && Validation'
        },
        'title': f'All floors, -Co peaks (10keV), #DeltaE<150*#sqrt 2',
        'smooth': smooth,
    },
    'no0v1740_v150_smooth': 
    {
        'cuts': {
            'run': f'{ecuts}&&{cut_350}&&{cut_v150}',
            '0v': f'{cut_0v1740} && PCA && MultipletValidation', 
            'mc': f'PCA && MultipletValidation',
            'data': f'PCA && Validation'
        },
        'title': f'All floors, -0v < 2440, #DeltaE<150*#sqrt 2',
        'smooth': smooth,
    },
    'no0v1740_smooth': 
    {
        'cuts': {
            'run': f'{ecuts}&&{cut_350}',
            '0v': f'{cut_0v1740} && PCA && MultipletValidation', 
            'mc': f'PCA && MultipletValidation',
            'data': f'PCA && Validation'
        },
        'title': f'All floors, -0v < 2440',
        'smooth': smooth,
    },
}

cut_roi = f'!(u>{ROI_MIN}/sqrt(2) && u<{ROI_MAX}/sqrt(2))'
run_configs_noroi = \
{
    'noco2_no0v1775_v150_smooth_noroi': 
    {
        'cuts': {
            'run': f'{ecuts}&&{cut_350}&&{cut_rmCo2_10}&&{cut_v150}&&{cut_roi}',
            '0v': f'{cut_0v1775} && PCA && MultipletValidation', 
            'mc': f'PCA && MultipletValidation',
            'data': f'PCA && Validation'
        },
        'title': f'All floors, -Co peaks (10keV), -0v < 2509, #DeltaE<150*#sqrt 2',
        'smooth': smooth,
    },
    'noco2_no0v1740_v150_smooth_noroi': 
    {
        'cuts': {
            'run': f'{ecuts}&&{cut_350}&&{cut_rmCo2_10}&&{cut_v150}&&{cut_roi}',
            '0v': f'{cut_0v1740} && PCA && MultipletValidation', 
            'mc': f'PCA && MultipletValidation',
            'data': f'PCA && Validation'
        },
        'title': f'All floors, -Co peaks (10keV), -0v < 2440, #DeltaE<150*#sqrt 2',
        'smooth': smooth,
    },
    'noco2_full0v_v150_smooth_noroi': 
    {
        'cuts': {
            'run': f'{ecuts}&&{cut_350}&&{cut_rmCo2_10}&&{cut_v150}&&{cut_roi}',
            '0v': f'PCA && MultipletValidation', 
            'mc': f'PCA && MultipletValidation',
            'data': f'PCA && Validation'
        },
        'title': f'All floors, -Co peaks (10keV), #DeltaE<150*#sqrt 2',
        'smooth': smooth,
    },
    'no0v1740_v150_smooth_noroi': 
    {
        'cuts': {
            'run': f'{ecuts}&&{cut_350}&&{cut_v150}&&{cut_roi}',
            '0v': f'{cut_0v1740} && PCA && MultipletValidation', 
            'mc': f'PCA && MultipletValidation',
            'data': f'PCA && Validation'
        },
        'title': f'All floors, -0v < 2440, #DeltaE<150*#sqrt 2',
        'smooth': smooth,
    },
    'no0v1740_smooth_noroi': 
    {
        'cuts': {
            'run': f'{ecuts}&&{cut_350}&&{cut_roi}',
            '0v': f'{cut_0v1740} && PCA && MultipletValidation', 
            'mc': f'PCA && MultipletValidation',
            'data': f'PCA && Validation'
        },
        'title': f'-0v < 2440',
        'smooth': smooth,
    }
}

hm = TempConfig(name='hmu', threshold=threshold, 
                umin=esummin/s2, umax=esummax/s2,
                vmin=vmin, vmax=vmax,
                binsize=1/s2, skip_roi=None, hist_type='h2uv')
hm_noroi = TempConfig(name='hmu', threshold=threshold, 
                umin=esummin/s2, umax=esummax/s2,
                vmin=vmin, vmax=vmax,
                binsize=1/s2, skip_roi=(ROI_MIN, ROI_MAX), hist_type='h2uv')

shifts = {
    'noshift': { 'fname_postfix' : 'noshift', 'shift': 0 },
    'ushiftp03': { 'fname_postfix' : 'ushiftp03', 'shift': 0.3, 'transform': ( '(esum + 0.0003*e1)/sqrt(2)', '(2*e1 + 0.0003*e1 - esum)/sqrt(2)') },
    'ushiftp05': { 'fname_postfix' : 'ushiftp05', 'shift': 0.5, 'transform': ( '(esum + 0.0005*e1)/sqrt(2)', '(2*e1 + 0.0005*e1 - esum)/sqrt(2)') },
    'ushiftp1': { 'fname_postfix' : 'ushiftp1',  'shift': 1,  'transform': ( '(esum + 0.001 *e1)/sqrt(2)', '(2*e1 + 0.001*e1 - esum)/sqrt(2)') },
    'ushiftp15': { 'fname_postfix' : 'ushiftp15', 'shift': 1.5, 'transform': ( '(esum + 0.0015*e1)/sqrt(2)', '(2*e1 + 0.0015*e1 - esum)/sqrt(2)') },
}
cfg_m2nat21_may23 = TempConfig(
    outdir =  fitdir,
    datapath = datapath,
    staging_config = staging_config,
    jags_fpath = f'{datapath}/BM_Nature2021/out_jags/Out.h5',
    run_configs = run_configs,
    shifts = shifts,
    hm = hm,
    smooth = smooth
)

global_vars_noroi = TempConfig(
    outdir =  fitdir,
    datapath = datapath,
    staging_config = staging_config,
    jags_fpath = '/global/cfs/projectdirs/cuore/syncData/CUORE/simulation/Nature2021/out_jags/Out.h5',
    run_configs = run_configs_noroi,
    shifts = shifts,
    hm = hm_noroi,
    smooth = smooth
)