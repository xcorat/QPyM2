import site
site.addsitedir('/global/homes/x/xcorat/Software//QPyM2/')
site.addsitedir('/global/homes/x/xcorat/Software//QPyMC/')

from qpym2.io import staging

datapath = '/global/cfs/projectdirs/cuore/scratch/xcorat/data'
# datapath = '/global/project/projectdirs/cuore/scratch/xcorat/data'
ares_outpath = f'{datapath}/BM_Nature2021/out_ares'

 
out_cols_mc = ['u', 'v', 'Ch1', 'Ch2', 'Dataset', 'PCA', 'MultipletValidation']
out_cols_data = ['u', 'v', 'Ch1', 'Ch2', 'Dataset', 'PCA', 'PCA8', 'Validation', 'Validation8']

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

write_config = {
    'outpath': f'{datapath}/MC/mcuvs_roi310',
    'fname_assembler': lambda fname, shift_name: f'{fname}_{shift_name}',
    'outtree_name': 'uvtree',
    'out_cols': out_cols_mc,
}

jags_fpath = f'{datapath}/BM_Nature2021/out_jags/Out.h5'


shifts = {
    'noshift': { 'fname_postfix' : 'noshift', 'shift': 0 },
    # 'ushiftp03': { 'fname_postfix' : 'ushiftp03', 'shift': 0.3, 'transform': ( '(esum + 0.0003*e1)/sqrt(2)', '(2*e1 + 0.0003*e1 - esum)/sqrt(2)') },
    # 'ushiftp05': { 'fname_postfix' : 'ushiftp05', 'shift': 0.5, 'transform': ( '(esum + 0.0005*e1)/sqrt(2)', '(2*e1 + 0.0005*e1 - esum)/sqrt(2)') },
    # 'ushiftp1': { 'fname_postfix' : 'ushiftp1',  'shift': 1,  'transform': ( '(esum + 0.001 *e1)/sqrt(2)', '(2*e1 + 0.001*e1 - esum)/sqrt(2)') },
    # 'ushiftp15': { 'fname_postfix' : 'ushiftp15', 'shift': 1.5, 'transform': ( '(esum + 0.0015*e1)/sqrt(2)', '(2*e1 + 0.0015*e1 - esum)/sqrt(2)') },
}

#  ---------------


hcfg = {
    'esum_max': 2700,
    'esum_min': 2390,
    'threshold': 340
}
hist_cuts = f"{read_config['esum_var']}>{hcfg['esum_min']} && {read_config['esum_var']}<{hcfg['esum_max']}"
hist_cuts += f" && {read_config['evar']}>{read_config['esum_var']}/2 && {read_config['esum_var']}-{read_config['evar']} > {hcfg['threshold']}"

cuts = "AllFilters && Multiplicity == 2 && " + hist_cuts 

from qpym2.io import bkg_model

def create_staging(shifts):
    pass
    table, _ = bkg_model.read_jagsh5(jags_fpath)
    default_truv = ('esum/sqrt(2)', '(2*e1 - esum)/sqrt(2)')

    table_co60 = table.filter(like='60Co', axis=0)


    table_others = table.drop(table_co60.index, axis=0)
    staging_other_transforms = {'noshift':  [
        ('filter', cuts), ('define', 'u', default_truv[0]), ('define', 'v', default_truv[1])]}

    print("creating staging files: ")
    ret = staging.create_staging_multi(table_others, read_config=read_config, write_config=write_config,
                            staging_multi_transforms=staging_other_transforms)
    
    staging_co60_transforms = {}
    for k, v in shifts.items():
        shift_tr = v.get('transform', default_truv)
        staging_co60_transforms[k] =  [('filter', cuts),
                                        ('define', 'u', shift_tr[0]),
                                        ('define', 'v', shift_tr[1])]

    staging.create_staging_multi(table_co60, read_config=read_config, write_config=write_config,
                            staging_multi_transforms=staging_co60_transforms)

    return ret

def create_staging_data(shifts):
    pass
    from qpym2.io.staging import write_blinded_datalist

    data_fname = '/global/project/projectdirs/cuore/scratch/xcorat/data//BM_Nature2021/Reduced_Nature21.root'
    datacuts = 'Multiplicity==2 && cutsTree.Included && cutsTree.FullMultiplet'
    nbb_fname = f'{datapath}/MC/ares/nat21_highstats/0vbb/0vbb_1M_*.root'

    cut_0v1740 = f'TotalEnergy>1740*sqrt(2)'
    cut_0v1775 = f'TotalEnergy>1775*sqrt(2)'

    print(shift_cfg)

    mccuts = f'AllFilters == 1'

    for k, v in shifts.items():
        shift_cfg = v
        shift_u = shift_cfg['shift']/2500
        print(shift_u)
        outpath_0v1740 = f'{datapath}/MC/mcuvs_roi300/data_0v1740_{shift_cfg["fname_postfix"]}.root'
        print('out :', outpath_0v1740)
        df = write_blinded_datalist(data_fname, nbb_fname, outpath_0v1740, hist_cuts, 150, 200, mccuts=f'{mccuts}&&{cut_0v1740}', datacuts=datacuts, out_cols=out_cols_data, data_treename='outTree', dscale=shift_u)
        outpath_0v1775 = f'{datapath}/MC/mcuvs_roi300/data_0v1775_{shift_cfg["fname_postfix"]}.root'
        print('out :', outpath_0v1775)
        df = write_blinded_datalist(data_fname, nbb_fname, outpath_0v1775, hist_cuts, 150, 200, mccuts=f'{mccuts}&&{cut_0v1775}', datacuts=datacuts, out_cols=out_cols_data, data_treename='outTree', dscale=shift_u)
        outpath_0vfull = f'{datapath}/MC/mcuvs_roi300/data_{shift_cfg["fname_postfix"]}.root'
        print('out :', outpath_0vfull)
        df = write_blinded_datalist(data_fname, nbb_fname, outpath_0vfull, hist_cuts, 150, 200, mccuts=f'{mccuts}', datacuts=datacuts, out_cols=out_cols_data, data_treename='outTree', dscale=shift_u)

    outpath_noroi = f'{datapath}/MC/mcuvs_full/data_noroi_noshift.root'
    cut_roi = f'!(TotalEnergy>2515 && TotalEnergy<2540)'
    datacuts_noroi = f'{datacuts} && {cut_roi}'
    df = write_blinded_datalist(data_fname, nbb_fname, outpath_noroi, hist_cuts, 0, 0, mccuts=f'{mccuts}&&{cut_0v1775}', datacuts=datacuts_noroi, out_cols=out_cols_data, data_treename='outTree', dscale=0)

def create_staging_0v(shifts):
    import pandas as pd
    ndbd_info = pd.Series({'name': 'ndbd', 'fname': '0vbb.root', 'mean': 0, 'std': 0, 'mode': 0})
    ndbd_table = pd.DataFrame([ndbd_info]).set_index('name')

    default_truv = ('esum/sqrt(2)', '(2*e1 - esum)/sqrt(2)')
    staging_multi_transforms = {}
    for k, v in shifts.items():
        shift_tr = v.get('transform', default_truv)
        staging_multi_transforms[k] =  [('filter', cuts),
                                        ('define', 'u', shift_tr[0]),
                                        ('define', 'v', shift_tr[1])]
        
    staging.create_staging_multi(ndbd_table, read_config=read_config_0v, write_config=write_config,
                             staging_multi_transforms=staging_multi_transforms)


    # from qpym2.io import staging#, write_UV_evlist_fakePCA
    # #cfg.rootpath + source.relpath + source.fname
    # fpath = read_config_0v['inpath_full']
        
    # outpath_0v = f'{datapath}/MC/mcuvs_roi310/0vbb_{shift_name}_test.root'
    # """ Selected: 2574857 """
    # staging.write_UV_evlist(outpath_0v, fpath, datacuts=f'AllFilters == 1&&{cuts}', out_cols=out_cols_mc, shift_u=0)

def create_staging_data():
    from qpym2.io.staging import write_blinded_datalist

    # data_fname = '/global/cfs/projectdirs/cuore/scratch/xcorat/data//BM_Nature2021/Reduced_Nature21.root'
    data_fname = '/global/cfs/projectdirs/cuore/syncData/CUORE/simulation/Nature2021/Reduced_Nature21.root'
    datacuts = 'Multiplicity==2 && cutsTree.Included && cutsTree.FullMultiplet'
    nbb_fname = f'{datapath}/MC/ares/nat21_highstats/0vbb/0vbb_1M_*.root'


    mccuts = f'AllFilters == 1'
    outpath_unb = f'{datapath}/MC/mcuvs_roi310/data_unblinded.root'
    print(datapath)
    # return
    df = write_blinded_datalist(data_fname, nbb_fname, outpath_unb, hist_cuts, 0, 0, mccuts=f'{mccuts}', datacuts=datacuts, out_cols=out_cols_data, data_treename='outTree', dscale=0)


if __name__ == '__main__':
    # create_staging()
    create_staging_data()
    # create_staging_0v(shifts=shifts)