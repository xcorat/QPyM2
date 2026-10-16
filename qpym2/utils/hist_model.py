
# def hist_model(esum_min, esum_max, threshold):
#     """ Return a hist model dictionary with the given parameters.
#     """
#     s2 = math.sqrt(2)
#     vmin, vmax = 0, (esum_max-2*threshold)/s2
#     return {
#         'esum_max': esum_max,
#         'esum_min': esum_min,
#         'threshold': threshold,
#         'vmin': vmin,
#         'vmax': vmax,
#     }

# def get_h2uv_cuts(esum_min=2400, esum_max=2700, threshold=350, inpars=('u', 'v')):
#     """ Return the hist cuts string for the given parameters.
#     """
#     _u, _v = inpars
#     ecuts = f"{_u}>{esum_min}/sqrt(2) && {_u}<{esum_max}/sqrt(2) && {_v}>0" # TODO: check e1>e2, 
#     ecuts += f'&& ({_u}-{_v})/sqrt(2) > {threshold}' # TODO: check e2>350
    # return ecuts

import math
from qpym2.utils.hists import HistType

def get_hm_cuts(hm, inpars=('u', 'v')):
    """ Return the hist cuts string for the given parameters.
    """
    _u, _v = inpars
    s2 = math.sqrt(2)

    if hm.hist_type == HistType.h2ee:
        esum_min = hm.e1min*2
        esum_max = hm.e2max*2
    elif hm.hist_type == HistType.h2uv:
        esum_min = hm.umin*s2
        esum_max = hm.umax*s2
        
    else:
        raise ValueError(f'Invalid hist type: {hm.hist_type}. Expected h2uv only.')

    ecuts = f"{_u}>{esum_min}/sqrt(2) && {_u}<{esum_max}/sqrt(2) && {_v}>0" # TODO: check e1>e2, 
    ecuts += f'&& ({_u}-{_v})/sqrt(2) > {hm.threshold}' # TODO: check e2>350

    if hm.skip_roi:
        roimin, roimax = hm.skip_roi
        roicuts =  f'!({_u}>{roimin}/sqrt(2) && {_u}<{roimax}/sqrt(2))'
        ecuts += f'&& {roicuts}'

    print(ecuts)
    return ecuts

def get_qroot_histcuts(esum_var='TotalEnergy', evar='Energy',
                       esum_min=2390, esum_max=2710, threshold=340):
    """ Return the hist cuts string for the given parameters.
    """
    hist_cuts = f"{esum_var}>{esum_min} && {esum_var}<{esum_max}"
    hist_cuts += f" && {evar}>{esum_var}/2 && {esum_var}-{evar} > {threshold}"
    return hist_cuts