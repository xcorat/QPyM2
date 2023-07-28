from typing import Tuple
def residual_sampling_uv(rspars_mean, rspars_sample, invars=('u0', 'v0')) -> Tuple[str, str]:
    """ Get the uv transform for the shift type and value.

    Args:
        residual (qpym2.cfg.fit_pars.Residual): residual
        invars (tuple): input (u, v) to transform

    Returns:
        tuple: (u, v) transform
    """
    u, v = invars
    e1,e2 = f'({u}+{v})/sqrt(2)', f'({u}-{v})/sqrt(2)'  
    a, b, c = rspars_sample - rspars_mean

    tru = f'{u} + {a}*sqrt(2) + {b}*({u}) + {c}/sqrt(2)*(({u})*({u})+({v})*({v}))'
    trv = f'{v}               + {b}*({v}) + {c}/sqrt(2)*({u})*({v})'

    return (tru, trv)