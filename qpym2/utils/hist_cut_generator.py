""" This file contains utility functions to generate cuts based on the 
histogram model parameters
"""
def for_qroot(esum_var='TotalEnergy', evar='Energy',
              esum_min=2400, esum_max=2570, threshold=200, shift_padding=10):
    """ Return the hist cuts string for the given parameters.
    shift_padding: extra keV below esum_min to include in staging, so events
    that shift up into the ROI after energy-scale correction are not lost.
    """
    hist_cuts = f"{esum_var}>{esum_min - shift_padding} && {esum_var}<{esum_max}"
    hist_cuts += f" && {evar}>{esum_var}/2 && {esum_var}-{evar} > {threshold - shift_padding}"
    return hist_cuts