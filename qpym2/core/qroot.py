"""
QRoot-layer utilities: ROOT filter strings and column definitions for staging.
"""


def get_staging_defs(evar, esum_var, ch1, ch2, ds, t):
    """Return the standard list of (alias, expression) defines for staging.

    Args:
        evar (str): Branch name for the higher-energy detector energy.
        esum_var (str): Branch name for the summed energy.
        ch1, ch2, ds, t (str): Branch names for channel, dataset, time.

    Returns:
        list of (str, str): Alias/expression pairs for RDataFrame.Define.
    """
    defs = [
        ('e1',   evar),
        ('esum', esum_var),
        ('ch1',  ch1),
        ('ch2',  ch2),
        ('ds',   ds),
        ('t',    t),
    ]
    defs.append(('u0', f'{esum_var}/sqrt(2)'))
    defs.append(('v0', f'(2*{evar} - {esum_var})/sqrt(2)'))
    return defs


def get_roi_cuts(esum_var='TotalEnergy', evar='Energy',
                 esum_min=2390, esum_max=2710, threshold=340,
                 shift_padding=0):
    """Return a ROOT filter string selecting events within the analysis ROI.

    Boundaries are loosened by ``shift_padding`` on each side to ensure events
    that fall within the ROI after an energy-scale shift are not discarded at
    staging time.

    Applies:
      - total energy window: ``esum_min - padding < esum_var < esum_max + padding``
      - per-detector threshold: ``evar > esum_var/2`` (i.e. e1 >= e2)
        and ``esum_var - evar > threshold - padding`` (i.e. e2 above threshold)

    Args:
        esum_var (str): Branch name for the summed energy. Default ``'TotalEnergy'``.
        evar (str): Branch name for the higher-energy detector. Default ``'Energy'``.
        esum_min (float): Lower bound on summed energy (keV).
        esum_max (float): Upper bound on summed energy (keV).
        threshold (float): Per-detector energy threshold (keV).
        shift_padding (float): Extra margin (keV) added to all boundaries.

    Returns:
        str: ROOT filter expression.
    """
    cuts = f"{esum_var}>{esum_min - shift_padding} && {esum_var}<{esum_max + shift_padding}"
    cuts += f" && {evar}>{esum_var}/2 && {esum_var}-{evar}>{threshold - shift_padding}"
    return cuts
