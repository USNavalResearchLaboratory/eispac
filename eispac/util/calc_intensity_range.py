
__all__ = ['calc_intensity_range']

import numpy as np

def calc_intensity_range(intensity, pmin=1, pmax=99, lower=1.0E-2, imin=10):
    """Compute an intensity range for an image based on the histogram.
    Uses np.percentile.

    Parameters
    ----------
    intensity : array-like
        The intensity array. Can have any dimensions.
    pmin : float, optional
        Lower cutoff percentile. Must be a value between 0 and 100, inclusive.
        Default is 1.
    pmax : float, optional
        Upper cutoff percentile. Must be a value between 0 and 100, inclusive.
        Default is 99.
    lower : float, optional
        Scale factor to calculate the minimum intensity range value using
        imin = lower*imax. Default is 0.01
    imin : float, optional
        Limit on the minimum intensity range value. Default is 10

    Returns
    -------
    irange : list
        two-element list of [imin, imax]
    """

    if pmax > 100: pmax = 100
    if pmin < 0: pmax = 0

    irange = np.percentile(intensity, (pmin,pmax))
    irange = [irange[1]*lower, irange[1]]

    if irange[0] < imin: irange[0] = imin

    return irange
