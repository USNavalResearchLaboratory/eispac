
__all__ = ['calc_intensity_range']

import numpy as np

def calc_intensity_range(intensity, pmin=1, pmax=99, lower=1.0E-2, imin=10):
    """
    Compute an intensity range for an image based on the histogram. Uses np.percentile.

    Inputs:
      intensity -> the intensity array
      pmin -> lower cutoff
      pmax -> upper cutoff
      lower -> imin = lower*imax
      imin -> specify imin

    Outputs:
      irange -> two-element list of [imin, imax]

    """

    if pmax > 100: pmax = 100
    if pmin < 0: pmax = 0

    irange = np.percentile(intensity, (pmin,pmax))
    irange = [irange[1]*lower, irange[1]]
    
    if irange[0] < imin: irange[0] = imin

    return irange


    
