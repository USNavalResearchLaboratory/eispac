
__all__ = ['scale_intensity']

import numpy as np

def scale_intensity(intensity, irange, log=True):

    if log:
        scaled = np.clip(intensity, irange[0], irange[1])
        scaled = np.log10(scaled)
        imin = np.log10(irange[0])
        imax = np.log10(irange[1])
        scaled = (scaled-imin)/(imax-imin)
    else:
        scaled = np.clip(intensity, irange[0], irange[1])
        imin = irange[0]
        imax = irange[1]
        scaled = (scaled-imin)/(imax-imin)        
        scaled = np.log10(scaled)

    return scaled
