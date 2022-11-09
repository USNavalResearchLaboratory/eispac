
__all__ = ['scale_intensity']

import numpy as np

def scale_intensity(intensity, irange, log=True):
    """Scale an array of intensities to have values in the range of 0 to 1.
    Can be useful for plotting.

    Parameters
    ----------
    intensity : array_like
        The intensity array. Can have any dimensions.
    irange : array_like
        Min and max intensity values used for the scaling. Intensity values
        outside this range will be clipped to the nearest bound.
    log : bool, optional
        If set to True, will scale the log of the intensity values. Default is
        True.

    Returns
    -------
    scaled : array_like
        Array of scaled intensity values. Will have the same dimensions as the
        input intensity array. The caling is calculated using the equation
        scaled = (intensity - imin) / (imax - imin), where imin and imax are
        specified by the input "irange" parameter.
    """

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
