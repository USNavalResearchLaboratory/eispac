"""
Functions for computing pixel offsets on CCD
"""
__all__ = ['ccd_offset']

import numpy as np

def ccd_offset(wavelength):
    """Calculate the spatial offset of a line relative to He II 256 Å

    Spatial offset of the specified wavelength relative
    to He II 256 Å. If you see a feature in the He II 256 image at
    coordinate :math:`Y`, then the corrected coordinate
    :math:`Y^{\prime}` for any other wavelength :math:`\lambda` is,

    .. math::
        Y^{\prime} = Y - o(\lambda),

    where :math:`o(\lambda)` is the CCD offset.
    Note that the spatial coordinate system for EIS is evalulated for
    He II 256 Å by cross-correlating with SOT.

    The value of ``ccd_offset`` represents the number of pixels
    that a feature seen in a wavelength sits above the He II 256 Å
    image on the CCD.

    The following goes into the calculation of EIS CCD offset:

    - the tilt of the grating is assumed to be linear with wavelength
    - the tilt for SW was derived by [young09]_
    - the tilt for LW was *assumed* to be the same as for SW
    - the offset between SW and LW was measured by co-aligning images
      from Fe VIII 185.21 Å and Si VII 275.35 Å

    .. note:: This routine is a (nearly) verbatim translation of the IDL version of
              eis_ccd_offset written by Peter Young. The above documentation has
              been adapted from that routine.

    Parameters
    ----------
    wavelength : array_like
        Wavelength of the spectral line(s) of interest in units of [Å]

    Returns
    -------
    offset : `~numpy.ndarray`
        The spatial offset between the specified wavelength and He II 256.32.
        The value represents how many pixels above He II 256.32 the specified
        wavelength sits on the EIS detectors.

    References
    ----------
    .. [young09] Young, P.R., et al., 2009, A&A, 495, `587 <https://doi.org/10.1051/0004-6361:200810143>`_
    """
    # Calculate the pixel offset using an equation based on the wavelengths of
    # specific reference lines (Fe VIII 185.21, Si 275.35, and He II 256.32)

    # Note_1: The y-scale of EIS is always given as 1 [arcsec]/[pixel]. Since
    #         this offset is an instrumental effect specific to EIS, there is
    #         not much point in generalizing this function for use with other
    #         y-scale values.
    # Note_2: The value of 18.5 accounts for the base offset between the
    #         shortwave and longwave CCDs.
    grating_tilt = -0.0792  # [pixels]/[angstrom]
    wavelength = np.atleast_1d(wavelength)
    # Calculate the offset for all lines in the shortwave band
    offset = grating_tilt * (wavelength-185.21) + 18.5 + grating_tilt * (275.35 - 256.32)
    # Find and calculate the offset for all lines in the longwave band
    i_ = np.where(wavelength > 230)
    offset[i_] = grating_tilt * (wavelength[i_] - 256.32)

    return offset
