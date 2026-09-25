__all__ = ['scale_guess']

import sys
import numpy as np

# function to scale parameters from a fit template to the data
def scale_guess(x, y, param, n_gauss, n_poly):
    """Scale inital guess of multigaussian model parameters to data values

    Parameters
    ----------
    x : array_like
        Independent variable values. For EIS data, this will usually
        correspond to wavelength values.
    y : array_like
        Observed data values. ForEIS data, this will be either raw counts
        or calibrated intensity measurements
    param : array_like
        Model fit parameters. There must be 3*n_gauss + n_poly param values.
        For each Gaussian component, the parameters are assumed to have the
        order of [peak, centroid, width]. Polynomial background terms (if any)
        should be at the end and in INCREASING order (e.g. c0, c1, c2, etc.)
    n_gauss : int, optional
        Number of Gaussian components. Default is "1"
    n_poly : int, optional
        Number of background polynomial terms. Common values are:
        0 (no background), 1 (constant), and 2 (linear). Default is "0"

    Returns
    -------
    newparam : array_like
        Array of scaled model parameters.
    """

    # Check inputs
    n_param = len(param)
    if n_param != 3*n_gauss+n_poly:
        print(' ! input parameter sizes do not match ... stopping')
        sys.exit()

    # Copy the input data
    newparam = param.copy()

    # Get background from data (mean of 3 lowest values)
    bkg_data = np.mean(np.sort(y)[0:3])

    # Get background from guess
    # TO-DO: check scaling of higher order terms
    if n_poly > 0:
        bkg_p = param[3*n_gauss::]
        # bkg_guess = np.polyval(bkg_p,x) #outdated API. HIGHEST order first
        bkg_guess = np.polynomial.polynomial.polyval(x, bkg_p) # LOWEST order first
        # scale background
        scale = bkg_data/np.median(np.sort(bkg_guess)[0:3])
        newparam[3*n_gauss::] = bkg_p*scale
        # compute new background
        bkg_p = newparam[3*n_gauss::]
        # new_bkg = np.polyval(bkg_p,x)
        new_bkg = np.polynomial.polynomial.polyval(x, bkg_p)
    else:
        new_bkg = np.zeros(len(x))

    # Compute new peaks
    for n in range(n_gauss):
        gauss_p = param[3*n:3*n+3]
        peak = gauss_p[0]
        cent = gauss_p[1]
        indx = np.abs(x-cent).argmin()
        new_peak = y[indx] - new_bkg[indx]
        if new_peak < 0: 
            new_peak = 0.0
        newparam[3*n] = new_peak

    return newparam
