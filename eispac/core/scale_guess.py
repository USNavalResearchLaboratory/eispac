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
        following order: [peak, centroid, width]
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

    # check inputs
    n_param = len(param)
    if n_param != 3*n_gauss+n_poly:
        print(' ! input parameter sizes do not match ... stopping')
        sys.exit()

    # copy the input data
    newparam = param.copy()

    # get background from data
    bkg_data = np.mean(np.sort(y)[0:3])

    # get background from guess
    if n_poly > 0:
        bkg = param[3*n_gauss::]
        bkg_guess = np.median(np.sort(np.polyval(bkg,x))[0:3])
        # scale background
        scale = bkg_data/bkg_guess
        newparam[3*n_gauss::] = bkg*scale
        # compute new background
        bkg = newparam[3*n_gauss::]
        new_bkg = np.polyval(bkg,x)
    else:
        new_bkg = np.zeros(len(x))

    # compute new peaks
    for n in range(n_gauss):
        p = param[3*n:3*n+3]
        peak = p[0]
        cent = p[1]
        indx = np.abs(x-cent).argmin()
        new_peak = y[indx]-new_bkg[indx]
        if new_peak < 0: new_peak = 0.0
        newparam[3*n] = new_peak

    return newparam
