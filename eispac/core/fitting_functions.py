__all__ = ['multigaussian', 'multigaussian_deviates']

import sys
import numpy as np

# -----------------------------------------------------------------
# model for fitting EIS line profiles using mpfit.py
# -----------------------------------------------------------------

# series of Guassian model functions with polynomial background
def multigaussian(param, x, n_gauss=1, n_poly=0):
    """Multigaussian model with a polynomial background

    Parameters
    ----------
    param : array_like
        Model fit parameters. There must be 3*n_gauss + n_poly param values.
        For each Gaussian component, the parameters are assumed to have the
        order of [peak, centroid, width] for each Gaussian, followed
        by the polynomial terms (if any) in INCREASING order (e.g. c0, c1, c2)
    x : array_like
        Independent variable values to evaluate the function at. For EIS data,
        this will usually correspond to wavelength values.
    n_gauss : int, optional
        Number of Gaussian components. Default is "1"
    n_poly : int, optional
        Number of background polynomial terms. Common values are:
        o (no background), 1 (constant), and 2 (linear). Default is "0"

    Returns
    -------
    f : array_like
        Function evaluated at each x value
    """

    # check inputs
    n_param = len(param)
    if n_param != 3*n_gauss+n_poly:
        print(' ! input parameter sizes do not match ... stopping')
        sys.exit()

    # compute Gaussians
    nx = len(x)
    f  = np.zeros(nx)
    for n in range(n_gauss):
        p = param[3*n:3*n+3]
        arg = ((x-p[1])/p[2])**2
        f = f + p[0]*np.exp(-arg/2.0)

    # compute polynomial background
    if n_poly > 0:
        p = param[3*n_gauss::]
        # f = f + np.polyval(p,x) #outdated: p starts with HIGHEST order term
        f = f + np.polynomial.polynomial.polyval(x, p) #LOWEST order first

    return f

# -----------------------------------------------------------------
# deviates for fitting EIS line profiles using mpfit.py
# -----------------------------------------------------------------

# computes deviates between fit model and data
def multigaussian_deviates(param, x=None, y=None, error=None,
                           n_gauss=None, n_poly=None, fjac=None, debug=False):
    """Computes the deviates between a multigaussian model fit and input data.

    Parameters
    ----------
    param : array_like
        Model fit parameters. There must be 3*n_gauss + n_poly param values.
        For each Gaussian component, the parameters are assumed to have the
        order of [peak, centroid, width] for each Gaussian, followed
        by the polynomial terms (if any) in INCREASING order (e.g. c0, c1, c2)
    x : array_like
    x : array_like
        Independent variable values to evaluate the function at. For EIS data,
        this will usually correspond to wavelength values.
    y : array_like
        Measured values. For EIS data, this will usually correspond to intensity
        observations.
    error : array_like
        Error values for each measurment.
    n_gauss : int
        Number of Gaussian components.
    n_poly : int
        Number of background polynomial terms. Common values are:
        o (no background), 1 (constant), and 2 (linear).
    fjac : None
        Used by mpfit. When fjac == None, partial derivatives will NOT be calculated.
        This is the default for mpfit.
    debug : boolean, optional
        Toggles 'debugging mode'. If set to 'True', then the code will print an
        error statement and exit when it encounters invalid inputs or empty data.
        Default it 'False', which will result in just a negative status flag.


    Returns
    -------
    [status, deviates]
    """

    # check inputs
    n_param = len(param)
    if n_param != 3*n_gauss+n_poly and debug == True:
        print(' ! input parameter sizes do not match ... stopping')
        sys.exit()
    elif n_param != 3*n_gauss+n_poly:
        deviates = np.zeros(len(y)) - 1
        return [-3, deviates]

    # check data
    match, = np.where(error < 0.0)
    nbad  = len(match)
    ndata = len(y)
    if nbad == ndata and debug == True:
        print(' ! no good data ... stopping')
        sys.exit()
    elif nbad == ndata:
        deviates = np.zeros(ndata) - 1
        return [-2, deviates]

    # compute model function
    model = multigaussian(param, x, n_gauss, n_poly)

    # compute deviates
    if nbad == 0:
        deviates = (y-model)/error
    else:
        error[match] = 1.0
        deviates = (y-model)/error
        deviates[match] = 0.0

    status = 0

    return [status, deviates]
