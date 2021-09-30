__all__ = ['calc_velocity']

import sys
import numpy as np
import astropy.constants as const

def calc_velocity(observed_wave, rest_wave):
    """Calculate the Doppler velocity of a given line profile

    Parameters
    ----------
    observed_wave : list, tuple, or numpy array
        Observed centroid wavelength values for a given spectral line in units
        of [Angstrom]
    rest_wave : float or str
        Rest wavelength value in units of [Angstrom]. If given a string, will
        check to see if it is a spectral line ID and, if it is, will try to
        extract the rest wavelength.

    Returns
    -------
    velocity : numpy array
        Doppler velocity in units of [km/s]. Computed using the equation,
            vel = c*(obs_wave - rest_wave)/rest_wave
    """
    # Check input values
    if isinstance(observed_wave, (list, tuple)):
        obs_wave = np.array(observed_wave)
    elif isinstance(observed_wave, np.ndarray):
        obs_wave = observed_wave
    else:
        print('Error: invalid observed_wave datatype. Please input a list or'
             +' array with valid wavelength values.', file=sys.stderr)
        return 0.0

    if isinstance(rest_wave, (int, float)):
        wave_0 = float(rest_wave)
    elif isinstance(rest_wave, str):
        if rest_wave.isnumeric():
            wave_0 = float(rest_wave)
        else:
            # Assume an input line ID (e.g. "Fe XII 195.119")
            try:
                wave_0 = float(rest_wave.split()[2])
            except:
                print('Error: invalid line ID input to rest_wave. Please use a'
                     +' line ID with a format like "Fe XII 195.119".',
                     file=sys.stderr)
                return 0.0
    else:
        print('Error: invalid rest_wave datatype. Please input a float or string'
             +' value with the line rest wavelength.', file=sys.stderr)
        return 0.0

    velocity = const.c.to('km/s').value*(obs_wave-wave_0)/wave_0
    median_vel = np.median(velocity)
    velocity = velocity - median_vel # Rough correction for S/C temp.

    return velocity
