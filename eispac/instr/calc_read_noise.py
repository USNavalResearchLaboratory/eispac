__all__ = ['calc_read_noise']

import sys
import numpy as np

def calc_read_noise(input_wave):
    """Calculate the read noise counts for a given wavelength spectrum

    Parameters
    ----------
    input_wave : float, list, tuple, or `~numpy.ndarray`
        Wavelength values to compute the read noise at

    Returns
    -------
    read_noise_counts : float or `~numpy.ndarray`
        read noise in units of [photon counts]. Will have the same dimensions as
        the input_wave array.
    """
    # Check input values
    if isinstance(input_wave, (list, tuple)):
        use_wave = np.array(input_wave)
    elif isinstance(input_wave, np.ndarray):
        use_wave = input_wave
    else:
        print('Error: invalid input_wave datatype. Please input a list or'
             +' array with valid wavelength values', file=sys.stderr)
        return 0.0


    # (1) read noise [e] = read noise [DN] x gain = 2.29 DN * 6.3 e/DN = 14.427e
    # Note: The value of 2.29 was determined by measuring the dispersion in the
    #       signal with the shutter closed.
    read_noise_per_e = 14.427 # units of [electrons]

    # (2) Electrons per photon = (12398.5 [AA/eV*ph] / lambda [AA]) / (3.65 [eV/e])
    # i.e. One electron is produced for every 3.65 eV in incident energy.
    e_per_ph = (12398.5/use_wave)/(3.65) # units of [electrons/photon]

    # (3) read noise [ph] = read noise [e] / electrons per photon
    read_noise_counts = read_noise_per_e/e_per_ph

    return read_noise_counts
