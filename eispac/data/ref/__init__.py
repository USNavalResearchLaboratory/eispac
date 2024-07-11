"""
Template files for fitting EIS spectra
"""
import os
import pathlib
import numpy as np

import eispac

__all__ = ['load_chianti_lines']    

def load_chianti_lines(sort='wave'):
    """Load reference list of spectral lines from CHIANTI

    Parameters
    ----------
    sort : str, optional
        Array column to sort by. Choose from "wave" or "id". Default is "Wave".

    Returns
    -------
    line_arr : numpy recarray
        Array containing spectral lines in both EIS wavebands as computed using
        CHIANTI (version 8?). The array has two columns, "wave" which contains 
        the rest wavelength of the line (in units of [angstrom]) and "id"
        which gives the spectral line ID (e.g. Fe XII). 
    """

    ref_dir = pathlib.Path(os.path.dirname(eispac.__file__)) / "data" / "ref"
    line_id_filepath = ref_dir / 'eis_chianti_lookup_table.txt'
    line_arr = np.genfromtxt(line_id_filepath, delimiter='  ', 
                            dtype=['float', 'U8'], names=['wave', 'id'])
    if sort.lower().startswith('wave'):
        line_arr = np.sort(line_arr, order=['wave'])
    elif sort.lower().startswith('id'):
        line_arr = np.sort(line_arr, order=['id'])

    return line_arr