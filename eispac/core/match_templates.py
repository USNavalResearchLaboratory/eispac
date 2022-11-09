__all__ = ['match_templates']

import os
import sys
import glob
import pathlib
import numpy as np
import h5py
from eispac.core.eiscube import EISCube
from eispac.core.read_wininfo import read_wininfo
import eispac.data

def match_templates(eis_obs):
    """Generate a list of all template files that match an EIS file or window.

    Parameters
    ----------
    eis_obs : `~eispac.core.EISCube` object, str, or `pathlib.Path`
        EIS data to use for searching. If given an EISCube, will only find
        templates that match the selected spectral window. If given the filepath
        to an EIS level-1 HDF5 file, will find all templates that match each
        window in the observation.

    Returns
    -------
    matched_templates : list or list of lists
        List of template files that match the selected spectral window. If a
        filepath was input, a list of lists will be returned instead. Each
        sublist contains the matched templates for the corresponding window
        (e.g. matched_templates[0] would contain the list of templates matching
        the first window in the data file).
    """
    # Validate input
    if isinstance(eis_obs, (str, pathlib.Path)):
        wininfo = read_wininfo(eis_obs)
        iwin_list = list(range(wininfo.size))
    elif isinstance(eis_obs, EISCube):
        wininfo = eis_obs.meta['wininfo']
        iwin_list = [eis_obs.meta['iwin']]
    else:
        print('Error: missing or invalid eis_obs. Please input a filepath or'
             +' EISCube class instance', file=sys.stderr)
        return None

    # Get master list of all templates included in eispac
    all_templates = np.array(eispac.data.fit_template_filenames())
    num_templates = len(all_templates)
    template_waves = np.zeros(num_templates)

    # Extract line wavelengths
    # Note: filenames have forms like "fe_12_195_119.2c.template.h5"
    #       or "s__11_188_675.3c.template.h5"
    for t in range(num_templates):
        line_id_str = all_templates[t].name.split('.')[0]
        line_id_str.replace('__', '_') # single character elements are padded
        line_id_parts = line_id_str.split('_')
        template_waves[t] = float(line_id_parts[2]+'.'+line_id_parts[3])

    # Loop over each spectral window and find all matching templates
    matched_templates = []
    for w in iwin_list:
        temp_file_list = []
        loc_templates = np.where((template_waves > wininfo[w].wvl_min)
                                &(template_waves < wininfo[w].wvl_max))[0]
        if len(loc_templates) > 0:
            for f in loc_templates:
                temp_file_list.append(all_templates[f])
        matched_templates.append(temp_file_list)

    if len(iwin_list) == 1:
        return matched_templates[0]
    else:
        return matched_templates
