__all__ = ['read_template']

import sys
import pathlib
import warnings

import numpy as np
import h5py

from eispac.core.eisfittemplate import EISFitTemplate

def read_template(filename):
    """Load an `EISFitTemplate` from an HDF5 or TOML template file

    Parameters
    ----------
    filename : str or `pathlib.Path`
        Path to a HDF5 template file provided with eispac or a 
        user-made TOML template file.

    Returns
    -------
    cls : `EISFitTemplate` class instance
        Object containing the fit template
    """
    # Note: this is just a convience function to avoid breaking the old API
    #       Please see the code in the class method
    return EISFitTemplate.read_template(filename)
    
    # # TODO: Add loading custom template from a .toml file

    # # NOTE: return None here rather than allow h5py to handle
    # # exception so that spectral fitting pipeline can error
    # # more gracefully
    # # FIXME: replace with proper exception handling and logging
    # if not isinstance(filename, (str, pathlib.Path)):
    #     warnings.warn('Error: Template filepath must be either a string or pathlib.Path')
    #     return None
    
    # filename = pathlib.Path(filename)
    # if not filename.is_file():
    #     warnings.warn(f'Error: Template filepath {filename} does not exist')
    #     return None

    # file_type = filename.suffix

    # if file_type.lower() in ['.h5', '.hdf5']:
    #     # Load a standard format HDF5 file (probably packaged with eispac)
    #     with h5py.File(filename, 'r') as f_temp:
    #         # Template
    #         template = {}
    #         for key in f_temp['template']:
    #             val = f_temp['template/'+key]
    #             if key == 'line_ids':
    #                 val = np.char.decode(val).flatten() # convert bytes to unicode
    #             elif len(val) > 1:
    #                 val = np.array(val)
    #             else:
    #                 val = val[0]
    #             template[key] = val
    #         # Parinfo
    #         nstr = len(f_temp['parinfo/value'])
    #         parinfo = []
    #         for istr in range(nstr):
    #             parameter = {}
    #             for key in f_temp['parinfo']:
    #                 val = f_temp['parinfo/'+key][istr]
    #                 if key == 'tied':
    #                     val = np.char.decode(val) # convert bytes to unicode
    #                 parameter[key] = val
    #             parinfo.append(parameter)
    # elif file_type.lower() == '.toml':
    #     # Load a custom template stored in a TOML file
    #     with open(filename, 'rb') as f_temp:
    #         toml_dict = tomllib.load(f_temp)
    #         # Ensure top-level keys are all lower case
    #         toml_dict = {KEY.lower(): VALUE for KEY, VALUE in toml_dict.items()}

    #         template = toml_dict.get('template', None)
    #         parinfo = toml_dict.get('parinfo', None)
    #         # note: this parinfo will be an alternative dict of lists

    # return EISFitTemplate(filename, template, parinfo)