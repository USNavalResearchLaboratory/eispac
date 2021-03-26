__all__ = ['read_fit']

import os
import sys
import pathlib
import numpy as np
import h5py
import eispac.core.fitting_functions as fit_fns
from eispac.core.eisfitresult import EISFitResult

# if __name__ == '__main__':
#     # Import local versions of submodules
#     print('Notice: Loading local version of eispac submodules')
#     import fitting_functions as fit_fns
#     from eisfitresult import EISFitResult
# else:
#     # Import from installed package
#     import eispac.fitting_functions as fit_fns
#     from eispac.eisfitresult import EISFitResult

def walk_and_load(hdf5_file, hdf5_path, verbose=False):
    """Helper function for loading EISFitResult data from a HDF5 file.

    Walks the internal structure of an HDF5 file and recursively loads data it
    finds along the way.

    Parameters
    ----------
    hdf5_file : h5py.File object
        HDF5 file object currently open in read mode
    hdf5_path : str
        String giving the internal path to a HDF5 group or dataset
    verbose : bool, optional
        If set to True, will print the name of each data variable it read.
        Default is False.

    Returns
    -------
    output : dict, list, scalar, or numpy array
        Loaded data values. In the case of an HDF5 group, will return a dict
        (with possible nested dicts) with the discovered values.
    """
    if isinstance(hdf5_file[hdf5_path], h5py.Group):
        # Create a dict (or list) and recursively loop over each item
        grp_keys = hdf5_file[hdf5_path].keys()
        # if hdf5_path == 'parinfo' or hdf5_path == 'funcinfo':
        #     # Special exception for restoring the 'parinfo' & 'funcinfo' lists
        if all([key.isdigit() for key in grp_keys]):
            # Restoring lists such as parinfo and funcinfo
            grp_list = list()
            for key in grp_keys:
                new_path = hdf5_path+'/'+key
                sub_val = walk_and_load(hdf5_file, new_path, verbose=verbose)
                grp_list.append(sub_val)
            output = grp_list
        else:
            grp_dict = {}
            for key in grp_keys:
                new_path = hdf5_path+'/'+key
                sub_val = walk_and_load(hdf5_file, new_path, verbose=verbose)
                grp_dict[key] = sub_val
            output = grp_dict
    elif isinstance(hdf5_file[hdf5_path], h5py.Dataset):
        # Read dataset from the file
        output = np.array(hdf5_file[hdf5_path])
        if verbose:
            print('   ', hdf5_path)
        if (isinstance(output.dtype, np.bytes_) or str(output.dtype) == 'object'
            or str(output.dtype).startswith('S', 1)):
            # Convert objects, byte or ascii strings to uncode
            output = output.astype(np.unicode_) # strings to unicode
        if output.size == 1:
            # Extract value from a single-element array (0- or 1-D arrays only)
            output = output.item()

    return output

# function to read fit dictionary
def read_fit(filename, verbose=False):
    """Load an EISFitResult object from an HDF5 file

    Parameters
    ----------
    filename : str or pathlib.Path object
        String or path to the fit result file that should be loaded.
    verbose : bool, optional
        If set to True, will print the name of each data variable read in.
        Default is False.

    Returns
    -------
    fit_result : EISFitResult object
        Copy of the fit results loaded from the file.
    """

    # Input type validation (value checks are implemented later)
    if not isinstance(filename, (str, pathlib.Path)):
        print('Error: Please input a valid filepath as '
                +'either a string or pathlib.Path object', file=sys.stderr)
        return None

    # Parse filename and determine the directory and filename
    abs_filepath = pathlib.Path(filename).resolve()
    input_name = str(abs_filepath.name)
    input_dir = abs_filepath.parent
    if str(input_dir) == '.':
        input_dir = pathlib.Path().cwd()

    fit_filepath = input_dir.joinpath(input_name)
    if not fit_filepath.is_file():
        print('Error: fit result file does not exist, ' + str(filename), file=sys.stderr)
        return None

    # Initialize the output EISFitResult object
    fit_result = EISFitResult(empty=True)

    # Loop over each EISFitResult attribute and load data structure from the file
    print('Reading fit result from ', str(filename))
    with h5py.File(filename, 'r') as fit_file:
        top_key_list = list(fit_file.keys())
        fit_key_list = list(fit_file['fit'].keys())
        for attr_name in fit_file.keys():
            attr_val = walk_and_load(fit_file, attr_name, verbose=verbose)
            setattr(fit_result, attr_name, attr_val)

    # Restore the fit function
    fit_result.fit_func = getattr(fit_fns, fit_result.func_name)

    # Make sure the .fit['Line_ids'] is ALWAYS an array (for code consistency)
    fit_result.fit['line_ids'] = np.atleast_1d(fit_result.fit['line_ids'])

    # If version number is missing, try to guess it
    if 'eispac_version' not in top_key_list:
        if 'data_units' in top_key_list:
            fit_result.eispac_version = '0.9.1'
        else:
            fit_result.eispac_version = '0.8.0'

    # Add 'wave_range' to fits saved before 2021-02-19
    if 'wave_range' not in fit_key_list:
        fit_result.fit['wave_range'] = np.zeros(2)
        fit_result.fit['wave_range'][0] = fit_result.template['data_x'][0]
        fit_result.fit['wave_range'][1] = fit_result.template['data_x'][-1]

    # Restore recarray of .meta['wininfo'] (stored in HDF5 as dict of arrays)
    if 'wininfo' in fit_result.meta.keys():
        wi_dict = fit_result.meta['wininfo']
        # key_names = list(wi_dict.keys()) # give alphabetical ordering...
        key_names = ['iwin', 'line_id', 'wvl_min', 'wvl_max', 'nl', 'xs']
        key_dtypes = [wi_dict[key].dtype for key in key_names]
        num_recs = len(wi_dict[key_names[0]])
        rec_list = [tuple([wi_dict[key][num] for key in key_names]) for num in range(num_recs)]
        wininfo_rec = np.rec.fromrecords(rec_list, names=key_names, formats=key_dtypes)
        fit_result.meta['wininfo'] = wininfo_rec

    return fit_result

if __name__ == '__main__':

    filename = './data/test/eis_20190404_131513_fe_12_195_119_1c.fit.h5'
    fit_res = read_fit(filename)

    for key in fit_res.fit.keys():
        if np.size(fit[key]) > 1:
            print('{:12} {:12} {:12}'.format(key,str(fit[key].dtype),str(np.shape(fit[key]))))
        else:
            print('{:12} {:12} {:12}'.format(key,str(fit[key].dtype),str(np.size(fit[key]))))
