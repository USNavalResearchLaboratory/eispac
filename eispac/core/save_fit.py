__all__ = ['lineid_to_name', 'save_fit']

import sys
import shutil
import pathlib
import numpy as np
import h5py
from eispac.core.eisfitresult import EISFitResult

# try:
#     # Import local versions of submodules
#     from eisfitresult import EISFitResult
#     print('Notice: Loading local version of eispac submodules')
# except:
#     # Import from installed package
#     from eispac.eisfitresult import EISFitResult

# function to convert line ids to string name
def lineid_to_name(lineid, component=None):

    # Decode byte strings (if needed)
    if not isinstance(lineid, str):
        lineid = lineid.decode('utf-8')

    # get components of string
    split = lineid.split(' ')
    element = split[0].lower()
    ion = split[1]
    wave = split[2].split('.')
    # fix element symbol
    if len(element) == 1: element = element[0]+'_'
    # convert roman numeral to integer
    roman = {'I':1, 'V':5, 'X':10, 'L':50}
    nion = len(ion)
    numbr = 0
    for i in range(nion):
        if (i+1) == nion or roman[ion[i]] >= roman[ion[i+1]]:
            numbr += roman[ion[i]]
        else:
            numbr -= roman[ion[i]]
    # assemble string name
    name = element+'_'+str(numbr)+'_'+str(wave[0])+'_'+str(wave[1])
    # if component is not None: name += '_c'+str(component)
    return name

def walk_and_save(hdf5_group, data_obj, data_key, print_str, verbose=False):
    """Helper function for saving EISFitResult data in a HDF5 file.

    Walks through dictionaries or lists and recursively saves data it
    finds along the way.

    Parameters
    ----------
    hdf5_group : h5py.File or h5py.Group objects
        HDF5 file or group object currently open in write mode
    data_obj : dict, list, numpy array, or scalar
        Data or container with data to be saved
    verbose : bool, optional
        If set to True, will print the name of each data variable saved.
        Default is False.

    Returns
    -------
    None
    """

    if isinstance(data_obj, (dict, list)):
        # Create a subgroup and recursively loop over each subgroup key
        subgrp = hdf5_group.create_group(str(data_key))
        if isinstance(data_obj, dict):
            subgrp_index = data_obj.keys()
        else:
            # For lists, the keys are just index numbers (e.g. 0, 1, 2,...)
            subgrp_index = range(len(data_obj))
        for subkey in subgrp_index:
            new_print_str = print_str+'/'+str(subkey)
            walk_and_save(subgrp, data_obj[subkey], subkey,
                          new_print_str, verbose=verbose)
    else:
        # Save data to group (also converts numpy fixed-width
        # Note: numpy arrays of fixed-width unicode strings will be
        #       converted to ascii, for better storage in HDF5 files
        if verbose:
            print('   ', print_str)
        temp_data = data_obj
        if isinstance(data_obj, (np.ndarray, np.str_)):
            if str(data_obj.dtype).startswith('U', 1):
                temp_data = data_obj.astype(np.string_)
        hdf5_group.create_dataset(str(data_key), data=temp_data)

# function to save fit dictionary
def save_fit(fit_result, save_dir=None, verbose=False):
    """Save an EISFitResult object to an HDF5 file (or files)

    Parameters
    ----------
    fit_result : EISFitResult object
        Fit parameter results from eispac.fit_spectra()
    save_dir : str or pathlib.Path object, optional
        Directory where the fit results should be saved. If set to None, the
        results will be saved in the same folder as the source data. If set to,
        'cwd', the results will be saved to the current working directory.
        Default is None.
    verbose : bool, optional
        If set to True, will print the name of each data variable saved.
        Default is False.

    Returns
    -------
    output_filepath : list or pathlib.Path object
        Path object pointing to the output file. If there are more one than one
        spectral lines fit, multiple copies of the same output file will be saved
        but with different filenames.
    """

    # Validate inputs
    if not isinstance(fit_result, EISFitResult):
        print("Error: Please input a valid EISFitResult object.", file=sys.stderr)
        return None

    # Parse filename and determine the directory and filename
    data_filepath = pathlib.Path(fit_result.meta['filename_head']).resolve()
    data_name = str(data_filepath.name)
    if save_dir is None:
        output_dir = data_filepath.parent
    elif isinstance(save_dir, pathlib.Path):
        output_dir = save_dir
    elif isinstance(save_dir, str):
        if save_dir.casefold() == 'cwd':
            output_dir = pathlib.Path().cwd()
        else:
            output_dir = pathlib.Path(save_dir)
    else:
        print("Error: Please input a valid save directory or set "
             +"save_dir='cwd' to save files to the current directory", file=sys.stderr)
        return None

    if not output_dir.is_dir():
        print("Error: sav_dir is invalid or missing!", file=sys.stderr)
        return None

    # get file name string
    file_prefix = data_name.split('.')[0]
    line_name = lineid_to_name(fit_result.fit['line_ids'][0])
    try:
        template_id = pathlib.Path(fit_result.meta['filename_template'])
        template_id = template_id.split('.')[1]
    except:
        template_id = str(len(fit_result.fit['line_ids']))+'c'
    output_name = file_prefix+'.'+line_name+'.'+template_id+'-0.fit.h5'
    output_filepath = output_dir.joinpath(output_name)
    print('Saving EIS fit results...')
    print('   Directory: '+str(output_dir))
    print('   Filenames: '+output_name)

    # Convert the .meta['wininfo'] record array to a dict of 1D arrays
    # Note, this is needed because HDF5 does not play nice with recarrays
    if 'wininfo' in fit_result.meta.keys():
        wi_rec = fit_result.meta['wininfo'].copy()
        wininfo_dict = {name:wi_rec[name] for name in wi_rec.dtype.names}
        fit_result.meta['wininfo'] = wininfo_dict

    # Loop over each EISFitResult attribute and save the structure as found
    # TODO: consider setting dtype directly (w/ h5py special types for strings)
    real_main_comp = fit_result.fit['main_component']
    fit_result.fit['main_component'] = 0
    with h5py.File(output_filepath, 'w') as fit_file:
        for attr_name in vars(fit_result):
            if attr_name == 'fit_func':
                # Don't try saving python functions, it does not work.
                continue
            attr = getattr(fit_result, attr_name)
            walk_and_save(fit_file, attr, attr_name, attr_name, verbose=verbose)

    # reset real main component number
    fit_result.fit['main_component'] = real_main_comp

    # Ensure wininfo stays as a recarray in the current session
    if 'wininfo' in fit_result.meta.keys():
        fit_result.meta['wininfo'] = wi_rec

    # If there are multiple line_ids, make a copy of the file with a new name
    # and return a list of filepaths (consistent with old IDL workflow)
    num_line_ids = len(fit_result.fit['line_ids'])
    if num_line_ids > 1:
        list_output = [output_filepath]
        for line_num in range(1, num_line_ids):
            new_line_name = lineid_to_name(fit_result.fit['line_ids'][line_num])
            new_out_name = file_prefix+'.'+new_line_name+'.'+template_id+'-'+str(line_num)+'.fit.h5'
            list_output.append(output_dir.joinpath(new_out_name))
            print('              '+new_out_name) # 14 spaces to align filenames
            shutil.copy(output_filepath, list_output[line_num])
            with h5py.File(list_output[line_num], 'r+') as new_file:
                main_comp = new_file['fit/main_component']
                main_comp[...] = line_num
        return list_output
    else:
        # If only one line, just return the filepath directly (no list)
        return output_filepath

if __name__ == '__main__':

    fake_fit = EISFitResult(empty=True)
    fake_meta = {'filename_head':'eis_fake_fit.h5'}
    fake_fit.meta = fake_meta
    fake_fit.fit['line_ids'] = np.array(['Fe XII 195.119', 'Fe XII 195.179'])

    fit_filepath = save_fit(fake_fit, save_dir='cwd')
    print(fit_filepath)
