__all__ = ['export_fits']

import sys
import copy
import shutil
import pathlib
from datetime import datetime, timedelta
import numpy as np
import h5py
from astropy.io import fits
import sunpy.coordinates as coords
from eispac.core.eisfitresult import EISFitResult
from eispac.core.save_fit import lineid_to_name

# function to save fit inten to fits files
def export_fits(fit_result, save_dir=None, verbose=False):
    """Save fit line intensites, velocites, and widths to a .fits file

    Parameters
    ----------
    fit_result : `~eispac.core.EISFitResult` object
        Fit parameter results from eispac.fit_spectra()
    save_dir : str or `pathlib.Path` object, optional
        Directory where the fit results should be saved. If set to None, the
        results will be saved in the same folder as the source data. If set to,
        'cwd', the results will be saved to the current working directory.
        Default is None.
    verbose : bool, optional
        If set to True, will print additional information to the console.
        Default is False.

    Returns
    -------
    output_filepath : list of `pathlib.Path` objects
        Path object pointing to the output files. Each parameter is saved in a
        seperate file and, if there are more than one spectral lines fit, in the
        template, each line will be saved to its own set of files.
    """

    # Validate inputs
    if not isinstance(fit_result, EISFitResult):
        print("Error: Please input a valid EISFitResult object.",
              file=sys.stderr)
        return None
    elif fit_result.meta is None or 'mod_index' not in fit_result.meta.keys():
        print("Error: Missing mod_index containing pointing information.",
              file=sys.stderr)
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
             +"save_dir='cwd' to save files to the current directory",
             file=sys.stderr)
        return None

    if not output_dir.is_dir():
        print("Error: sav_dir is invalid or missing!", file=sys.stderr)
        return None

    # get file name string
    file_prefix = data_name.split('.')[0] # e.g. "eis_YYYYMMDD_HHMMSS"
    line_name = lineid_to_name(fit_result.fit['line_ids'][0])
    try:
        tmplt_id = pathlib.Path(fit_result.meta['filename_template'])
        tmplt_id = tmplt_id.split('.')[1]
    except:
        tmplt_id = str(len(fit_result.fit['line_ids']))+'c'
    output_name = file_prefix+'.'+line_name+'.'+tmplt_id+'-0.inten.fits'
    output_filepath = output_dir.joinpath(output_name)
    print('Saving fit EIS intensities, velocities, and widths to fits files...')
    print('   Directory: '+str(output_dir))

    # Fetch index from the meta structure, cut out spectral data and update
    hdr_dict = copy.deepcopy(fit_result.meta['mod_index'])
    for KEY in ['cname3', 'crval3', 'crpix3', 
                'cdelt3', 'ctype3', 'cunit3', 'naxis3']:
        void = hdr_dict.pop(KEY, None)
    hdr_dict['naxis'] = 2
    code_ver = fit_result.eispac_version
    date_fit = fit_result.date_fit
    hdr_dict['history'] = 'fit using eispac '+code_ver+' on '+date_fit
    data_hdr = fits.Header(hdr_dict)

    # Create the header for the bintable containing the errors
    err_hdr = copy.deepcopy(data_hdr)
    void = err_hdr.pop('bunit', None)
    void = err_hdr.pop('measurement', None)

    total_num_vals = int(data_hdr['naxis2'])
    err_tform = str(total_num_vals)+'D'
    err_tdim = '('+str(data_hdr['naxis1'])+')'
    err_hdr['naxis1'] = total_num_vals*8
    err_hdr['naxis2'] = data_hdr['naxis1']
    err_hdr['tfields'] = data_hdr['naxis2']
    err_hdr['ttype1'] = 'errors'
    err_hdr['tdim1'] = err_tdim
    err_hdr['tform1'] = err_tform

    # Create the bintable table and header for step_date_obs and step_exptime 
    col_date_obs = fits.Column(name='step_date_obs', format='24A', coord_type='UTC',
                              array=fit_result.meta['date_obs'].astype('<U24'))
    col_exptime = fits.Column(name='step_exptime', format='E', unit='s',
                              array=fit_result.meta['duration'])
    extra_hdu = fits.BinTableHDU.from_columns([col_date_obs, col_exptime])

    # Add duplicate timestamps (to conform with the FITS-4 standard)
    hdr_dict['date-obs'] = hdr_dict['date_obs']
    hdr_dict['date-beg'] = hdr_dict['date_beg']
    hdr_dict['date-avg'] = hdr_dict['date_avg']
    hdr_dict['date-end'] = hdr_dict['date_end']
    data_hdr = fits.Header(hdr_dict) # update primary header ONLY

    # Loop over all of the parameters and create the fits files
    # If there are multiple line_ids, output each line to its own set of files
    # and return a list of filepaths (consistent with old IDL workflow)
    params = ['int', 'vel', 'wid']
    param_full_names = ['intensity', 'velocity', 'width']
    num_params = len(params)
    num_line_ids = len(fit_result.fit['line_ids'])
    output_files = []

    for i in range(num_line_ids):
        line_id = fit_result.fit['line_ids'][i]
        if 'NO' in line_id:
            print(f' LINE ID = {line_id}, skipping')
            continue
        line_name = lineid_to_name(line_id)

        for p in range(num_params):
            output_name = file_prefix+'.'+line_name+'.'+tmplt_id+'-'+str(i)+'.'+params[p]+'.fits'
            output_files.append(output_dir.joinpath(output_name))
            if i==0 and p == 0:
                print('   Filenames: '+output_name)
            else:
                print('              '+output_name) # 14 spaces for alignment

            # Fetch data arrays
            if param_full_names[p] == 'intensity':
                data_array = fit_result.fit['int'][:,:,i]
                err_array = fit_result.fit['err_int'][:,:,i]
                data_hdr['bunit'] = fit_result.fit['param_units'][0]
            elif param_full_names[p] == 'velocity':
                data_array = fit_result.fit['vel'][:,:,i]
                err_array = fit_result.fit['err_vel'][:,:,i]
                data_hdr['bunit'] = 'km/s'
            elif param_full_names[p] == 'width':
                data_array = fit_result.fit['params'][:,:,2+3*i]
                err_array = fit_result.fit['perror'][:,:,2+3*i]
                data_hdr['bunit'] = 'Angstrom'

            # Update header information
            data_hdr['line_id'] = line_id
            data_hdr['measrmnt'] = param_full_names[p]
            err_hdr['tunit1'] = data_hdr['bunit']

            # Create fits HDUs and save to a file
            # NB: Errors are stored in a binary table with n_pxls total rows,
            #     each with n_steps values. This is reloaded into a np.recarray
            main_hdu = fits.PrimaryHDU(data_array, header=data_hdr)
            err_col = fits.Column(name='errors', format=err_tform, dim=err_tdim,
                                  unit=err_hdr['tunit1'], array=err_array)
            err_hdu = fits.BinTableHDU.from_columns([err_col], header=err_hdr)
            hdu_list = fits.HDUList([main_hdu, err_hdu, extra_hdu])
            # hdu_list.writeto(output_files[-1], output_verify='silentfix',
            hdu_list.writeto(output_files[-1], output_verify='fix',
                             overwrite=True)

    return output_files
