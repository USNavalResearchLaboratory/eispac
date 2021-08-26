__all__ = ['export_fits']

import sys
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
    """Save a fit line intensity to a .fits file

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
        spectral lines fit, each line intensity will be saved to its own file.
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
    output_name = file_prefix+'.'+line_name+'.'+template_id+'-0.inten.fits'
    output_filepath = output_dir.joinpath(output_name)
    print('Saving fit EIS intensities to fits files...')
    print('   Directory: '+str(output_dir))
    print('   Filenames: '+output_name)

    # Get the corrected center coords (using both the AIA and CCD offsets)
    nx_steps = fit_result.n_steps
    ny_pxls = fit_result.n_pxls
    pointing = fit_result.meta['pointing']

    x_center = pointing['xcen'] + pointing['offset_x']
    y_center = pointing['ycen'] + pointing['offset_y']
    y_center = y_center - np.mean(fit_result.meta['ccd_offset'])
    x1 = x_center - pointing['x_scale']*nx_steps/2.0
    x2 = x_center + pointing['x_scale']*nx_steps/2.0
    y1 = y_center - pointing['y_scale']*ny_pxls/2.0
    y2 = y_center + pointing['y_scale']*ny_pxls/2.0

    # Determine the proper start and end time for the entire raster
    date_obs = fit_result.meta['index']['date_obs']
    median_exp = np.median(fit_result.meta['duration'])
    date_end = datetime.fromisoformat(date_obs) + timedelta(seconds=nx_steps*median_exp)
    date_end = date_end.isoformat(timespec='milliseconds')

    # Determine the effective location of EIS (i.e. Earth)
    hg_coords = coords.get_body_heliographic_stonyhurst('earth', time=date_obs)

    # Create a clean and updated fits header
    output_hdr = fits.Header()
    output_hdr['naxis'] = 2
    output_hdr['naxis1'] = nx_steps
    output_hdr['naxis2'] = ny_pxls

    output_hdr['date_obs'] = date_obs
    output_hdr['date_end'] = date_end

    output_hdr['crval1'] = x1
    output_hdr['crpix1'] = 0
    output_hdr['cdelt1'] = pointing['x_scale']
    output_hdr['crota1'] = fit_result.meta['index']['crota1']
    output_hdr['ctype1'] = 'HPLN-TAN'
    output_hdr['cunit1'] = 'arcsec'

    output_hdr['crval2'] = y1
    output_hdr['crpix2'] = 0
    output_hdr['cdelt2'] = pointing['y_scale']
    output_hdr['crota2'] = fit_result.meta['index']['crota2']
    output_hdr['ctype2'] = 'HPLT-TAN'
    output_hdr['cunit2'] = 'arcsec'

    output_hdr['line_id'] = fit_result.fit['line_ids'][0]
    # output_hdr['fit_file'] = fit_file
    output_hdr['fovx'] = x2-x1
    output_hdr['fovy'] = ny_pxls
    output_hdr['xcen'] = x1 + 0.5*(x2-x1)
    output_hdr['ycen'] = y1 + 0.5*(y2-y1)

    output_hdr['hgln_obs'] = hg_coords.lon.deg
    output_hdr['hglt_obs'] = hg_coords.lat.deg
    output_hdr['dsun_obs'] = hg_coords.radius.m

    output_hdr['history'] = 'fit using eispac '+fit_result.eispac_version+' on '+fit_result.date_fit


    # Save the first (and maybe only) .fits file
    fits.writeto(output_filepath, fit_result.fit['int'][:,:,0],
                 header=output_hdr, output_verify='fix', overwrite=True)

    # If there are multiple line_ids, output each intensity to its own file
    # and return a list of filepaths (consistent with old IDL workflow)
    num_line_ids = len(fit_result.fit['line_ids'])
    if num_line_ids > 1:
        list_output = [output_filepath]
        for line_num in range(1, num_line_ids):
            line_id = fit_result.fit['line_ids'][line_num]
            if 'NO' in line_id:
                print(f' LINE ID = {line_id}, skipping')
                continue
            new_line_name = lineid_to_name(line_id)
            new_out_name = file_prefix+'.'+new_line_name+'.'+template_id+'-'+str(line_num)+'.inten.fits'
            list_output.append(output_dir.joinpath(new_out_name))
            print('              '+new_out_name) # 14 spaces to align filenames

            # Update line_id in header and save another .fits file
            output_hdr['line_id'] = line_id
            fits.writeto(list_output[-1], fit_result.fit['int'][:,:,line_num],
                         header=output_hdr, output_verify='silentfix', overwrite=True)
        return list_output
    else:
        # If only one line, just return the filepath directly (no list)
        return output_filepath
