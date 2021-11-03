__all__ = ['read_cube']

import os
import sys
import pathlib
from datetime import datetime, timedelta
import numpy as np
import h5py
import astropy.wcs
import astropy.units as u
import sunpy.coordinates as coords
from eispac.core.eiscube import EISCube
from eispac.core.read_wininfo import read_wininfo
from eispac.instr.calc_read_noise import calc_read_noise

def read_cube(filename=None, window=0, apply_radcal=True, radcal=None,
              abs_errs=True, count_offset=None, debug=False):
    """Load a single window of EIS data from an HDF5 file into an EISCube object

    Parameters
    ----------
    filename : str or pathlib.Path object
        Name of either the data or head HDF5 file for a single EIS observation
    window : int, float, or str, optional
        Requested spectral window number or the value of any wavelength within
        the requested window. Default is '0'
    apply_radcal : bool, optional
        If set to True, will apply the pre-flight radiometric calibration curve
        found in the HDF5 header file and set units to erg/(cm^2 s sr). If set
        to False, will simply return the data in units of photon counts. Default
        is True.
    radcal : array_like, optional
        User-inputted radiometric calibration curve to be applied to the data.
    abs_errs : bool, optional
        If set to True, will calulate errors based on the absolute value of the
        counts. This allows for reasonable errors to be estimated for valid
        negative count values that are the result of the dark count subtraction
        method (not bad or filled data). Default it True.
    count_offset : int or float, optional
        Constant value to add to the count array before error estimate or
        calibration. Could be useful for testing data processing methods.
        Default is None (no count offset).
    debug : bool, optional
        If set to True, will return a dictionary with the raw counts and metadata
        instead of an EISCube class instance. Useful for examining data files that
        fail to load properly.

    Returns
    -------
    output_cube : EISCube class instance
        An EISCube class instance containing the requested window
    """
    ############################################################################
    ### Read data and header information from hdf5 files
    ############################################################################
    # Input type validation (value checks are implemented later)
    if not isinstance(filename, (str, pathlib.Path)):
        print('Error: Please input a valid filepath as '
             +'either a string or pathlib.Path object', file=sys.stderr)
        return None
    if not isinstance(window, (int, float, str)):
        print('Error: Please input a valid window or wavelength number '
             +'as either an integer, float, or string', file=sys.stderr)
        return None

    # Initialize "meta" dictionary to contain ALL of the extra EIS information.
    # We may want to add some of the values as attributes of the EISCube instead.
    meta = dict()

    # Parse filename and determine the directory and filename
    abs_filepath = pathlib.Path(filename).resolve()
    input_name = str(abs_filepath.name)
    input_dir = abs_filepath.parent
    if str(input_dir) == '.':
        input_dir = pathlib.Path().cwd()

    # Determine data and header filenames, regardless of which one was inputted
    data_filepath = input_dir.joinpath(input_name.replace('.head.h5', '.data.h5'))
    head_filepath = input_dir.joinpath(input_name.replace('.data.h5', '.head.h5'))

    # Check for data and head files. Exit if either file does not exist.
    if not data_filepath.is_file():
        print('Error: Data file does not exist, ' + str(data_filepath),
              file=sys.stderr)
        return None
    else:
        meta['filename_data'] = str(data_filepath)
        print('Data file,\n   ' + str(data_filepath))
    if not head_filepath.is_file():
        print('Error: Header file does not exist, ' + str(head_filepath),
              file=sys.stderr)
        return None
    else:
        meta['filename_head'] = str(head_filepath)
        print('Header file,\n   ' + str(head_filepath))

    # Read in min and max wavelength for each window in the file so we can
    # search for the requested window. Note: wininfo is a recarray with the
    # following fields: 'iwin', 'line_id', 'wvl_min', 'wvl_max', 'nl', & 'xs'
    wininfo = read_wininfo(head_filepath)
    num_win = wininfo.size
    meta['wininfo'] = wininfo

    # Locate the requested data window. Exit if it does not exist.
    if int(window) < 25:
        # Interpret values < 25 as window number
        if window >=0 and window < num_win:
            meta['iwin'] = window
            meta['iwin_str'] = f'win{window:02d}'
            print(f'Found window {window}')
        else:
            print(f'Error: Window {window} does not exist! The input data'
                 +f' file contains window numbers between 0 and {num_win}',
                  file=sys.stderr)
            return None
    else:
        # Interpret values > 25 as wavelength
        wvl = float(window)
        p = (wininfo['wvl_max'] - wvl)*(wvl - wininfo['wvl_min'])
        iwin = np.where(p >= 0)[0]
        if len(iwin) == 1:
            meta['iwin'] = iwin[0]
            meta['iwin_str'] = f'win{iwin[0]:02d}'
            print(f'Found a wavelength {wvl:.2f} [Angstroms] in window {iwin[0]}')
        else:
            print(f'Error: Wavelength not found! The input data file does'
                 +f' not contain a window that observes {window} Angstroms',
                  file=sys.stderr)
            return None

    # Read in the photon counts from data file
    with h5py.File(data_filepath, 'r') as f_data:
        lv_1_counts = np.array(f_data['level1/'+meta['iwin_str']])
        lv_1_count_units = f_data['level1/intensity_units'][0]
        lv_1_count_units = lv_1_count_units.decode('utf-8')

    # Read in metadata and instrumental correction factors from head file
    with h5py.File(head_filepath, 'r') as f_head:
        # Read index information (example data says it is from the level0 FITS file)
        index = {}
        for key in f_head['index']:
            val = np.array(f_head['index/'+key])
            if type(val[0]) == np.bytes_:
                val = val.astype(np.unicode_) # convert bytes to unicode
            if val.size == 1:
                val = val[0]
            index[key] = val

        meta['index'] = index

        # Read general EIS pointing information (includes average corrections)
        pointing = {}
        for key in f_head['pointing']:
            val = np.array(f_head['pointing/'+key])
            if type(val[0]) == np.bytes_:
                val = val.astype(np.unicode_) # convert bytes to unicode
            if val.size == 1:
                val = val[0]
            pointing[key] = val

        meta['pointing'] = pointing

        # Read calibration data
        meta['wave'] = np.array(f_head['wavelength/'+meta['iwin_str']])
        meta['radcal'] = np.array(f_head['radcal/'+meta['iwin_str']+'_pre'])
        meta['slit_width'] = np.array(f_head['instrumental_broadening/slit_width'])
        slit_width_units = f_head['instrumental_broadening/slit_width_units'][0]
        meta['slit_width_units'] = slit_width_units.decode('utf-8')
        meta['ccd_offset'] = np.array(f_head['ccd_offsets/'+meta['iwin_str']])

        # Read wavelength-dependent correction factor
        meta['wave_corr'] = np.array(f_head['wavelength/wave_corr'])
        meta['wave_corr_t'] = np.array(f_head['wavelength/wave_corr_t'])
        meta['wave_corr_tilt'] = np.array(f_head['wavelength/wave_corr_tilt'])

        # Read time and duration information
        try:
            meta['date_obs'] = np.array(f_head['times/date_obs']).astype(np.unicode_)
            meta['date_obs_format'] = np.array(f_head['times/time_format']).astype(np.unicode_)[0]
        except KeyError:
            print('WARNING: complete date_obs information for each raster step'
                 +' is missing in the HDF5 header file!')
        meta['duration'] = np.array(f_head['exposure_times/duration'])
        step_duration_units = f_head['exposure_times/duration_units'][0]
        meta['duration_units'] = step_duration_units.decode('utf-8')

    if debug == True:
        print('DEBUG MODE ON: returning dictionary with raw counts and metadata')
        return {'data':lv_1_counts, 'data_units':lv_1_count_units, 'meta':meta}

    # Check for user-inputted radcal curve
    if apply_radcal or radcal is not None:
        apply_radcal = True
        if radcal is not None:
            # Confirm dimensions are compatiable
            radcal_array = np.array(radcal)
            num_wave = len(meta['wave'])
            if len(radcal_array) != num_wave:
                print(f'Error: Input radcal array has the incorrect number of'
                     +f' elements. For the selected window, please input an'
                     +f' array with {num_wave} elements.',
                     file=sys.stderr)
                return None
        else:
            # Just use the pre-flight radcal curve
            radcal_array = meta['radcal']
    else:
        radcal_array = None

    ############################################################################
    ### Apply pointing corrections and create output EISCube
    ############################################################################

    ### (1) Apply the AIA offset corrections
    x_center = pointing['xcen'] + pointing['offset_x']
    y_center = pointing['ycen'] + pointing['offset_y']

    ### (2) Compute mean ccd offset for the current window and apply to y_center
    # Note_1: 'ccd_offsets' are in units of [pixels] while 'y_center' is in [arcsec].
    #         However, this should not normally be problem since the EIS y-scale
    #         is 1 [arcsec]/[pixel]
    # Note_2: technically, the ccd_offset varies with wavelength. However, within
    #         a single window, the difference is on the order of 0.05 arcsec
    mean_ccd_offset = np.mean(meta['ccd_offset'])
    y_center = y_center - mean_ccd_offset

    ### (3) Apply wave correction and get median base wavelength value and delta
    counts_shape = lv_1_counts.shape
    if len(counts_shape) > 3:
        print('ERROR: EIS observations with multiple rasters in a single file'
             +' are not currently supported.', file=sys.stderr)
        return None
    ny_pxls = counts_shape[0] # num pixels along the slit (y-axis)
    nx_steps = counts_shape[1] # num raster steps (x-axis)
    n_wave = counts_shape[2] # num wavelength values (i.e. EIS window width)
    corrected_wave = np.zeros(counts_shape)
    for i in range(ny_pxls):
        for j in range(nx_steps):
            corrected_wave[i,j,:] = meta['wave'] - meta['wave_corr'][i,j]

    base_wave = np.median(corrected_wave[:,:,0])
    wave_delt = np.median(np.diff(corrected_wave, axis=2))

    ### (4) Calculate reference pixel coords and fov
    x1 = x_center - pointing['x_scale']*nx_steps/2.0
    x2 = x_center + pointing['x_scale']*nx_steps/2.0
    y1 = y_center - pointing['y_scale']*ny_pxls/2.0
    y2 = y_center + pointing['y_scale']*ny_pxls/2.0

    ### (5) Extract timestamps and calculate the mid point of the observation
    date_obs = index['date_obs']
    date_end = index['date_end']
    date_diff = datetime.fromisoformat(date_end) - datetime.fromisoformat(date_obs)
    date_avg = datetime.fromisoformat(date_obs) + date_diff/2.0
    date_avg = date_avg.isoformat(timespec='milliseconds') # convert to string

    ### (6) Fetch the observer location in heliographic coords
    hg_coords = coords.get_body_heliographic_stonyhurst('earth', time=date_obs)

    ### (7) Create a new header dict updated values
    # Note: the order of axes here is the same as the original fits index
    output_hdr = dict()

    output_hdr['naxis'] = 3
    output_hdr['date_obs'] = date_obs
    output_hdr['date_beg'] = date_obs
    output_hdr['date_avg'] = date_avg
    output_hdr['date_end'] = date_end
    output_hdr['timesys'] = 'UTC'

    output_hdr['telescop'] = 'Hinode'
    output_hdr['instrume'] = 'EIS'
    output_hdr['stud_acr'] = index['stud_acr']
    output_hdr['obstitle'] = index['obstitle']
    output_hdr['target'] = index['target']
    output_hdr['sci_obj'] = index['sci_obj']
    output_hdr['line_id'] = wininfo['line_id'][meta['iwin']]
    output_hdr['measrmnt'] = 'intensity'
    output_hdr['bunit'] = 'unknown' # units of primary observable
    output_hdr['slit_id'] = index['slit_id']
    output_hdr['slit_ind'] = index['slit_ind']
    output_hdr['nraster'] = index['nraster']

    output_hdr['naxis1'] = nx_steps
    output_hdr['cname1'] = 'Solar-X'
    output_hdr['crval1'] = x1
    output_hdr['crpix1'] = 1
    output_hdr['cdelt1'] = pointing['x_scale']
    output_hdr['ctype1'] = 'HPLN-TAN'
    output_hdr['cunit1'] = 'arcsec'

    output_hdr['naxis2'] = ny_pxls
    output_hdr['cname2'] = 'Solar-Y'
    output_hdr['crval2'] = y1
    output_hdr['crpix2'] = 1
    output_hdr['cdelt2'] = pointing['y_scale']
    output_hdr['ctype2'] = 'HPLT-TAN'
    output_hdr['cunit2'] = 'arcsec'

    output_hdr['naxis3'] = n_wave
    output_hdr['cname3'] = 'Wavelength'
    output_hdr['crval3'] = base_wave
    output_hdr['crpix3'] = 1
    output_hdr['cdelt3'] = wave_delt
    output_hdr['ctype3'] = 'WAVE'
    output_hdr['cunit3'] = 'Angstrom'

    output_hdr['fovx'] = x2 - x1
    output_hdr['fovy'] = y2 - y1
    output_hdr['xcen'] = x1 + 0.5*(x2-x1)
    output_hdr['ycen'] = y1 + 0.5*(y2-y1)

    output_hdr['hgln_obs'] = hg_coords.lon.deg
    output_hdr['hglt_obs'] = hg_coords.lat.deg
    output_hdr['dsun_obs'] = hg_coords.radius.m

    # Calculate and append extra keys to meta dict for user convenience
    meta['mod_index'] = output_hdr
    meta['aspect'] = pointing['y_scale']/pointing['x_scale']
    meta['aspect_ratio'] = pointing['y_scale']/pointing['x_scale'] # DEPRICATED
    meta['extent_arcsec'] = [x1, x2, y1, y2] # [left, right, bottom, top]
    meta['notes'] = []

    try:
        # Create the WCS object
        # NB: For some reason, the order of axes in the WCS is reversed relative
        #     to the data array inside an NDCube. We should take care to fully
        #     document this for our users.
        clean_wcs = astropy.wcs.WCS(output_hdr, fix=True)
        clean_wcs = clean_wcs.swapaxes(0,1) # swap x with y
        clean_wcs = clean_wcs.swapaxes(0,2) # now swap y with wavelength

        # Add a user-supplied constant value to the count array
        if count_offset is not None:
            lv_1_counts = lv_1_counts + count_offset

        # Calculate errors due to Poisson noise and read noise
        # Also create data mask to flag invalid values
        read_noise = calc_read_noise(corrected_wave)
        if abs_errs == True:
            data_mask = lv_1_counts <= -100 # mask missing data
            lv_1_count_errs = np.sqrt(np.abs(lv_1_counts) + read_noise**2)
        else:
            data_mask = lv_1_counts <= 0 # mask ALL negative or zero values
            clean_lv_1_counts = lv_1_counts.copy()
            clean_lv_1_counts[data_mask] = 0.0
            lv_1_count_errs = np.sqrt(clean_lv_1_counts + read_noise**2)

        lv_1_count_errs[data_mask] = -100 # EIS missing data err value (in IDL)
        if apply_radcal:
            cube_data = lv_1_counts*radcal_array
            cube_errs = lv_1_count_errs*radcal_array
            data_units = 'erg / (cm2 s sr)'
        else:
            cube_data = lv_1_counts
            cube_errs = lv_1_count_errs
            data_units = 'photon'

        meta['mod_index']['bunit'] = data_units
        output_cube = EISCube(cube_data, wcs=clean_wcs, uncertainty=cube_errs,
                              wavelength=corrected_wave, radcal=radcal_array,
                              meta=meta, unit=data_units, mask=data_mask)
    except:
        print('Error: Failed to initialize WCS or EISCube instance due to bad'
             +' header data. Please report this issue to the eispac development'
             + 'team', file=sys.stderr)
        return None

    return output_cube

if __name__ == '__main__':

    filename = '../data/eis_20190404_131513.data.h5'
    wave = 195.119
    eis = read_cube(filename, wave)

    print(f'date_obs = ' + eis.meta['index']['date_obs'])
    shape = eis.data.shape
    print(f'data shape = {shape}')
