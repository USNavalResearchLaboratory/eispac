__all__ = ['read_cube']

import os
import sys
import pathlib
import numpy as np
import h5py
import astropy.wcs
import astropy.units as u
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
    # We may want to add some of the values as attributes in the EISCube instead.
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

    # Read in min and max wavelength for each window in the file
    #   so we can search for the requested window
    wininfo = read_wininfo(head_filepath)
    # with h5py.File(head_filepath, 'r') as f_head:
    #     num_win = f_head['/wininfo/nwin'][0]
    #     dt = np.dtype([('iwin', 'i4'), ('line_id', 'U64'),
    #                    ('wvl_min', 'f'), ('wvl_max','f'),
    #                    ('nl', 'i4'), ('xs','i4')])
    #     wininfo = np.recarray((num_win,), dtype=dt)
    #     # wininfo = {'iwin': np.zeros(num_win, dtype='i4'),
    #     #            'line_id': np.zeros(num_win, dtype='U64'),
    #     #            'wvl_min': np.zeros(num_win, dtype='f'),
    #     #            'wvl_max': np.zeros(num_win, dtype='f'),
    #     #            'nl': np.zeros(num_win, dtype='i4'),
    #     #            'xs':np.zeros(num_win, dtype='i4')}
    #     for iwin in range(num_win):
    #         line_id = f_head[f'/wininfo/win{iwin:02d}/line_id'][0]
    #         wvl_min = f_head[f'/wininfo/win{iwin:02d}/wvl_min'][0]
    #         wvl_max = f_head[f'/wininfo/win{iwin:02d}/wvl_max'][0]
    #         win_nl = f_head[f'/wininfo/win{iwin:02d}/nl'][0]
    #         win_xs = f_head[f'/wininfo/win{iwin:02d}/xs'][0]
    #         wininfo[iwin].iwin = iwin
    #         wininfo[iwin].line_id = line_id.decode('utf-8')
    #         wininfo[iwin].wvl_min = wvl_min
    #         wininfo[iwin].wvl_max = wvl_max
    #         wininfo[iwin].nl = win_nl # NEW
    #         wininfo[iwin].xs = win_xs # NEW
    #         # wininfo['iwin'][iwin] = iwin
    #         # wininfo['line_id'][iwin] = line_id.decode('utf-8')
    #         # wininfo['wvl_min'][iwin] = wvl_min
    #         # wininfo['wvl_max'][iwin] = wvl_max
    #         # wininfo['nl'][iwin] = win_nl
    #         # wininfo['xs'][iwin] = win_xs

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
    x_center = pointing['xcen']+pointing['offset_x']
    y_center = pointing['ycen']+pointing['offset_y']

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

    ### (4) Check values of cdelt and reference pixel coords?
    # Should we also check crota1 and crota2 here too?

    ### (5) Create a clean header dict with only the values needed by astropy.wcs.WCS
    # Note: for some reason, the order of axes in the WCS is reversed relative
    #       to the input data array. We should be careful to document this for users
    basic_hdr = dict()

    basic_hdr['crpix1'] = 0
    basic_hdr['crval1'] = base_wave
    basic_hdr['cdelt1'] = wave_delt
    basic_hdr['naxis1'] = n_wave
    basic_hdr['ctype1'] = 'Wavelength'
    basic_hdr['cunit1'] = 'Angstrom'

    basic_hdr['crpix2'] = (nx_steps+1)/2.0
    basic_hdr['crval2'] = x_center
    basic_hdr['cdelt2'] = pointing['x_scale']
    basic_hdr['naxis2'] = nx_steps
    basic_hdr['ctype2'] = 'HPLN-TAN' # 'Solar-X'
    basic_hdr['cunit2'] = 'arcsec'

    basic_hdr['crpix3'] = (ny_pxls+1)/2.0
    basic_hdr['crval3'] = y_center
    basic_hdr['cdelt3'] = pointing['y_scale']
    basic_hdr['naxis3'] = ny_pxls
    basic_hdr['ctype3'] = 'HPLT-TAN' # 'Solar-Y'
    basic_hdr['cunit3'] = 'arcsec'

    # Calculate and append extra keys to meta dict for user convenience
    meta['aspect_ratio'] = pointing['y_scale']/pointing['x_scale']
    meta['extent_arcsec'] = [0.0, 0.0, 0.0, 0.0] # [left, right, bottom, top]
    if nx_steps > 1: #proper raster
        meta['extent_arcsec'][0] = pointing['solar_x'][0]+pointing['offset_x']
        meta['extent_arcsec'][1] = pointing['solar_x'][-1]+pointing['offset_x']
    else: # single sit-and-stare slit observation
        meta['extent_arcsec'][0] = pointing['solar_x']+pointing['offset_x']
        meta['extent_arcsec'][1] = pointing['solar_x']+pointing['offset_x']+1
    meta['extent_arcsec'][2] = pointing['solar_y'][0]+pointing['offset_y']-mean_ccd_offset
    meta['extent_arcsec'][3] = pointing['solar_y'][-1]+pointing['offset_y']-mean_ccd_offset
    meta['notes'] = []

    try:
        # Create the WCS object
        clean_wcs = astropy.wcs.WCS(basic_hdr, fix=True)

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
