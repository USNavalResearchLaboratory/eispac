__all__ = ['rebin_cube']

import os
import sys
import pathlib
import numpy as np
import h5py
import astropy.wcs
import astropy.units as u
from ndcube import NDCube
from eispac.core.read_cube import EISCube

# try:
#     # Import local versions of submodules
#     from read_cube import EISCube
#     print('Notice: Loading local version of eispac submodules')
# except:
#     # Import from installed package
#     from eispac.core.read_cube import EISCube

class EISCube(NDCube):
    # NOTE: this is based on the example given at
    # https://docs.astropy.org/en/stable/nddata/subclassing.html#slicing-an-additional-property
    def __init__(self, *args, **kwargs):
        # Remove wavelength attribute, if given, and pass it to the setter.
        self.wavelength = kwargs.pop('wavelength') if 'wavelength' in kwargs else None
        super().__init__(*args, **kwargs)

    @property
    def wavelength(self):
        return self._wavelength

    @wavelength.setter
    def wavelength(self, value):
        self._wavelength = value

    def _slice(self, item):
        # slice all normal attributes
        kwargs = super()._slice(item)
        # The arguments for creating a new instance are saved in kwargs. So we
        # need to add another keyword "wavelength" and add the sliced wavelengths
        kwargs['wavelength'] = self.wavelength[item]
        return kwargs # these must be returned

    def total_intensity(self):
        """
        Sum the intensity over the entire spectral window.
        """
        data = np.sum(self.data, axis=2)
        wcs = self.wcs.dropaxis(0)
        return NDCube(data, wcs)

def rebin_cube(datacube):
    """Rebin an EISCube in order to improve the signal-to-noise ratio

    Parameters
    ----------
    datacube : EISCube or simple NDCube instance
        Datacube containing EIS observations and a WCS array

    Returns
    -------
    output_cube : EISCube class instance
        An EISCube class instance containing the rebinned data
    """
    ############################################################################
    ### Extract meta dictionary and other information
    ############################################################################
    input_meta = datacube.meta
    input_wcs = datacube.wcs

    # ############################################################################
    # ### Read data and header information from hdf5 files
    # ############################################################################
    # # Input type validation (value checks are implemented later)
    # if not isinstance(filename, (str, pathlib.Path)):
    #     sys.exit('Error: Please input a valid filepath as '
    #             +'either a string or pathlib.Path object')
    # if not isinstance(window, (int, float, str)):
    #     sys.exit('Error: Please input a valid window or wavelength number '
    #             +'as either an integer, float, or string')
    #
    # # Initialize "meta" dictionary to contain ALL of the extra EIS information.
    # # We may want to add some of the values as attributes in the EISCube instead.
    # meta = dict()
    #
    # # Parse filename and determine the directory and filename
    # abs_filepath = pathlib.Path(filename).resolve()
    # input_name = str(abs_filepath.name)
    # input_dir = abs_filepath.parent
    # if str(input_dir) == '.':
    #     input_dir = pathlib.Path().cwd()
    #
    # # Determine data and header filenames, regardless of which one was inputted
    # data_filepath = input_dir.joinpath(input_name.replace('.head.h5', '.data.h5'))
    # head_filepath = input_dir.joinpath(input_name.replace('.data.h5', '.head.h5'))
    #
    # # Check for data and head files. Exit if either file does not exist.
    # if not data_filepath.is_file():
    #     sys.exit('Error: Data file does not exist, ' + str(data_filepath))
    # else:
    #     meta['filename_data'] = str(data_filepath)
    #     print('Data file,\n   ' + str(data_filepath))
    # if not head_filepath.is_file():
    #     sys.exit('Error: Header file does not exist, ' + str(head_filepath))
    # else:
    #     meta['filename_head'] = str(head_filepath)
    #     print('Header file,\n   ' + str(head_filepath))
    #
    # # Read in min and max wavelength for each window in the file
    # #   so we can search for the requested window
    # with h5py.File(head_filepath, 'r') as f_head:
    #     num_win = f_head['/wininfo/nwin'][0]
    #     # dt = np.dtype([('iwin', 'i4'), ('line_id', 'U64'),
    #     #                ('wvl_min', 'f'), ('wvl_max','f'),
    #     #                ('nl', 'i4'), ('xs','i4')])
    #     # wininfo = np.recarray((num_win,), dtype=dt)
    #     wininfo = {'iwin': np.zeros(num_win, dtype='i4'),
    #                'line_id': np.zeros(num_win, dtype='U32'),
    #                'wvl_min': np.zeros(num_win, dtype='f'),
    #                'wvl_max': np.zeros(num_win, dtype='f'),
    #                'nl': np.zeros(num_win, dtype='i4'),
    #                'xs':np.zeros(num_win, dtype='i4')}
    #
    #     for iwin in range(num_win):
    #         line_id = f_head[f'/wininfo/win{iwin:02d}/line_id'][0]
    #         wvl_min = f_head[f'/wininfo/win{iwin:02d}/wvl_min'][0]
    #         wvl_max = f_head[f'/wininfo/win{iwin:02d}/wvl_max'][0]
    #         win_nl = f_head[f'/wininfo/win{iwin:02d}/nl'][0]
    #         win_xs = f_head[f'/wininfo/win{iwin:02d}/xs'][0]
    #         wininfo['iwin'][iwin] = iwin
    #         wininfo['line_id'][iwin] = line_id.decode('utf-8')
    #         wininfo['wvl_min'][iwin] = wvl_min
    #         wininfo['wvl_max'][iwin] = wvl_max
    #         wininfo['nl'][iwin] = win_nl
    #         wininfo['xs'][iwin] = win_xs
    #
    # meta['wininfo'] = wininfo
    #
    # # Locate the requested data window. Exit if it does not exist.
    # if int(window) < 25:
    #     # Interpret values < 25 as window number
    #     if window >=0 and window < num_win:
    #         meta['iwin'] = window
    #         meta['iwin_str'] = f'win{window:02d}'
    #         print(f'Found window {window}')
    #     else:
    #         sys.exit(f'Error: Window {window} does not exist! The input data '
    #                 +f'file contains window numbers between 0 and {num_win}')
    # else:
    #     # Interpret values > 25 as wavelength
    #     wvl = float(window)
    #     p = (wininfo['wvl_max'] - wvl)*(wvl - wininfo['wvl_min'])
    #     iwin = np.where(p >= 0)[0]
    #     if len(iwin) == 1:
    #         meta['iwin'] = iwin[0]
    #         meta['iwin_str'] = f'win{iwin[0]:02d}'
    #         print(f'Found a wavelength {wvl:.2f} [Angstroms] in window {iwin[0]}')
    #     else:
    #         sys.exit(f'Error: Wavelength not found! The input data file does '
    #                 +f'not contain a window that observes {window} Angstroms')
    #
    # # Read in the photon counts from data file
    # with h5py.File(data_filepath, 'r') as f_data:
    #     lv_1_counts = np.array(f_data['level1/'+meta['iwin_str']])
    #     # lv_1_count_units = np.array(f_data['level1/intensity_units'])
    #     lv_1_count_units = f_data['level1/intensity_units'][0]
    #     lv_1_count_units = lv_1_count_units.decode('utf-8')
    #
    #
    # # Read in metadata and instrumental correction factors from head file
    # with h5py.File(head_filepath, 'r') as f_head:
    #     # Read index information (example data says it is from the level0 FITS file)
    #     index = {}
    #     for key in f_head['index']:
    #         val = np.array(f_head['index/'+key])
    #         if type(val[0]) == np.bytes_:
    #             val = val.astype(np.unicode_) # convert bytes to unicode
    #         if val.size == 1:
    #             val = val[0]
    #         index[key] = val
    #
    #     meta['index'] = index
    #
    #     # Read general EIS pointing information (includes average corrections)
    #     pointing = {}
    #     for key in f_head['pointing']:
    #         val = np.array(f_head['pointing/'+key])
    #         if type(val[0]) == np.bytes_:
    #             val = val.astype(np.unicode_) # convert bytes to unicode
    #         if val.size == 1:
    #             val = val[0]
    #         pointing[key] = val
    #
    #     meta['pointing'] = pointing
    #
    #     # Read calibration data
    #     meta['wave'] = np.array(f_head['wavelength/'+meta['iwin_str']])
    #     meta['radcal'] = np.array(f_head['radcal/'+meta['iwin_str']+'_pre'])
    #     meta['slit_width'] = np.array(f_head['instrumental_broadening/slit_width'])
    #     slit_width_units = f_head['instrumental_broadening/slit_width_units'][0]
    #     meta['slit_width_units'] = slit_width_units.decode('utf-8') # NEW
    #     meta['ccd_offset'] = np.array(f_head['ccd_offsets/'+meta['iwin_str']])
    #
    #     # Read wavelength-depandent correction factor
    #     meta['wave_corr'] = np.array(f_head['wavelength/wave_corr'])
    #     meta['wave_corr_t'] = np.array(f_head['wavelength/wave_corr_t']) # NEW
    #     meta['wave_corr_tilt'] = np.array(f_head['wavelength/wave_corr_tilt']) # NEW
    #
    #     # Read time and duration information
    #     try:
    #         meta['timestamp'] = np.array(f_head['times/date_obs']).astype(np.unicode_)
    #         meta['timestamp_format'] = np.array(f_head['times/time_format']).astype(np.unicode_)[0]
    #     except KeyError:
    #         print('WARNING: timestamps are missing in the HDF5 header file!')
    #     meta['duration'] = np.array(f_head['exposure_times/duration']) # NEW
    #     step_duration_units = f_head['exposure_times/duration_units'][0]
    #     meta['duration_units'] = step_duration_units.decode('utf-8') # NEW
    #
    # ############################################################################
    # ### Apply pointing corrections and create output EISCube
    # ############################################################################
    #
    # ### (1) Apply the AIA offset corrections
    # x_center = pointing['xcen']+pointing['offset_x']
    # y_center = pointing['ycen']+pointing['offset_y']
    #
    # ### (2) Compute mean ccd offset for the current window and apply to y_center
    # # Note_1: 'ccd_offsets' are in units of [pixels] while 'y_center' is in [arcsec].
    # #         However, this should not normally be problem since the EIS y-scale
    # #         is 1 [arcsec]/[pixel]
    # # Note_2: technically, the ccd_offset varies with wavelength. However, within
    # #         a single window, the difference is on the order of 0.05 arcsec
    # mean_ccd_offset = np.mean(meta['ccd_offset'])
    # y_center = y_center - mean_ccd_offset
    #
    # ### (3) Apply wave correction and get median base wavelength value and delta
    # counts_shape = lv_1_counts.shape
    # ny_pxls = counts_shape[0] # num y pixels along the slit
    # nx_steps = counts_shape[1] # num raster steps in the x direction
    # n_wave = counts_shape[2] # num wavelength values (i.e. EIS window width)
    # corrected_wave = np.zeros(counts_shape)
    # for i in range(ny_pxls):
    #     for j in range(nx_steps):
    #         corrected_wave[i,j,:] = meta['wave'] - meta['wave_corr'][i,j]
    #
    # base_wave = np.median(corrected_wave[:,:,0])
    # wave_delt = np.median(np.diff(corrected_wave, axis=2))
    #
    # ### (4) Check values of cdelt and reference pixel coords?
    # # Should we also check crota1 and crota2 here too?
    #
    # ### (5) Create a clean header dict with only the values needed by astropy.wcs.WCS
    # # Note: for some reason, the order of axes in the WCS is reversed relative
    # #       to the input data array. We should be careful to document this for users
    # basic_hdr = dict()
    #
    # basic_hdr['crpix1'] = 1
    # basic_hdr['crval1'] = base_wave
    # basic_hdr['cdelt1'] = wave_delt
    # basic_hdr['naxis1'] = n_wave
    # basic_hdr['ctype1'] = 'Wavelength'
    # basic_hdr['cunit1'] = 'Angstrom'
    #
    # basic_hdr['crpix2'] = (nx_steps+1)/2.0
    # basic_hdr['crval2'] = x_center
    # basic_hdr['cdelt2'] = pointing['x_scale']
    # basic_hdr['naxis2'] = nx_steps
    # # basic_hdr['ctype2'] = 'Solar-X'
    # basic_hdr['ctype2'] = 'HPLN-TAN'
    # basic_hdr['cunit2'] = 'arcsec'
    #
    # basic_hdr['crpix3'] = (ny_pxls+1)/2.0
    # basic_hdr['crval3'] = y_center
    # basic_hdr['cdelt3'] = pointing['y_scale']
    # basic_hdr['naxis3'] = ny_pxls
    # # basic_hdr['ctype3'] = 'Solar-Y'
    # basic_hdr['ctype3'] = 'HPLT-TAN'
    # basic_hdr['cunit3'] = 'arcsec'
    #
    # # Create the WCS object
    # clean_wcs = astropy.wcs.WCS(basic_hdr, fix=True)
    #
    # # Calculate intensity values and then create the EISCube
    # lv_1_inten = lv_1_counts*meta['radcal']
    # lv_1_inten = u.Quantity(lv_1_inten, 'erg / (cm2 s sr)')
    #
    # loc_neg_counts = np.where(lv_1_counts < 0)
    # lv_1_counts[loc_neg_counts] = 0.0
    # lv_1_inten_errs = np.sqrt(lv_1_counts)*meta['radcal']
    #
    # inten_mask = lv_1_inten < 0*lv_1_inten.unit
    #
    # output_cube = EISCube(lv_1_inten, clean_wcs, wavelength=corrected_wave,
    #                       uncertainty=lv_1_inten_errs, mask=inten_mask,
    #                       meta=meta)

    return output_cube

if __name__ == '__main__':

    # filename = './data/eis_20190404_131513.data.h5'
    # wave = 195.119
    # eis = read_cube(filename, wave)
    #
    # print(f'date_obs = ' + eis.meta['index']['date_obs'])
    # shape = eis.data.shape
    # print(f'data shape = {shape}')
