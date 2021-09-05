__all__ = ['EISFitResult', 'create_fit_dict']

import copy
from datetime import datetime
import numpy as np
from scipy.ndimage import shift as shift_img
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from eispac import __version__ as eispac_version
import eispac.core.fitting_functions as fit_fns
from eispac.core.read_template import create_funcinfo
from eispac.util.rot_xy import rot_xy
from eispac.instr.ccd_offset import ccd_offset

class EISFitResult:
    """Object containing the results from fitting one or more EIS window spectra

    Parameters
    ----------
    wave : array_like
        Wavelength values of the data being fit. Used mostly to
        determine the correct dimensions for the output arrays.
    template : dict
        Fit template parameters and metadata contained in the 'template' attribute
        of an EISFitTemplate object.
    parinfo : dict
        Fit parameter initial values and constraints in the format expectd by
        mpfit. Normally found in the 'parinfo' attribute of an EISFitTemplate object
    func_name : str, optional
        String name of function that will be fit to the data. Must be one of the
        functions defined in the eispac.fitting_functions submodule. Default is
        'multigaussian'.

    Attributes
    ----------
    template : dict
        Full copy of the input template dictionary
    parinfo : list of dicts
        Full copy of the input parinfo list
    funcinfo : dict
        Function component information generated from the template dict
    n_pxls : int
        Number of pixels along the y-axis
    n_steps : int
        Number of raster steps or sit-and-stare exposures in the observation
    n_wave : int
        Number of wavelength values in the window
    n_gauss : int
        Number of Guassian functions in the fit
    n_poly : int
        Number of terms in the polynomial background
    func_name : str
        Name of the function fit to the data
    fit_func : function
        Copy of the actual Python function named by 'fit_func'
    fit : dict
        Dictionary of output fit parameters
    """

    def __init__(self, wave=None, template=None, parinfo=None,
                 func_name='multigaussian',
                 data_units='unknown', radcal='unknown', empty=False):
        self.date_fit = datetime.now().replace(microsecond=0).isoformat()
        self.eispac_version = eispac_version
        self.meta = dict()
        self.template = copy.deepcopy(template)
        self.parinfo = copy.deepcopy(parinfo)
        self.funcinfo = None
        self.n_pxls = 11
        self.n_steps = 22
        self.n_wave = 33
        self.n_gauss = 1
        self.n_poly = 1
        self.fit_module = 'unknown'
        self.fit_method = 'unknown'
        self.func_name = func_name
        self.fit_func = None
        self.fit = None
        self.data_units = data_units
        self._current_radcal = radcal

        # Determine input data shape, create fit dictionary, and copy misc. info
        if not empty:
            num_dims = wave.ndim
            dims_size = wave.shape
            self.n_wave = dims_size[-1] # last dim is always wavelength
            if num_dims == 1:
                self.n_pxls = 1
                self.n_steps = 1
            elif num_dims == 2:
                self.n_pxls = dims_size[0]
                self.n_steps = 1
            elif num_dims == 3:
                self.n_pxls = dims_size[0]
                self.n_steps = dims_size[1]

            self.n_gauss = template['n_gauss']
            self.n_poly = template['n_poly']
            self.funcinfo = create_funcinfo(self.template)
            self.fit_func = getattr(fit_fns, func_name)

            fit = create_fit_dict(self.n_pxls, self.n_steps, self.n_wave,
                                  self.n_gauss, self.n_poly,
                                  data_units=self.data_units)
            line_ids = template['line_ids']
            # if isinstance(line_ids, np.ndarray):
            #     line_ids = line_ids.tolist() # multiple entries = array of strings
            # else:
            #     line_ids = [line_ids.decode('utf-8')] # single entry = byte array
            fit['line_ids'] = line_ids
            fit['main_component'] = template['component']
            fit['wave_range'][0] = template['data_x'][0]
            fit['wave_range'][1] = template['data_x'][-1]
            self.fit = fit
        else:
            # Initialize an empty EISFitResult (used by eispac.read_fit)
            self.fit = create_fit_dict(11, 22, 33, 1, 1)

    @property
    def radcal(self):
        return self._current_radcal

    @radcal.setter
    def radcal(self, input_array):
        print('Error: Please use the .apply_radcal() and .remove_radcal()'
             +' methods to modify or change the radiometric calibration.')

    def _validate_component_num(self, component):
        """Helper function for validating functional component number.

        Note: unless there is an error, the output is a numpy array
        """
        num_comp = self.n_gauss + 1 if self.n_poly > 0 else self.n_gauss
        if component is not None:
            if np.size(component) == 1:
                use_comp = np.array([component], dtype=int)
            else:
                use_comp = np.array(component, dtype=int)
            use_comp = use_comp.round(0).astype(int)
            for c in range(use_comp.size):
                if use_comp[c] < 0:
                    # Converted a negative index to the correct positive value
                    use_comp[c] = num_comp + use_comp[c]
                elif use_comp[c] >= num_comp:
                    print('Error: input component number is too large!')
                    use_comp = 'error'
        else:
            use_comp = None

        return use_comp

    def _validate_coords(self, coords):
        """Helper function for validating pixel coordinate pairs.
        """
        if coords is not None:
            if np.size(coords) == 2:
                use_coords = [int(coords[0]), int(coords[1])]
                if use_coords[0] < 0:
                    use_coords[0] = self.n_pxls + use_coords[0]
                if use_coords[1] < 0:
                    use_coords[1] = self.n_steps + use_coords[1]
                if use_coords[0] >= self.n_pxls or use_coords[1] >= self.n_steps:
                    print('Error: requested coordinates are outside the range'
                          +' of available results!')
                    use_coords = 'error'
            else:
                print('Error: please input a valid coordinate pair or'
                      +' "set coords=None"')
                use_coords = 'error'
        else:
            use_coords = None

        return use_coords

    def _validate_param_name(self, param_name, casefold):
        """Helper function for validating a single parameter name.
        """
        if param_name is not None and not isinstance(param_name, str):
            print('Error: please set param_name to either a single string value'
                  +' or None')
            use_name = "0"
        else:
            use_name = param_name

        if casefold:
            use_name = use_name.lower()

        return use_name

    def _get_param_filter(self, component=None, param_name=None):
        print('Sorry, _get_param_filter() is not implemented at this time.')
        return None

    def get_params(self, component=None, param_name=None, coords=None,
                   casefold=False):
        """Extract parameters values by component number, name, or pixel coords

        Parameters
        ----------
        component : int or list, optional
            Integer number (or list of ints) of the functional component(s).
            If set to None, will return the total combined fit profile.
            Default is None.
        param_name : str, optional
            String name of the requested parameter. If set to None, will not
            filter based on paramater name. Default is None
        coords : list or tupple, optional
            (Y, X) coordinates of the requested datapoint. If set to None, will
            instead return the parameters at all locations. Default is None
        casefold : bool, optional
            If set to True, will ignore case when extracting parameters by
            name. Default is False.

        Returns
        -------
        param_vals : numpy array
            Parameter values
        param_errs : numpy array
            Estimated parameter errors
        """
        # Validate input values
        use_comp = self._validate_component_num(component)
        use_coords = self._validate_coords(coords)
        use_name = self._validate_param_name(param_name, casefold)
        if (isinstance(use_coords, str)) or (isinstance(use_comp, str)) or (use_name == '0'):
            return None, None

        # Now, create create index array and apply filters
        num_params = self.fit['param_names'].size
        p_name_arr = self.fit['param_names']
        if casefold:
            p_name_arr = p_name_arr.lower()
        p_ind = np.full(num_params, False, dtype=bool)
        if (param_name is None) and (component is None):
            p_ind[:] = True

        if param_name is not None:
            loc_name = np.char.startswith(p_name_arr, use_name)
        else:
            loc_name = np.full(num_params, True, dtype=bool)

        if component is not None:
            for c in range(use_comp.size):
                loc_comp = np.where((self.fit['component'] == use_comp[c])
                                    & (loc_name == True))
                p_ind[loc_comp] = True
        elif param_name is not None:
            p_ind[loc_name] = True

        # Finally, extract the parameters at the given coords
        if coords is not None:
            param_vals = self.fit['params'][coords[0], coords[1], p_ind]
            param_errs = self.fit['perror'][coords[0], coords[1], p_ind]
        else:
            param_vals = self.fit['params'][:, :, p_ind]
            param_errs = self.fit['perror'][:, :, p_ind]

        return param_vals, param_errs

    def get_fit_profile(self, component=None, coords=None, num_wavelengths=None,
                        wave_range='auto', use_mask=True):
        """Calculate the fit intensity profile (total or component) at a location.

        Parameters
        ----------
        component : int or list, optional
            Integer number (or list of ints) of the functional component(s).
            If set to None, will return the total combined fit profile.
            Default is None.
        coords : list or tupple, optional
            (Y, X) coordinates of the requested datapoint. If set to None, will
            instead return the parameters at all locations. Default is None
        num_wavelengths : int, optional
            Number of wavelength values to compute the fit intensity at. These
            values will be equally spaced and span the entire fit window. If set
            to None, will use the observed wavelength values. Default is None.

        Returns
        -------
        fit_wave : numpy array
            Wavelength values
        fit_inten : numpy array
            Fit intensity values
        """
        # Validate input values
        use_comp = self._validate_component_num(component)
        use_coords = self._validate_coords(coords)
        if (isinstance(use_coords, str)) or (isinstance(use_comp, str)):
            return None, None

        # Validate wave range for interpolated results
        # Note: a range of [0, 0] will be replaced with the obs. window range
        if str(wave_range).lower() == 'auto':
            use_range = self.fit['wave_range']
        elif str(wave_range).lower() == 'none':
            use_range = [0, 0]
        elif len(wave_range) == 2:
            use_range = wave_range
        else:
            print('Warning: incorrect number of wave_range values.'
                 +' Returning full spectral window instead.')
            use_range = [0, 0]

        if num_wavelengths is not None and use_range[0] == 0:
            if coords is not None:
                # single spectrum
                # Note: the paren are needed here to unpack use_coords
                use_range[0] = self.fit['wavelength'][(*use_coords, 0)]
                use_range[1] = self.fit['wavelength'][(*use_coords, -1)]
            else:
                # full image
                use_range[0] = np.median(self.fit['wavelength'][:,:,0])
                use_range[1] = np.median(self.fit['wavelength'][:,:,-1])

        # Determine numbers and types of components in output profile
        full_num_comp = self.n_gauss + 1 if self.n_poly > 0 else self.n_gauss
        if use_comp is None:
            num_gauss = self.n_gauss
            num_poly = self.n_poly
        else:
            if (self.n_poly > 0) and (full_num_comp-1 in use_comp):
                num_gauss = use_comp.size - 1
                num_poly = self.n_poly
            else:
                num_gauss = use_comp.size
                num_poly = 0

        # Extract the parameters for either the total or component profile
        # Note: all poly terms are considered part of a single background component
        param_vals, param_errs = self.get_params(component=use_comp, coords=use_coords)

        # Create output wavelength array and then calculate the fit profile
        if coords is not None:
            # Single spectrum
            if num_wavelengths is None:
                fit_wave = self.fit['wavelength'][use_coords[0], use_coords[1], :]
            else:
                fit_wave = np.linspace(use_range[0], use_range[-1], num_wavelengths)
            if self.fit['status'][use_coords[0], use_coords[1]] > 0:
                fit_inten = self.fit_func(param_vals, fit_wave, num_gauss, num_poly)
            else:
                fit_inten = np.zeros_like(fit_wave)
        else:
            # Full image
            if num_wavelengths is None:
                fit_wave = self.fit['wavelength'][:,:,:]
                fit_inten = np.zeros((self.n_pxls, self.n_steps, self.n_wave))
            else:
                fit_wave = np.linspace(use_range[0], use_range[-1], num_wavelengths)
                fit_inten = np.zeros((self.n_pxls, self.n_steps, num_wavelengths))
            # Loop over locations and calculate each fit profile
            for ii in range(self.n_pxls):
                for jj in range(self.n_steps):
                    if self.fit['status'][ii, jj] > 0:
                        fit_inten[ii,jj,:] = self.fit_func(param_vals[ii,jj,:],
                                                           fit_wave[ii,jj,:],
                                                           num_gauss, num_poly)

        # Apply data masking as needed. Remember: 0 = False = good data (unmasked)
        if use_mask == True and num_wavelengths is None and 'mask' in self.fit.keys():
            if coords is not None:
                # Single spectrum
                mask_arr = self.fit['mask'][use_coords[0], use_coords[1], :]
                fit_wave = np.ma.array(fit_wave, mask=mask_arr)
                fit_inten = np.ma.array(fit_inten, mask=mask_arr)
            else:
                # Full image
                mask_arr = self.fit['mask'][:,:,:]
                fit_wave = np.ma.array(fit_wave, mask=mask_arr)
                fit_inten = np.ma.array(fit_inten, mask=mask_arr)

        return fit_wave, fit_inten

    def apply_radcal(self, input_radcal=None):
        """Apply a radiometric calibration curve (user-inputted or preflight)

        Parameters
        ----------
        input_radcal : array_like, optional
            User-inputted radiometric calibration curve. If set to None, will
            use the preflight radcal curve from the .meta dict. Default is None

        Returns
        -------
        output_fit : EISFitResult class instance
            A new EISFitResult class instance containing the calibrated params
        """
        if input_radcal is None:
            # Preflight radcal from HDF5 header file
            new_radcal = self.meta['radcal']
            radcal_wave = self.meta['wave']
        else:
            # User-inputted radcal curve
            new_radcal = np.array(input_radcal)
            radcal_wave = self.meta['wave'] #TODO: allow user to input wavelengths
            if len(new_radcal) != self.n_wave:
                print('Error: input_radcal must have the same number of elements'
                     +' as length of the wave dimension.')
                return self

        output_radcal = new_radcal
        if self.data_units not in ['ph', 'photon']:
            if str(self.radcal) == 'unknown':
                print('Error: Data currently has an unknown radcal applied.'
                     +' Unable to apply new calibration.')
                return self
            elif np.all(self.radcal == new_radcal):
                print('Error: input_radcal is identical to current radcal.'
                     +' No calculation is required.')
                return self
            else:
                print('Warning: Data currently has a different radcal applied.'
                     +' Old calibration curve will be removed.')
                new_radcal = new_radcal/self.radcal

        if str(radcal_wave) == 'unknown':
            print('Error: unknown wavelength locations of radcal curve.'
                 +' Unable to apply radcal curve.')
            return self

        # Make interpolated function of the radcal curve
        interp_radcal = interp1d(radcal_wave, new_radcal, kind='linear',
                                 fill_value='extrapolate')

        # Get interpolated radcal values for the peaks and constant background term
        loc_peaks = np.arange(self.n_gauss)*3
        loc_cen = np.arange(self.n_gauss)*3+1
        loc_const = self.n_gauss*3
        peak_radcal = interp_radcal(self.fit['params'][:,:,loc_cen])
        mean_waves = np.mean(self.fit['wavelength'], axis=2)
        const_radcal = interp_radcal(mean_waves)

        # Calculate new peak and background values
        new_peaks = self.fit['params'][:,:,loc_peaks]*peak_radcal
        new_const = self.fit['params'][:,:,loc_const]*const_radcal

        # Create a new EISFitResult object with the correct output
        output_res = EISFitResult(wave=self.fit['wavelength'], template=self.template,
                                  parinfo=self.parinfo, func_name=self.func_name,
                                  data_units='erg / (cm2 s sr)', radcal=output_radcal)

        new_fit = copy.deepcopy(self.fit)
        new_fit['param_units'] = output_res.fit['param_units']
        new_fit['params'][:,:,loc_peaks] = new_peaks
        new_fit['params'][:,:,loc_const] = new_const
        output_res.fit = new_fit
        output_res.meta = self.meta

        return output_res

    def remove_radcal(self):
        """Remove the applied radiometric calibration and convert data to counts

        Returns
        -------
        output_cube : EISCube class instance
            A new EISCube class instance containing the photon count data
        """
        if self.data_units in ['ph', 'photon']:
            print('Error: Data is already in units of photon counts.'
                 +' No calculation required.')
            return self
        elif str(self.radcal) == 'unknown':
            print('Error: Data currently has an unknown radcal applied.'
                 +' Unable to remove calibration.')
            return self

        radcal_wave = self.meta['wave']
        if str(radcal_wave) == 'unknown':
            print('Error: unknown wavelength locations of radcal curve.'
                 +' Unable to remove radcal curve.')
            return self

        # Make interpolated function of the radcal curve
        interp_radcal = interp1d(radcal_wave, self.radcal, kind='linear',
                                 fill_value='extrapolate')

        # Get interpolated radcal values for the peaks and constant background term
        loc_peaks = np.arange(self.n_gauss)*3
        loc_cen = np.arange(self.n_gauss)*3+1
        loc_const = self.n_gauss*3
        peak_radcal = interp_radcal(self.fit['params'][:,:,loc_cen])
        mean_waves = np.mean(self.fit['wavelength'], axis=2)
        const_radcal = interp_radcal(mean_waves)

        # Calculate new peak and background values
        new_peaks = self.fit['params'][:,:,loc_peaks]/peak_radcal
        new_const = self.fit['params'][:,:,loc_const]/const_radcal

        # Create a new EISFitResult object with the correct output
        output_res = EISFitResult(wave=self.fit['wavelength'], template=self.template,
                                  parinfo=self.parinfo, func_name=self.func_name,
                                  data_units='photon', radcal=None)

        new_fit = copy.deepcopy(self.fit)
        new_fit['param_units'] = output_res.fit['param_units']
        new_fit['params'][:,:,loc_peaks] = new_peaks
        new_fit['params'][:,:,loc_const] = new_const
        output_res.fit = new_fit
        output_res.meta = self.meta

        return output_res

    def rot_fov(self, end_time):
        # Return pointing information for the raster rotated to the input time.
        pointing = self.meta['pointing']
        xcen = pointing['xcen'] + pointing['offset_x']
        ycen = pointing['ycen'] + pointing['offset_y']
        new = rot_xy(xcen, ycen, pointing['ref_time'], end_time)
        fov = {'ref_time': end_time, 'xcen': new.Tx.value, 'ycen': new.Ty.value,
               'fovx': pointing['fovx'], 'fovy': pointing['fovy']}
        return fov

    def plot_fov(self, end_time, color='red', lw=1, ls='-'):
        # Return a patch of the raster FOV for plotting on an image.
        fov = self.rot_fov(end_time)
        x1 = fov['xcen'] - fov['fovx']/2
        y1 = fov['ycen'] - fov['fovy']/2
        rect = plt.Rectangle((x1, y1), fov['fovx'], fov['fovy'], fc='none', ec=color, lw=lw, ls=ls)
        return rect

    def get_aspect_ratio(self):
        # Aspect_ratio is given as y_scale/x_scale (NB: y_scale is 1 arcsec)
        return self.meta['aspect_ratio']

    def shift2wave(self, array, wave=195.119):
        # Shift an array from this fit to the desired wavelength
        this_wave = self.fit['wave_range'].mean()
        disp = np.zeros(len(array.shape))
        disp[0] = ccd_offset(wave) - ccd_offset(this_wave)
        array = shift_img(array, disp)
        return array

def create_fit_dict(n_pxls, n_steps, n_wave, n_gauss, n_poly, data_units='unknown'):
    """Dictionary to hold the fit parameters returned by fit_spectra()

    Parameters
    ----------
    n_pxls : int
        Number of pixels along each data slit
    n_steps : int
        Number of steps in the raster or sit-and-stare observation set
    n_wave : int
        number of wavelength points
    n_gauss :
        Number of Gaussian functions in the combinted fit profile
    n_poly :
        Degree of background polynomial

    Returns
    -------
    output : dict
        Empty dictionary with the correct dimensions and keys for a fit result.
    """

    n_param = 3*n_gauss + n_poly

    output = {'line_ids': np.zeros(n_gauss, dtype='<U32'),
              'main_component': 0,
              'n_gauss': n_gauss,
              'n_poly': n_poly,
              'wave_range': np.zeros(2),
              'status': np.zeros((n_pxls, n_steps)),
              'chi2': np.zeros((n_pxls, n_steps)),
              'mask': np.zeros((n_pxls, n_steps, n_wave), dtype='int'),
              'wavelength': np.zeros((n_pxls, n_steps, n_wave)),
              'int': np.zeros((n_pxls, n_steps, n_gauss)),
              'err_int': np.zeros((n_pxls, n_steps, n_gauss)),
              # 'peak': np.zeros((n_pxls, n_steps, n_gauss)),
              # 'err_peak': np.zeros((n_pxls, n_steps, n_gauss)),
              # 'centroid': np.zeros((n_pxls, n_steps, n_gauss)),
              # 'err_centroid': np.zeros((n_pxls, n_steps, n_gauss)),
              # 'width': np.zeros((n_pxls, n_steps, n_gauss)),
              # 'err_width': np.zeros((n_pxls, n_steps, n_gauss)),
              # 'background': np.zeros((n_pxls, n_steps, n_poly)),
              # 'err_background': np.zeros((n_pxls, n_steps, n_poly)),
              'params': np.zeros((n_pxls, n_steps, n_param)),
              'perror': np.zeros((n_pxls, n_steps, n_param)),
              'component': np.zeros(n_param, dtype='int'),
              'param_names': np.zeros(n_param, dtype='<U32'),
              'param_units': np.zeros(n_param, dtype='<U32')}

    # Create string labels for each component parameter
    # TODO: standardize name system to better match Astropy.modeling
    num_comp = n_gauss + 1 if n_poly > 0 else n_gauss
    for s in range(num_comp):
        if s == n_gauss and n_poly > 0:
            output['component'][3*n_gauss:] = s
            output['param_names'][3*n_gauss:] = ['c0', 'c1', 'c2', 'c3', 'c4'][0:n_poly]
            poly_AA = ['', ' /Angstrom ', ' /Angstrom^2', ' /Angstrom^3', ' /Angstrom^4']
            output['param_units'][3*n_gauss:] = [apu+data_units for apu in poly_AA][0:n_poly]
        else:
            output['component'][3*s:3*s+3] = s
            output['param_names'][3*s:3*s+3] = ['peak', 'centroid', 'width']
            output['param_units'][3*s:3*s+3] = [data_units, 'Angstrom', 'Angstrom']

    return output
