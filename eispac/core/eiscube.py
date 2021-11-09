__all__ = ['EISCube']

import sys
import copy
import numpy as np
import astropy.units as u
from astropy.convolution import convolve, CustomKernel
from astropy.coordinates import SkyCoord
from ndcube import __version__ as ndcube_ver
from ndcube import NDCube

class EISCube(NDCube):
    # NOTE: this is based on the example given at
    # https://docs.astropy.org/en/stable/nddata/subclassing.html#slicing-an-additional-property
    def __init__(self, *args, **kwargs):
        # Extract extra attributes, if given, and initialize them correctly.
        input_wave = kwargs.pop('wavelength') if 'wavelength' in kwargs else None
        if input_wave is None:
            input_wave = np.zeros_like(args[0], dtype=float)
        self.wavelength = input_wave
        input_radcal = kwargs.pop('radcal') if 'radcal' in kwargs else 'unknown'
        self._current_radcal = input_radcal
        kwargs['copy'] = True # Try to ensure meta is not leaked between cutouts
        super().__init__(*args, **kwargs)

    @property
    def wavelength(self):
        return self._wavelength

    @wavelength.setter
    def wavelength(self, input_array):
        self._wavelength = input_array

    @property
    def radcal(self):
        return self._current_radcal

    @radcal.setter
    def radcal(self, input_array):
        print('Error: Please use the .apply_radcal() and .remove_radcal()'
             +' methods to modify or change the radiometric calibration.')

    def _slice(self, item):
        kwargs = super()._slice(item) # slice all normal attributes
        old_meta = kwargs.pop('meta')
        kwargs['meta'] = copy.deepcopy(old_meta) # no refs, please!
        # The arguments for creating a new instance are saved in kwargs. So we
        # need to add additional keywords with our sliced, extra properties
        kwargs['wavelength'] = self.wavelength[item]
        kwargs['radcal'] = self.radcal

        # Update the 'mod_index' (used for exporting to .fits after fitting)
        # Reminder: 'mod_index' uses a fits image axes order of [X, Y, Wave]
        m_idx = copy.deepcopy(kwargs['meta']['mod_index'])
        ll_wcs = kwargs['wcs'].low_level_wcs
        wcs_arr_shape = ll_wcs.array_shape # axis order of [Y, X, wave]
        new_ori = kwargs['wcs'].array_index_to_world(0,0,0)
        ax_shape = [1, 1, 1] # true length of [Y, X, Wave] axes

        if isinstance(new_ori, list):
            # 3D subcube or 2D [spatial, wave] slice (one axis removed)
            x1 = new_ori[-1].Tx.to('arcsec').value
            y1 = new_ori[-1].Ty.to('arcsec').value
            w1 = new_ori[0].to('Angstrom').value
        elif isinstance(new_ori, SkyCoord):
            # 2D slice at a given wavelength value or 1D spatial profile
            x1 = new_ori.Tx.to('arcsec').value
            y1 = new_ori.Ty.to('arcsec').value
            lost_unit = ll_wcs.dropped_world_dimensions['world_axis_units']
            lost_value = ll_wcs.dropped_world_dimensions['value']
            w1 = u.Quantity(lost_value[0], lost_unit[0]).to('Angstrom').value
        elif isinstance(new_ori, u.Quantity):
            # Single spectrum at a selected location
            w1 = new_ori.to('Angstrom').value
            lost_units = ll_wcs.dropped_world_dimensions['world_axis_units']
            lost_values = ll_wcs.dropped_world_dimensions['value']
            x1 = u.Quantity(lost_values[0], lost_units[0]).to('arcsec').value
            y1 = u.Quantity(lost_values[1], lost_units[1]).to('arcsec').value

        # ndcube >= 2.0 drops all shallow (length 1) axes from .array_shape
        # UNLESS an axis was sliced with explicit start:stop values of i:i+1
        # Unfortunatly, this means there is no reliable method to identify the
        # true length of each axis and which axis was dropped (if any).
        # Therefore, we must resort to checking the type of each input
        # slice parameter in "item" and try to manaully match axes with lengths
        if ndcube_ver >= '2.0.0':
            # First, extract and rearrange old shape in order of [Y, X, Wave]
            old_shape = (m_idx['naxis2'], m_idx['naxis1'], m_idx['naxis3'])
            ax_i = 0 # true axis index
            wcs_i = 0 # sliced wcs axis index
            # Note: we must take care when slicing a 2D slice of a cube
            for s_i in range(len(item)):
                while ax_i <= 2:
                    if old_shape[ax_i] == 1:
                        # skip previously dropped axes
                        ax_i += 1
                    else:
                        if isinstance(item[s_i], slice):
                            ax_shape[ax_i] = wcs_arr_shape[wcs_i]
                            wcs_i += 1
                        else:
                            ax_shape[ax_i] = 1
                        ax_i += 1
                        break
        else:
            # Works just fine in ndcube 1.4.2 (for all kinds of slices)
            ax_shape = wcs_arr_shape

        x2 = x1 + ax_shape[1]*m_idx['cdelt1']
        y2 = y1 + ax_shape[0]*m_idx['cdelt2']
        x_shift = round(abs(x1 - m_idx['crval1'])/m_idx['cdelt1'])
        y_shift = round(abs(y1 - m_idx['crval2'])/m_idx['cdelt2'])
        w_shift = round(abs(w1 - m_idx['crval3'])/m_idx['cdelt3'])
        m_idx['naxis1'] = ax_shape[1] # X-axis
        m_idx['naxis2'] = ax_shape[0] # Y-axis
        m_idx['naxis3'] = ax_shape[2] # Wavelength axis
        m_idx['crpix1'] = 1.0 - x_shift
        m_idx['crpix2'] = 1.0 - y_shift
        m_idx['crpix3'] = 1.0 - w_shift
        m_idx['fovx'] = x2 - x1
        m_idx['fovy'] = y2 - y1
        m_idx['xcen'] = x1 + 0.5*(x2-x1)
        m_idx['ycen'] = y1 + 0.5*(y2-y1)
        kwargs['meta']['mod_index'] = m_idx
        kwargs['meta']['extent_arcsec'] = [x1, x2, y1, y2] # [L, R, Top, Bot]

        return kwargs # these must be returned

    def crop_by_coords(self, *args, **kwargs):
        """REMOVED IN NDCube 2.0"""
        print('Error: crop_by_coords() was removed in NDCube 2.0. Please use'
             +' the .crop() or .crop_by_values() methods instead. See the'
             +' NDCube documentation for more information.', file=sys.stderr)
        return None

    def apply_radcal(self, input_radcal=None):
        """Apply a radiometric calibration curve (user-inputted or preflight)

        Parameters
        ----------
        input_radcal : array_like, optional
            User-inputted radiometric calibration curve. If set to None, will
            use the preflight radcal curve from the .meta dict. Default is None

        Returns
        -------
        output_cube : EISCube class instance
            A new EISCube class instance containing the calibrated data
        """
        if input_radcal is None:
            # Preflight radcal from HDF5 header file
            new_radcal = self.meta['radcal']
        else:
            # User-inputted radcal curve
            new_radcal = np.array(input_radcal)
            if len(new_radcal) != self.data.shape[-1]:
                print('Error: input_radcal must have the same number of elements'
                     +' as the last dimension in the data array.')
                return self

        output_radcal = new_radcal
        if self.unit != u.photon:
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

        new_data = self.data.copy()*new_radcal
        new_errs = self.uncertainty.array.copy()*new_radcal
        new_meta = copy.deepcopy(self.meta)
        new_meta['mod_index']['bunit'] = 'erg / (cm2 s sr)'
        new_meta['notes'].append('Applied radcal to convert photon counts to intensity')
        # wcs_mask = (np.array(tuple(reversed(self.wcs.array_shape))) <= 1).tolist()

        output_cube = EISCube(new_data, wcs=self.wcs, uncertainty=new_errs,
                              wavelength=self.wavelength, radcal=output_radcal,
                              meta=new_meta, unit='erg / (cm2 s sr)',
                              # mask=self.mask, missing_axes=wcs_mask)
                              mask=self.mask)

        return output_cube

    def remove_radcal(self):
        """Remove the applied radiometric calibration and convert data to counts

        Returns
        -------
        output_cube : EISCube class instance
            A new EISCube class instance containing the photon count data
        """
        if self.unit == u.photon:
            print('Error: Data is already in units of photon counts.'
                 +' No calculation required.')
            return self
        elif str(self.radcal) == 'unknown':
            print('Error: Data currently has an unknown radcal applied.'
                 +' Unable to remove calibration.')
            return self

        new_data = self.data.copy()/self.radcal
        new_errs = self.uncertainty.array.copy()/self.radcal
        new_meta = copy.deepcopy(self.meta)
        new_meta['mod_index']['bunit'] = 'photon'
        new_meta['notes'].append('Removed radcal to convert intensity to photon counts')
        # wcs_mask = (np.array(tuple(reversed(self.wcs.array_shape))) <= 1).tolist()

        output_cube = EISCube(new_data, wcs=self.wcs, uncertainty=new_errs,
                              wavelength=self.wavelength, radcal=None,
                              meta=new_meta, unit='photon',
                              # mask=self.mask, missing_axes=wcs_mask)
                              mask=self.mask)

        return output_cube

    def sum_spectra(self, wave_range=None, units=u.Angstrom):
        """Sum the data along the spectral axis.

        Parameters
        ----------
        wave_range : list of ints, floats, or Quantity instances
            Wavelength range to sum over. Values can be input as either
            [min, max] or [center, half width]. Units can be specified using
            either Astropy units instances or by inputting a pair of ints or
            floats and then also using the "units" keyword. If wave_range is set
            to None, then entire spectra will be summed over. Default is None.
        units : str or Quantity instance
            Units to be used for the wavelength range if wave_range is given a
            list of ints or floats. Will be ignored if either wave_range is None
            or is given a list with Astropy units. Default is 'Angstrom'.

        Returns
        -------
        output_cube : NDCube class instance
            A new NDCube class instance containing the summed data
        """
        if wave_range is None:
            # Sum over entire wavelength axis and return an NDCube
            try:
                new_wcs = self.wcs.dropaxis(0)
            except:
                new_wcs = copy.deepcopy(self[:,:,0].wcs)
            sum_data = np.sum(self.data, axis=2)
            new_meta = copy.deepcopy(self.meta)
            new_meta['notes'].append('Summed over entire wavelength axis.')
            return NDCube(sum_data, new_wcs, meta=new_meta)

        # Validate input wavelength range
        if isinstance(wave_range, (list, tuple)):
            use_range = [0, 0]
            range_units = ['unknown', 'unknown']
            print('Summing EISCube spectra over a select wavelength range.')
            if len(wave_range) != 2:
                print('Error: invalid number of wave_range values. Please input'
                     +' a list or tuple with exactly two elements.',
                     file=sys.stderr)
                return None
        else:
            print('Error: invalid wave_range type. Please input either None or'
                 +' a list (or tuple) with two elements.', file=sys.stderr)
            return None

        for w in range(2):
            if isinstance(wave_range[w], u.Quantity):
                # Parse an astropy.units.Quantity and convert as needed
                # Note: this will overwrite any inputs to the "units" kwarg
                if wave_range[w].unit == u.pix:
                    use_range[w] = wave_range[w].value
                    range_units[w] = u.pix
                elif wave_range[w].unit.physical_type == 'length':
                    use_range[w] = wave_range[w].to('Angstrom').value
                    range_units[w] = u.Angstrom
                else:
                    print('Error: invalid wavelength unit. Please input a pixel'
                    +' or length unit.', file=sys.stderr)
                    return None
            else:
                # Assume default or user inputted units (still convert if needed)
                input_units = u.Unit(units)
                if input_units == u.pix:
                    use_range[w] = float(wave_range[w])
                    range_units[w] = u.pix
                elif input_units.physical_type == 'length':
                    u_scale = input_units.to('Angstrom')
                    use_range[w] = float(wave_range[w])*u_scale
                    range_units[w] = u.Angstrom
                else:
                    print('Error: invalid wavelength unit. Please input a pixel'
                         +' or length unit.', file=sys.stderr)
                    return None

        # Check for consistent units
        if range_units[0] != range_units[1]:
            print('Error: mismatched units. Please input the same units for'
                 +' both wave_range elements or use the "units" keyword',
                 file=sys.stderr)
            return None

        # If given values of [center, half width], compute the actual range
        if use_range[1] < use_range[0]:
            temp_center = use_range[0]
            temp_half_wid = use_range[1]
            use_range[0] = temp_center - temp_half_wid
            use_range[1] = temp_center + temp_half_wid

        # Get indices to be summed over
        w_indices = [0, -1]
        if range_units[0] == u.pix:
            # Round pixels values to nearest whole indice
            w_indices[w] = int(round(use_range[w]))
        elif range_units[0] == u.Angstrom:
            # Find the closest pixel location on the average wavelength axis
            try:
                # Note: the corrected wavelength has units of [Angstrom]
                w_coords = np.mean(self.wavelength, axis=(0,1))
            except KeyError:
                print('Error: missing or invalid corrected wavelength array.')
                return None
            for w in range(2):
                abs_w_diff = np.abs(w_coords - use_range[w])
                w_indices[w] = np.argmin(abs_w_diff)

        try:
            new_wcs = self.wcs.dropaxis(0)
        except:
            new_wcs = copy.deepcopy(self[:,:,0].wcs)
        sum_data = np.sum(self.data[:,:,w_indices[0]:w_indices[1]+1], axis=2)
        new_meta = copy.deepcopy(self.meta)
        new_meta['notes'].append('Summed wavelength axis over the range of '
                                +str(use_range)+' '+str(range_units[0]))
        return NDCube(sum_data, new_wcs, meta=new_meta)

    def smooth_cube(self, width=3, **kwargs):
        """Smooth the data along one or more spatial axes.

        Parameters
        ----------
        width : list or single value of ints, floats, or Quantity instances
            Number of pixels or angular distance to smooth over. If given a
            single value, only the y-axis will be smoothed. Floats and angular
            distances will be converted to the nearest whole pixel value.
            If a width value is even, width + 1 will be used instead.
            Default is width = 3
        **kwargs : keywords or dict
            Keyword arguments to be passed to the astropy.convolution.convolve()
            function.

        Returns
        -------
        output_cube : EISCube class instance
            A new EISCube class instance containing the smoothed data
        """
        # Validate input width
        num_dims = len(self.dimensions)
        wid_list = [1]*num_dims # NB: a width of 1 results in no smoothing
        if isinstance(width, (list, tuple)):
            # Note: we assume the last dim is always wavelength
            wid_list[0] = width[0]
            if num_dims > 2:
                wid_list[1] = width[1]
                print('Warning: smoothing over the x-axis can yield unexpected'
                     +' results due to the time interval between observations.'
                     +' Use with care.')

            if len(width) >= num_dims:
                print('Warning: smoothing over the wavelength axis is not'
                     +' supported. Only widths for the Y & X axes will be used')
        elif isinstance(width, (int, float, u.Quantity)):
            wid_list[0] = width # Only smooth along y-axis
        else:
            print('Error: invalid width data type. Please input an int, float,'
                 +' or astropy.units.Quantity instance', file=sys.stderr)
            return None

        coord_ax = ['y', 'x', 'w']
        for w in range(len(wid_list)-1):
            # Parse a astropy.units.Quantity and convert to units of pixels
            if isinstance(wid_list[w], u.Quantity):
                if wid_list[w].unit == u.pix:
                    wid_list[w] = wid_list[w].value
                elif not wid_list[w].unit.physical_type == 'angle':
                    print('Error: invalid width unit. Please input a pixel or'
                         +' angular unit.', file=sys.stderr)
                    return None
                else:
                    try:
                        # Note: y & x scales are in units of [arcsec]/[pixel]
                        ax_scale = self.meta['pointing'][coord_ax[w]+'_scale']
                    except KeyError:
                        print('Error: missing '+coord_ax[w]+'-axis scale.')
                        return None
                    angular_wid_str = str(wid_list[w])
                    wid_list[w] = wid_list[w].to('arcsec').value / ax_scale
                    print('Note: on the '+coord_ax[w]+'-axis, '+angular_wid_str
                         +' is equivalent to '+str(wid_list[w])+' pixels.')

            # Round to nearest pixel and add 1 to even values
            wid_list[w] = int(round(wid_list[w]))
            if wid_list[w] % 2 == 0:
                wid_list[w] = wid_list[w] + 1

        # Create smoothing kernel with normalized weights (i.e. sum to 1)
        # Note: Using a 2D or 3D kernel allows us to smooth everything at once
        sm_weights = np.ones(wid_list) / (wid_list[0]*wid_list[1])
        sm_kernel = CustomKernel(sm_weights)

        # Calculate smoothed data and uncertainty values
        sm_data = convolve(self.data, sm_kernel, **kwargs)
        if self.uncertainty is not None:
            sm_errs = np.sqrt(convolve(self.uncertainty.array**2,
                                       sm_kernel, **kwargs))
        else:
            sm_errs = none
        sm_data_mask = np.logical_or(np.isnan(sm_data), sm_data < 0)

        # Pack everything up in a new EISCube
        old_radcal = self.radcal
        new_meta = copy.deepcopy(self.meta)
        new_meta['notes'].append('Smoothed using pixel widths of '+str(wid_list))
        # wcs_mask = (np.array(tuple(reversed(self.wcs.array_shape))) <= 1).tolist()

        output_cube = EISCube(sm_data, wcs=self.wcs, uncertainty=sm_errs,
                              wavelength=self.wavelength, radcal=old_radcal,
                              meta=new_meta, unit=self.unit,
                              # mask=sm_data_mask, missing_axes=wcs_mask)
                              mask=sm_data_mask)

        return output_cube
