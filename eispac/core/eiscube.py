__all__ = ['EISCube']

import sys
import numpy as np
import astropy.units as u
from astropy.convolution import convolve, CustomKernel
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
        # slice all normal attributes
        kwargs = super()._slice(item)
        # The arguments for creating a new instance are saved in kwargs. So we
        # need to add additional keywords with our sliced, extra properties
        kwargs['wavelength'] = self.wavelength[item]
        kwargs['radcal'] = self.radcal
        return kwargs # these must be returned

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
        new_meta = self.meta.copy()
        new_meta['notes'].append('Applied radcal to convert photon counts to intensity')
        wcs_mask = (np.array(tuple(reversed(self.wcs.array_shape))) <= 1).tolist()

        output_cube = EISCube(new_data, wcs=self.wcs, uncertainty=new_errs,
                              wavelength=self.wavelength, radcal=output_radcal,
                              meta=new_meta, unit='erg / (cm2 s sr)',
                              mask=self.mask, missing_axes=wcs_mask)
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
        new_meta = self.meta.copy()
        new_meta['notes'].append('Removed radcal to convert intensity to photon counts')
        wcs_mask = (np.array(tuple(reversed(self.wcs.array_shape))) <= 1).tolist()

        output_cube = EISCube(new_data, wcs=self.wcs, uncertainty=new_errs,
                              wavelength=self.wavelength, radcal=None,
                              meta=new_meta, unit='photon',
                              mask=self.mask, missing_axes=wcs_mask)
        return output_cube

    def sum_spectra(self):
        """Sum the data along the spectral axis.

        Returns
        -------
        output_cube : NDCube class instance
            A new NDCube class instance containing the summed data
        """
        sum_data = np.sum(self.data, axis=2)
        new_wcs = self.wcs.dropaxis(0)
        new_meta = self.meta.copy()
        new_meta['notes'].append('Summed over the wavelength axis.')
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
        new_meta = self.meta.copy()
        new_meta['notes'].append('Smoothed using pixel widths of '+str(wid_list))
        wcs_mask = (np.array(tuple(reversed(self.wcs.array_shape))) <= 1).tolist()

        output_cube = EISCube(sm_data, wcs=self.wcs, uncertainty=sm_errs,
                              wavelength=self.wavelength, radcal=old_radcal,
                              meta=new_meta, unit=self.unit,
                              mask=sm_data_mask, missing_axes=wcs_mask)

        return output_cube
