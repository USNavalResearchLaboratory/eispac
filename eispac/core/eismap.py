"""
`~sunpy.map.Map` subclass for the EUV Imaging Spectrometer (EIS) on Hinode
"""
import sys
import pathlib

import numpy as np
import astropy.units as u
from astropy.io import fits
from astropy.nddata import StdDevUncertainty
from astropy.visualization import ImageNormalize, AsinhStretch, LinearStretch
import sunpy.map
from sunpy.map.mapbase import SpatialPair
from sunpy.time import parse_time

__all__ = ['EISMap']


class EISMap(sunpy.map.GenericMap):
    """
    EIS fit parameter map.

    The EUV Imaging Spectrometer (EIS) is part of the Hinode mission and was sponsored
    by the Japan Aerospace Exploration Agency (JAXA), the United Kingdom Space
    Agency (UKSA), and National Aeronautics and Space Administration (NASA) with
    contributions from ESA and Norway. Hinode was launched on September 22, 2006
    at 21:36 UTC from the Uchinoura Space Center in Japan and continues to
    operate. EIS observes two wavelength ranges in the extreme ultraviolet,
    171—212 Å and 245—291 Å with a spectral resolution of about 22 mÅ and a plate
    scale of 100 per pixel.

    This data structure is designed to hold the fit parameters derived from multi-gaussian
    spectral fits to level 1, wavelength-resolved EIS rasters. These maps can contain the
    intensity, doppler velocity, or line width.

    Notes
    -----
    Measurement errors are stored in a binary table. To load them correctly,
    you must pass the .fits file directly to `eispac.EISMap` instead of using
    sunpy.map.Map

    References
    ----------
    * `Instrument Paper: Culhane, J. L., Harra, L. K., James, A. M., et al.
       2007, Sol. Phys., 243, 19`_
    """

    def __init__(self, data, header=None, **kwargs):
        step_date_obs = None
        step_exptime = None
    
        # Check for a fits file containing a binary table with the errs
        # NB: errs are in a table with y-axis number of rows, each with x-axis
        #     number of values
        if isinstance(data, (str, pathlib.Path)):
            if pathlib.Path(data).suffix.lower().startswith('.fit'):
                with fits.open(data, mode='readonly') as hdul:
                    data = hdul[0].data
                    header = hdul[0].header
                    if len(hdul) >= 2:
                        data_errors = StdDevUncertainty(hdul[1].data['errors'])
                        kwargs['uncertainty'] = data_errors
                    if len(hdul) >= 3:
                        step_date_obs = parse_time(hdul[2].data['step_date_obs'])
                        step_exptime = hdul[2].data['step_exptime']

        # Get user-input step_date_obs and step_exptime
        # NB: this will overwrite any values from an input fits file
        user_step_date_obs = kwargs.pop('step_date_obs', None)
        user_step_exptime = kwargs.pop('step_exptime', None)
        if user_step_date_obs is not None and header is not None:
            if len(user_step_date_obs) == header['naxis1']:
                step_date_obs = parse_time(user_step_date_obs)
            else:
                print(f'WARNING incorrect number of "step_date_obs" values!'
                     +f' This EIS observation has {header["naxis1"]} steps.',
                      file=sys.stderr)
        if user_step_exptime is not None and header is not None:
            if len(user_step_exptime) == header['naxis1']:
                step_exptime = user_step_exptime
            else:
                print(f'WARNING: incorrect number of "step_exptime" values!'
                     +f' This EIS observation has {header["naxis1"]} steps.',
                      file=sys.stderr)

        # Initalize the map
        super().__init__(data, header, **kwargs)
        self._step_date_obs = step_date_obs
        self._step_exptime = step_exptime

        # Setup plot settings and get default data masks
        # This includes adjusting colormap and normalization depending on whether 
        # the map contains intensity, velocity, or line width data
        default_mask = None
        self.plot_settings['aspect'] = self.meta['cdelt2'] / self.meta['cdelt1']
        if self.meta['measrmnt'].lower().startswith('int'):
            self.plot_settings['cmap'] = 'Blues_r'
            self.plot_settings['norm'] = ImageNormalize(stretch=AsinhStretch())
            default_mask = self.data == 0
        elif self.meta['measrmnt'].lower().startswith('vel'):
            self.plot_settings['cmap'] = 'RdBu_r'
            # Autoscale color range to 3*std (rounded to nearest multiple of 5)
            vlim = 5*round(3*self.data.std()/5)
            self.plot_settings['norm'] = ImageNormalize(vmin=-vlim, vmax=vlim)
            if self.uncertainty is not None:
                # Note: velocities of 0 may be valid UNLESS the errors are NaN
                default_mask = (self.data == 0) & np.isnan(self.uncertainty.array)
        elif self.meta['measrmnt'].lower().startswith('wid'):
            self.plot_settings['cmap'] = 'viridis'
            default_mask = self.data == 0

        # Set the default mask (ignored if the user input their own mask)
        if self.mask is None:
            self.mask = default_mask

    @property
    def spatial_units(self):
        units = self.meta.get('cunit1', 'arcsec'), self.meta.get('cunit2', 'arcsec')
        return SpatialPair(u.Unit(units[0]), u.Unit(units[1]))

    @property
    def processing_level(self):
        return self.meta.get('lvl_num', 2)

    @property
    def waveunit(self):
        return u.Unit(self.meta.get('waveunit', 'angstrom'))

    @property
    def wavelength(self):
        line_id = self.meta.get('line_id')
        if line_id is not None:
            wave = float(line_id.split()[-1])
            return u.Quantity(wave, self.waveunit)

    @property
    def measurement(self):
        return self.meta.get('measrmnt', '')

    @property
    def observatory(self):
        return 'Hinode'

    @property
    def nickname(self):
        line_id = self.meta.get('line_id', '')
        return f'{self.observatory} {self.instrument} {line_id}'

    @property
    def date_start(self):
        # Try default key DATE-BEG. This is to future proof against
        # switching to DATE-BEG when constructing the L1 headers
        # NOTE: the DATE_OBS key is the beginning of the observation
        # so we can use this in case DATE_BEG is missing
        date_beg = self._get_date('date_beg') or super().date_start
        date_beg = date_beg or self._date_obs
        return date_beg

    @property
    def date_end(self):
        # Try default key DATE-END. This is to future proof against
        # switching to DATE-END when constructing the L1 headers
        return self._get_date('date_end') or super().date_end

    @property
    def date_average(self):
        return self._get_date('date_avg') or super().date_average

    @property
    def reference_date(self):
        """
        The reference date for the coordinate system

        According to Section 2.4 of
        `EIS Software Note 9 <https://solarb.mssl.ucl.ac.uk/SolarB/eis_docs/eis_notes/09_POINTING/eis_swnote_09.pdf>`_,
        the pointing keywords are defined at the start of the raster. As such, this property returns the
        time at the beginning of the raster.

        .. note:: This property is overridden because `sunpy.map.GenericMap` sets
                  this to be the ``.date_average`` which in this case is the midpoint
                  of the raster.
        """
        return self.date_start

    @property
    def duration(self):
        """Total duration of the observation in units of [min]
        """
        total_time = (np.datetime64(self.meta['date_end']) 
                    - np.datetime64(self.meta['date_obs']))
        total_time = total_time / np.timedelta64(60, 's') # convert to [min] 
        return total_time * u.Unit('min')
    
    @property
    def step_date_obs(self):
        """date_obs timestamp of each step along the x-axis
        """
        if self._step_date_obs is not None:
            return self._step_date_obs
        else:
            print(f'WARNING: exact "step_date_obs" values are unknown!' 
                 +f' Estimating based on the observation start and end times.',
                  file=sys.stderr)
            total_time = self.duration.to('s').value
            est_cad = total_time / self.meta['naxis1']
            est_date_obs = (np.datetime64(self.meta['date_obs']) 
                           + np.arange(self.meta['naxis1'])
                           * np.timedelta64(int(est_cad*1000), 'ms'))

            if self.meta['nraster'] == 1: 
                # Sit-and-stare timestamps inc left to right
                return parse_time(est_date_obs)
            else:
                # Normal raster timestamps inc from right to left (scan dir)
                return parse_time(np.flip(est_date_obs))
        
    @property
    def step_exptime(self):
        """Exposure time of each step along the x-axis
        """
        if self._step_exptime is not None:
            return self._step_exptime
        else:
            print(f'WARNING: exact "step_exptime" values are unknown!' 
                 +f' Estimating based on the observation start and end times.'
                 +f' Actual exposure times will be shorter due to on-board'
                 +f' processing and the specific observation plan.',
                  file=sys.stderr)
            total_time = self.duration.to('s').value
            est_avg_exptime = total_time / self.meta['naxis1']
            return np.zeros(self.meta['naxis1']) + est_avg_exptime
        
    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """
        Determines if header corresponds to an EIS image. Used to register
        EISMap with the sunpy.map.Map factory.
        """
        return str(header.get('instrume', '')).startswith('EIS')
