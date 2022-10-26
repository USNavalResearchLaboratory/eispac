"""
`~sunpy.map.Map` subclass for the EUV Imaging Spectrometer (EIS) on Hinode
"""
import pathlib

import numpy as np
import astropy.units as u
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
        # Check for a fits file containing a binary table with the errs
        # NB: errs are in a table with y-axis number of rows, each with x-axis
        #     number of values
        if isinstance(data, (str, pathlib.Path)):
            if pathlib.Path(data).suffix.lower().startswith('.fit'):
                import sunpy.io.fits
                hdu_list = sunpy.io.fits.read(data)
                data = hdu_list[0][0]
                header = hdu_list[0][1]
                if len(hdu_list) >= 2:
                    rec_errs = hdu_list[1][0]
                    ny = rec_errs.shape[0] # y-axis
                    nx = rec_errs[0][0].shape[0] # x-axis
                    errs = rec_errs['errors'].view(type=np.ndarray).reshape(ny,nx)
                    kwargs['uncertainty'] = errs

        super().__init__(data, header, **kwargs)

        # Setup plot settings
        self.plot_settings['aspect'] = self.meta['cdelt2'] / self.meta['cdelt1']
        # Adjust colormap and normalization depending on whether the map contains
        # intensity, velocity, or line width data
        if self.meta['measrmnt'].lower().startswith('int'):
            self.plot_settings['cmap'] = 'Blues_r'
            self.plot_settings['norm'] = ImageNormalize(stretch=AsinhStretch())
        elif self.meta['measrmnt'].lower().startswith('vel'):
            self.plot_settings['cmap'] = 'RdBu_r'
            # Autoscale color range to 3*std (rounded to nearest multiple of 5)
            vlim = 5*round(3*self.data.std()/5)
            self.plot_settings['norm'] = ImageNormalize(vmin=-vlim, vmax=vlim)
        elif self.meta['measrmnt'].lower().startswith('wid'):
            self.plot_settings['cmap'] = 'viridis'

    @property
    def spatial_units(self):
        units = self.meta.get('cunit1', 'arcsec'), self.meta.get('cunit2', 'arcsec')
        return SpatialPair(u.Unit(units[0]), u.Unit(units[1]))

    @property
    def processing_level(self):
        return self.meta.get('lvl_num', 3)

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
    def date(self):
        # NOTE: we override this property to prioritize date_average 
        # over DATE-OBS (or DATE_OBS). In GenericMap, this is reversed.
        # We do this because we want to make sure we are constructing our 
        # coordinate frames from DATE_AVG (the midpoint of the raster) and 
        # not DATE-OBS which is the beginning of the raster.
        time = self.date_average
        time = time or self._date_obs
        time = time or self.date_start
        time = time or self.date_end

        # This will fall back to the dafault behavior in GenericMap if the 
        # date property is still None, i.e. all other keys are still missing.
        time = time or super().date

        return time

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """
        Determines if header corresponds to an EIS image. Used to register
        EISMap with the sunpy.map.Map factory.
        """
        return str(header.get('instrume', '')).startswith('EIS')
