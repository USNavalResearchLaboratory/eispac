"""
`~sunpy.map.Map` subclass for the EUV Imaging Spectrometer (EIS) on Hinode
"""
import pathlib

import numpy as np
import astropy.units as u
from astropy.visualization import ImageNormalize, AsinhStretch, LinearStretch
import sunpy.map
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
    171—212Å and 245—291Å with a spectral resolution of about 22mÅ and a plate
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

        # Validate important keywords and add them in, if missing
        self.meta['ctype1'] = self.meta.get('ctype1', 'HPLN-TAN')
        self.meta['ctype2'] = self.meta.get('ctype2', 'HPLT-TAN')
        self.meta['cunit1'] = self.meta.get('cunit1', 'arcsec')
        self.meta['cunit2'] = self.meta.get('cunit2', 'arcsec')
        self.meta['waveunit'] = self.meta.get('waveunit', 'angstrom')
        self.meta['date-beg'] = self.meta.get('date_beg', self.meta.get('date_obs'))
        self.meta['date-end'] = self.meta.get('date_end')
        self.meta['date-obs'] = self.meta.get('date_obs')
        self.meta['date-avg'] = self.meta.get('date_avg')
        # NOTE: this block can be removed once sunpy>=3.1 is only supported as
        # the .date_average property will always be constructed in this way if
        # date_start and date_end are present
        if self.meta['date-avg'] is None:
            if self.meta['date-beg'] is not None and self.meta['date-end'] is not None:
                timesys = self.meta.get('timesys', 'UTC').lower()
                start = parse_time(self.meta['date-beg'], scale=timesys)
                end = parse_time(self.meta['date-end'], scale=timesys)
                self.meta['date-avg'] = (start + (end - start)/2).isot

        # Setup plot settings
        self.plot_settings['aspect'] = self.meta['CDELT2'] / self.meta['CDELT1']
        if self.meta['measrmnt'].lower().startswith('int'):
            self.plot_settings['cmap'] = 'Blues_r'
            self.plot_settings['norm'] = ImageNormalize(stretch=AsinhStretch())
        elif self.meta['measrmnt'].lower().startswith('vel'):
            self.plot_settings['cmap'] = 'RdBu_r'
            self.plot_settings['norm'] = ImageNormalize(vmin=-40.0, vmax=40.0)
        elif self.meta['measrmnt'].lower().startswith('wid'):
            self.plot_settings['cmap'] = 'viridis'

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """
        Determines if header corresponds to an EIS image. Used to register
        EISMap with the sunpy.map.Map factory.
        """
        return str(header.get('instrume', '')).startswith('EIS')

    @property
    def observatory(self):
        return 'Hinode'

    @property
    def measurement(self):
        line_id = self.meta.get('line_id', '')
        quantity = self.meta.get('measrmnt', '')
        return  f'{line_id} {quantity}'

    @property
    def wavelength(self):
        line_id = self.meta.get('line_id')
        if line_id is not None:
            wave = float(line_id.split()[-1])
            return u.Quantity(wave, self.meta['waveunit'])

    @property
    def date(self):
        # Want to make sure we are constructing our coordinate frames
        # from DATE_AVG and not the DATE-OBS which is the beginning of
        # the observation
        # In sunpy 3.1, this may not be needed as we could just remove
        # date-obs and then .date will default to .date_average
        t = self.meta.get('date-avg')
        timesys = self.meta.get('timesys', 'UTC')
        if t is None:
            return super().date
        else:
            return parse_time(t, scale=timesys.lower())

    @property
    def processing_level(self):
        return self.meta.get('lvl_num', 3)
