__all__ = ['EISMap']

import pathlib
import numpy as np
import sunpy.map
from astropy.visualization import ImageNormalize, AsinhStretch, LinearStretch

"""
 A class to help sunpy maps properly plot EIS intensity fits files.

 v0.1 - 02-SEP-2021 : Original version, Will Barnes
"""

class EISMap(sunpy.map.GenericMap):
    """EIS fit parameter map.

    The EUV Imaging Spectrometer is part of the Hinode mission and was sponsored
    by the Japan Aerospace Exploration Agency (JAXA), the United Kingdom Space
    Agency (UKSA), and National Aeronautics and Space Administration (NASA) with
    contributions from ESA and Norway. Hinode was launched on September 22, 2006
    at 21:36 UTC from the Uchinoura Space Center in Japan and continues to
    operate. EIS observes two wavelength ranges in the extreme ultraviolet,
    171—212Å and 245—291Å with a spectral resolution of about 22mÅ and a plate
    scale of 100 per pixel.

    Notes
    -----
    Measurement errors are stored in a binary table. To load them correctly,
    you must pass the .fits file directly to eispac.EISMap instead of using
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
        self.meta['date-beg'] = self.meta['date_obs']
        self.meta['date-end'] = self.meta['date_end']

        # Want our coordinate frames to be constructed from DATE_AVG
        # and not DATE_OBS as this is the start of the observing run
        # del self.meta['date-obs'] # should we leave these for now, some routines expect them?
        # del self.meta['date_obs']
        # del self.meta['date_end']

        # Setup plot settings
        self.plot_settings['aspect'] = self.meta['CDELT2'] / self.meta['CDELT1']
        self.plot_settings['interpolation'] = 'kaiser'
        if self.meta['measrmnt'].lower().startswith('int'):
            self.plot_settings['cmap'] = 'Blues_r'
            self.plot_settings['norm'] = ImageNormalize(stretch=AsinhStretch())
        elif self.meta['measrmnt'].lower().startswith('vel'):
            self.plot_settings['cmap'] = 'RdBu_r'
            self.plot_settings['norm'] = ImageNormalize(stretch=LinearStretch())
        elif self.meta['measrmnt'].lower().startswith('wid'):
            self.plot_settings['cmap'] = 'viridis'
            self.plot_settings['norm'] = ImageNormalize(stretch=LinearStretch())

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an EIS image. Used to register
        EISMap with the sunpy.map.Map factory."""
        return str(header.get('instrume', '')).startswith('EIS')

    @property
    def observatory(self):
        return 'Hinode'

    @property
    def measurement(self):
        return self.meta.get('line_id', '')
