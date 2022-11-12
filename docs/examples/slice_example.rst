.. _ex-slice:

Slicing an EISCube
==================

Slicing an EISCube using real-world coordinates

.. code:: python

   >>> import eispac
   >>> import astropy.units as u
   >>> from astropy.coordinates import SkyCoord, SpectralCoord
   >>> from astropy.wcs.utils import wcs_to_celestial_frame

   >>> data_filename = 'eis_20190404_131513.data.h5'
   >>> data_cube = eispac.read_cube(data_filename, 195.119)
   >>> eis_frame = wcs_to_celestial_frame(data_cube.wcs)
   >>> lower_left = [SpectralCoord(195.0, unit=u.angstrom),
   ...               SkyCoord(Tx=48, Ty=225, unit=u.arcsec, frame=eis_frame)]
   >>> upper_right = [SpectralCoord(195.3, unit=u.AA),
   ...                SkyCoord(Tx=165, Ty=378, unit=u.arcsec, frame=eis_frame)]
   >>> data_cutout = data_cube.crop(lower_left, upper_right)

   >>> data_cube.dimensions
   [512, 87, 24] pix

   >>> data_cutout.dimensions
   [154, 48, 14] pix
