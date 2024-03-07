.. _sec-quick:

A Quick Guide to EISPAC
=======================

Once installed, EISPAC can be imported into any Python script or
interactive session with a simple ``import eispac`` statement.

The three most important functions are `~eispac.core.read_cube`,
`~eispac.core.read_template`, and `~eispac.core.fit_spectra`.
Details concerning their usage can be found below. Please also see :ref:`ex-advanced`
program for a more complete demonstration.

read_cube
---------

Example call signature,

::

   data_cube = read_cube('eis_20190404_131513.data.h5', 195.119)

`~eispac.core.read_cube` typically requires two arguments:

1. **filename** (str or `~pathlib.Path`)

   Name or path of either the data or head HDF5 file for a single EIS observation

2. **window** (int or float, optional)

   Requested spectral window number or the value of any wavelength within the
   requested window. Default is '0'

The return value is an `~eispac.core.EISCube` class instance which contains
calibrated intensities, corrected wavelengths, and all of the associated metadata.
``EISCube`` objects are a subclass of `~ndcube.NDCube` from the Sunpy-affiliated
package of the same name.

You can slice an ``NDCube`` object using either array indices or by using physical
coordinates and the `ndcube.NDCube.crop` method. Please see the `NDCube documentation
<https://docs.sunpy.org/projects/ndcube/en/stable/index.html>`_ for more details.

The `~eispac.core.EISCube` subclass extends `~ndcube.NDCube` and provides
additional features. First, an extra ``.wavelength`` attribute has been added
which contains a 3D array with the corrected wavelength values at all locations
within the cube. Slicing an ``EISCube`` will also appropriately slice the wavelength
array. Secondly, there are a few extra methods for your convenience, such as
``.sum_spectra()`` which returns a new, 2D ``NDCube`` containing the sum along
the wavelength axis.

A few more notes about `EISCube` objects,

- Data axes are in the same order as the array stored in the HDF5 file. Namely, the
  order is (slit_pixel, raster_step, wavelength), i.e. (solar-y, solar-x, wavelength).

- All of the EIS metadata and header information are stored as separate keys in
  the ``.meta`` dictionary attribute (e.g. the original FITS header is in
  ``.meta['index']`` while the updated pointing information is in
  ``.meta['pointing']``).

read_template
-------------

Example call signature,

::

   tmplt = read_template('fe_12_195_119.2c.template.h5')

`~eispac.core.read_template` is relatively simple. It takes a single argument
giving the filename of a template file and returns an `~eispac.core.EISFitTemplate`
class instance containing the the initial fit parameters. You can also load custom 
templates by first defining a TOML file with your template parameters (see the
:ref:`Custom Fit Templates <sec-template_toml>`: section for an example TOML
file). You may print a summary of the template and parameter contraints by typing
``print(TEMPLATE)``. For convenience, there is also a ``.central_wave`` attribute 
that contains the mean wavelength value within the template. This can be useful 
for loading the correct spectral window using `~eispac.core.read_cube`.

fit_spectra
-----------

This is the main fitting routine. Example call signature,

::

   fit_results = fit_spectra(data_cube, fit_template)

The simplest way to use `~eispac.core.fit_spectra()` is to give it two arguments:

1. **data_cube** (EISCube object or filepath)

   An EISCube class instance (or slice) containing one or more intensity profiles
   to be fit. Wavelength and error values will be extracted as needed from the EISCube.

2. **template** (EISFitTemplate object or filepath)

   An EISFitTemplate class instance containing both the general template information
   (in the ``.template`` attribute) and a ``.parinfo`` attribute with the fit parameter
   dictionary.

Alternatively, you may provide the function with individual data arrays like this,

::

   fit_result = fit_spectra(inten_arr, template_dict, parinfo=parinfo_dict,
                            wave=wavelength_arr, errs=error_arr)

The function will loop over the data according to its dimensionality. 3D data is
assumed to be a full EIS raster (or a sub region), 2D data is assumed to be a
single EIS slit, and 1D data is assumed to be a single profile.

`~eispac.core.fit_spectra` returns a `~eispac.core.EISFitResult` class instance
containing the fit parameter values. As well as assorted metadata, there are two
important class methods that may be of use ``.get_params()`` and ``.get_fit_profile()``.

The method ``.get_params()`` extracts parameters values by either component number,
name, or pixel coordinate (or any combination of the three). The arguments are:

- **component** (int or None, optional)

  Integer number (or list of ints) of the functional component(s). If set to None,
  will return the total combined fit profile. Default is None.

- **param_name** (str, optional)

  String name of the requested parameter. If set to None, will not filter based
  on paramater name. Default is None.

- **coords** (list or tuple, optional)

  (Y, X) coordinates of the requested datapoint. If set to None, will instead
  return the parameters at all locations. Default is None

- **num_wavelengths** (int, optional)

  Number of wavelength values to compute the fit intensity at. These values will
  be equally spaced and span the entire fit window. If set to None, will use the
  observed wavelength values. Default is None.

- **casefold** (bool, optional)

  If set to True, will ignore case when extracting parameters by name. Default is False.

Examples,

::

   c0_params = fit_res.get_fit_profile(component=0)
   widths = fit_res.get_fit_profile(param_name='width')

The ``.get_fit_profile()`` method may be used to generate the either the combined fit
intensity profile or the profile of a single component function. The method takes
up to three arguments:

- **component** (int or None, optional)

  Integer number (or list of ints) of the functional component(s). If set to None,
  will return the total combined fit profile. Default is None.

- **coords** (list or tuple, optional)

  (Y, X) coordinates of the requested datapoint. If set to None, will instead
  return the parameters at all locations. Default is None

- **num_wavelengths** (int, optional)

  Number of wavelength values to compute the fit intensity at. These values will
  be equally spaced and span the entire fit window. If set to None, will use the
  observed wavelength values. Default is None.

``get_fit_profile()`` returns two arrays, ``fit_wave`` & ``fit_inten``, which
contain the wavelengths and corresponding fit intensity values.

Examples,

::

   fit_x, fit_y = fit_res.get_fit_profile(coords=(5,5))
   c0_fit_x, c0_fit_y = fit_res.get_fit_profile(component=0, num_wavelengths=100)
