# A Quick and Dirty Guide to the Updated eispac Software
_version 0.42_

Once installed, the package can be imporated as you would expect using
`import eispac`

There are currently only three main scripts you will need to use - `read_cube`,
`read_template`, and `fit_spectra`. Details concerning their usage can be found
below. Please also see the example program `eis_fit_cube_example.py`.

## read_cube
Example call signature,
```
read_cube('eis_20190404_131513.data.h5', 195.119)
```

`read_cube()` takes two arguments:
1. _filename_ (str or pathlib path) - Name or path of either the data or head
   HDF5 file for a single EIS observation
2. _window_ (int or float, optional) - Requested spectral window number or the
   value of any wavelength within the requested window. Default is '0'

The return value is an `EISCube` class instance which contains calibrated intensities,
corrected wavelengths, and all of the associated metadata. `EISCube` objects are a
subclass of `NDCube` from the Sunpy-affiliated package of the same name.

You can slice an `NDCube` object using either array indices or by using physical
coordinates and the `crop_by_coords` method. Please see the `ndcube` documentation
at <https://docs.sunpy.org/projects/ndcube/en/stable/index.html> for more details.

The `EISCube` subclass extends `ndcube` by including two additional features.
First, an extra `wavelength` attribute has been added which contains a 3D array
with the corrected wavelength values at all locations within the cube. Slicing
an `EISCube` will also appropriately slice the wavelength array. Secondly, there
is now a `total_intensity()` method which returns a new, 2D `NDCube` containing
the sum along the wavelength axis.

A few more notes about `EISCube` objects,
* Currently, the data axes are in the same order as the array stored in the HDF5
  file and used in the previous version of "pyeis". Namely, the order is
  (slit_pixel, raster_step, wavelength), i.e. (solar-y, solar-x, wavelength).
* All of the EIS metadata and header information are stored as seperate keys in
  the `meta` dictionary attribute (e.g. the original FITS header is in `.meta['index']`
  while the updated pointing information is in `.meta['pointing']`).

## read_template
Example call signature,
```
read_template('fe_12_195_119.2c.template.h5')
```

`read_template()` is relatively simple. It takes a single argument giving the filename
of a template file and returns an `EISFitTemplate` class instance containing the
the initial fit parameters. Users may view the parameter values by using either the
`print_parinfo()` method or manually inspecting the `template` attribute. For
convenience, there is also a `central_wave` attribute that contains the mean
wavelength value within the template. This can be useful for loading the correct
spectral window using `read_cube`.

## fit_spectra
Example call signature,
```
fit_spectra(eis_raster, fit_template)
```

This is the main fitting routine. The simplest way to use `fit_spectra()` is to
give it two arguments:
1. _inten_ (EISCube object) - An EISCube class instance (or slice) containing one
   or more intensity profiles to be fit. Wavelength and error values will be
   extracted as needed from the EISCube.
2. _template_ (EISFitTemplate object) - An EISFitTemplate class instance containing
   both the general template information (in the `template` attribute) and a
   `parinfo` attribute with the fit parameter dictionary.

Alternatively, you may provide the function with individual data arrays using a
call signature like this,
```
fit_spectra(inten_arr, template_dict, parinfo=parinfo_dict,
            wave=wavelength_arr, errs=error_arr)
```

The function will loop over the data according to its dimensionality. 3D data is
assumed to be a full EIS raster (or a sub region), 2D data is assumed to be a
single EIS slit, and 1D data is assumed to be a single profile.

`fit_spectra()` returns a `EISFitResult` class instance containing the fit
parameter values. As well as assorted metadata. There are two important class
methods that may be of use `.get_params()` and `.get_fit_profile()`.

The method `.get_params()` extracts parameters values by either component number,
name, or pixel coordinate (or any combination of the three). The arguments are:

1. _component_ (int or None, optional) - Integer number (or list of ints) of the
   functional component(s). If set to None, will return the total combined fit
   profile. Default is None.
2. _param_name_ (str, optional) - String name of the requested parameter. If set
   to None, will not filter based on paramater name. Default is None.
3. _coords_ (list or tupple, optional) - (Y, X) coordinates of the requested
   datapoint. If set to None, will instead return the parameters at all locations.
   Default is None
4. _num_wavelengths_ (int, optional) - Number of wavelength values to compute the
   fit intensity at. These values will be equally spaced and span the entire fit
   window. If set to None, will use the observed wavelength values. Default is None.
5. _casefold_ (bool, optional) - If set to True, will ignore case when extracting
   parameters by name. Default is False.

Examples,
```
c0_params = fit_res.get_fit_profile(component=0)
widths = fit_res.get_fit_profile(param_name='width')
```

The `.get_fit_profile()` method may be used to generate the either the combined fit
intensity profile or the profile of a single component function. The method takes
up to three arguments:
1. _component_ (int or None, optional) - Integer number (or list of ints) of the
   functional component(s). If set to None, will return the total combined fit
   profile. Default is None.
2. _coords_ (list or tupple, optional) - (Y, X) coordinates of the requested
   datapoint. If set to None, will instead return the parameters at all locations.
   Default is None
3. _num_wavelengths_ (int, optional) - Number of wavelength values to compute the
   fit intensity at. These values will be equally spaced and span the entire fit
   window. If set to None, will use the observed wavelength values. Default is None.

`get_fit_profile()` returns two arrays, fit_wave & fit_inten, which contain the
wavelengths and corresponding fit intensity values.

Examples,
```
fit_x, fit_y = fit_res.get_fit_profile(coords=(5,5))
c0_fit_x, c0_fit_y = fit_res.get_fit_profile(component=0, num_wavelengths=100)
```
