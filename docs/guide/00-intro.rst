Introduction
============

The EIS Python Analysis Code (EISPAC) is a software package designed for easy
and accurate analysis of spectroscopic data from the EIS instrument on board
the Hinode spacecraft.

What Can EISPAC Do?
-------------------

-  Search and download EIS observations

-  Explore data contents

-  Fit spectral line using multigaussian fit templates

-  Compute measurements (line intensities, velocities, and widths)

-  Generate coordinate-aware, Sunpy Maps of the fit measurements


Subcomponents of EISPAC
-----------------------

There are three main subcomponents of EISPAC,

1. **An archive of Level-1 data files**

   These files are saved in the HDF5 file format and come in pairs of header
   and data files. The archive is updated regularly and can be found at
   https://eis.nrl.navy.mil/

2. **A set of GUI and command line tools**

   Can be used to quickly search and and download EIS data, view available fit
   templates, and batch process multiple files at once using parallel processing.

3. **The Python package itself**

   Provides classes and functions that can read the Level-1 HDF5 files, perform
   all of the necessary calibration and pointing adjustments, and create
   user-friendly Python objects that can be manipulated as needed. Also included
   are functions for fitting the intensity profiles with multi-Gaussian functions
   using template files and a Python port of the venerable ``MPFIT`` library [#]_.
   Installing the Python package also automatically installs and registers the
   GUI and command line tools with your Python environment.


Simple Example
--------------
Reading and fitting EIS data requires minimal time and effort. Below is a basic,
but complete, example:

.. code:: python

   import eispac

   if __name__ == '__main__':
       # input data and template files
       data_filepath = './eis_20190404_131513.data.h5'
       template_filepath = './fe_12_195_119.2c.template.h5'

       # read fit template
       tmplt = eispac.read_template(template_filepath)

       # Read spectral window into an EISCube
       data_cube = eispac.read_cube(data_filepath, tmplt.central_wave)

       # Fit the data using parallel processing
       fit_res = eispac.fit_spectra(data_cube, tmplt, ncpu='max')

.. rubric:: Citations

.. [#] Markwardt, C. B. 2009, in Astronomical Society of the Pacific Conference
   Series, Vol. 411, Astronomical Data Analysis Software and
   Systems XVIII, ed. D. A. Bohlender, D. Durand, & P. Dowler, 251
