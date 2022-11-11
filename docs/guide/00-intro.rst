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

# An archive of Level-1 data files. These files are saved in the HDF5 file
  format and come in pairs of header and data files. The archive is updated
  regularly and can be found at https://eis.nrl.navy.mil/

# A set of GUI and command line tools. Can be used to quickly search and
  and download EIS data, view available fit templates, and batch process
  multiple files at once using parallel processing.

# The Python package itself. Provides classes and functions that can read the
  Level-1 HDF5 files, perform all of the necessary calibration and pointing
  adjustments, and create user-friendly Python objects that can be manipulated
  as needed. Also included are functions for fitting the intensity profiles
  with multi-Gaussian functions using template files and a Python port of
  the venerable ``MPFIT`` library [#]_. Installing the Python package also
  automatically installs and registers the GUI and command line tools with
  your Python environment.

Requirements
------------

EISPAC depends on a number of Python packages that are commonly used in
scientific and solar research. Normally, the installation process should
automatically check and install missing dependencies, assuming your
environment is configured appropriately. If it does not, you may wish to
try installing the required packages individually first.

-  python >= 3.8

-  numpy >= 1.18

-  scipy >= 1.4

-  matplotlib >= 3.1

-  h5py >= 2.9

-  astropy >= 4.2.1

-  sunpy >= 4.0

-  ndcube >= 2.0

-  pyqt >= 5.9

-  parfive >= 1.5

-  python-dateutil >= 2.8

Standard Installation
---------------------

EISPAC is now available on PyPI. To install, just use the following command,

::

   >>> python -m pip install eispac

To upgrade the package, please use:

::

   >>> python -m pip install --upgrade eispac

pip should automatically install all package dependencies. If it does not, please
see the list of required packages above. Note: if you are using conda to manage your
Python packages, you may wish to install or update the dependencies manually first,
before installing eispac using pip (this is by no means required, but it can help
simplify updating packages).

Manual Installation
-------------------

1.  Download or clone "eispac" to a convenient location on your computer (it does not matter where).

::

   >>> git clone https://github.com/USNavalResearchLaboratory/eispac.git

2.  Open a terminal and navigate to the directory
3.  To install:

::

   >>> python -m pip install .

4.  To upgrade:

::
   >>> python -m pip install --upgrade .


The package should then be installed to the correct location for your current Python
environment. You can now import the package using `import eispac`.

Now, you should be all set to do some science!

.. rubric:: Citations

.. [#] Markwardt, C. B. 2009, in Astronomical Society of the Pacific Conference
   Series, Vol. 411, Astronomical Data Analysis Software and
   Systems XVIII, ed. D. A. Bohlender, D. Durand, & P. Dowler, 251
