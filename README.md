# EISPAC - EIS Python Analysis Code
_version 0.42_

This software provides a set of tools for analyzing Hinode / EIS data within a
python environment. The general approach is as follows:
1. Sets of level 1 HDF5 files are processed from the latest EIS level 0 fits files
   and made available online by the NRL EIS team at <https://eis.nrl.navy.mil/>.
   The HDF5 files come in pairs of "data" and "header" files which contain corrected
   count rates, the calibration curve needed to convert counts into intensity,
   and all of the associated metadata and pointing information.
2. This package provides python classes and functions that can read these hdf5
   files, perform all of the necessary calibration and pointing adjustments, and
   create user-friendly python objects that can be manipulated as needed. Also
   included are functions for fitting the intensity profiles using the same
   template files and underlying methodology that is used in the IDL SolarSoft
   environment.

Please note that this package is under active development. Some of the functionality
and organization of the code is subject to change, including the name of the package
itself. If you have any questions or suggestions for future improvements, please email
the development team

## Installation and Requirements
As long as the scripts are in your current working directory or somewhere on
your Python path, you should be able to run the individual programs and experiment
with the functions without issue. However, if you want to load the package as a
single module or use the example scripts, you will need to install it.

### Installing using pip
1. Unzip and copy the "eispac_develop" folder to convenient location on your
   computer (it does not matter where).
2. Open a terminal and navigate to the folder
3. To install: `python -m pip install .` To upgrade: `python -m pip install --upgrade .`

The package should then be installed to the correct location for your current Python
environment. You can now import the package using `import eispac`. Since the package
is not publically available or registered on PyPi,org, the only way to "update" the
package is by installing a new version over the old using the above method (you do not
need to uninstall the old version first, pip will automatically take care of that).

**Important note**: currently all programs, aside from the example programs, that
depend on other scripts within the package are configured to load local versions
of the required scripts whenever the program is run directly (such as when you open
the program in a Python IDE and hit "run") and will use the installed versions of
the scripts whenever imported. This allows for rapid development and testing of
individual programs without needing to update your installed package.

### Required Packages
pip should take care of all of the package dependencies. If it does not, here is
a list of the required packages (older package versions might still work).
* python >= 3.7
* numpy >= 1.18.1
* scipy >= 1.4.1
* matplotlib >= 3.1
* h5py >= 2.9
* astropy >= 3.1
* sunpy >= 1.0.3
* ndcube >= 1.2.1
* wget
* cURL

## Code Organization

As with the previous "pyEIS" version, there are currently three core directories:

1. **convert**: this is the IDL software used to read the level-1 files and write
   them to HDF5. If you're interested in how an HDF5 file is created, look here.
   Note: in addition in a working IDL environment, you will need some FITS files
   in order to actually run this software. The FITS files for the example dataset
   can be found at <https://github.com/hpwarren/pyEIS-test-data>.

2. **EAG**: initial documentation for the previous version of this code, before
   some general refactoring and the implementation of ndcube class objects. While
   this document is not currently up-to-date, it still contains useful information
   about some of the internal details of the code.

3. **eispac**: main python code directory containing all of the scripts required to
   read level 1 HDF5 files and fit templates and fit spectra using mpfit.

   Notable subdirectories:
   * `../eispac/data/`:  Example data containing a full EIS raster from 2019-04-04
     at 13:15:13.
   * `../eispac/examples/`: Example scripts showing how to load and fit the example data
   * `../eispac/templates/`: fit templates for specific spectral lines. These HDF5
     files are direct conversions of the ".genx" files used by the EIS software in IDL.

  * `mpfit.py`: really mpfit.py by Sergey Koposov, which is based on IDL code by Craig
    Markwardt.

Full, updated documentation of the code will be included in a future version. For now,
please refer to the `QUICK_GUIDE.md` text document and the (partial) doc strings of the
individual scripts. It should also be noted that `mpfit.py` was written by Mark Rivers
and Sergey Kopsov and is direct Python port of the `mpfit.pro` IDL procedure written
by Craig Markwardt. As such, much of the documentation online for the IDL version of
the code is still applicable to the Python version (please see the Python doc
for more information).

## TODO list
Here, in no particular order, is a list of some things that still need work.
* Update and expand documentation
* Add complete unit and integration tests
* Add a subclass of `NDCubeSequence` which can hold multiple spectral windows
* Improve memory efficiency of parallel methods in `fit_spectra()`
* Restructure project to use the Sunpy affiliated package template?
