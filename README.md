# EISPAC - EIS Python Analysis Code

This software provides a set of tools for analyzing Hinode/EIS data within a
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

Please note that this package is under active development. If you have any questions or suggestions
for future improvements, please email the development team

## Installation and Requirements

### Installing using pip

1.  Clone "eispac" to convenient location on your computer (it does not matter where).
```   
   > git clone https://github.com/USNavalResearchLaboratory/eispac.git
```   
2.  Open a terminal and navigate to the directory
3.  To install:
```
	> python -m pip install .
```
4.  To upgrade: 
```
	> python -m pip install --upgrade .
```

The package should then be installed to the correct location for your current Python
environment. You can now import the package using `import eispac`. Since the package is not
publically available or registered on PyPi,org, the only way to "update" the package is by
repeating the process above (you do not need to uninstall the old version first, pip will
automatically take care of that). A conda install script will be added in a future update

### Required Packages
pip should automatically install the package dependencies. If it does not, here is
a list of the required packages (older package versions might still work). Note: if 
you are using conda to manage your Python packages, you may wish to install or update 
the dependencies manually first, before installing eispac using pip.
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

### Getting Started

* `QUICK_GUIDE-cli.md`: A very brief description of some command line tools for downloading and
  fitting the data
  
* `QUICK_GUIDE.md`: A very brief description of the EISPAC objects.

* `users_guide/EISPAC_Users_Guide.pdf`: A more detailed guide to the software.

* `notebooks`: Tutorials using Juypter notebooks.

## Code Organization

There are currently three core directories:

1. **eispac**: main python code directory containing all of the programs required to
   read level 1 HDF5 files and fit templates and fit spectra using mpfit.

   Notable subdirectories:
   * `../eispac/core/`:  Main code directory. All functions here are loading into the
     top-level namespace (i.e. eispac.{function name})
   * `../eispac/data/`:  Example data containing a full EIS raster from 2019-04-04
     at 13:15:13.
   * `../eispac/examples/`: Example scripts showing how to load and fit the example data
   * `../eispac/templates/`: fit templates for specific spectral lines. These HDF5
     files are direct conversions of the ".genx" files used by some IDL users.

2. **scripts**: GUI and command line tools

2. **users_guide**: PDF with instructions and examples for how to use the code. Also
   included the LaTeX source code.

The `QUICK_GUIDE.md` text document also give a very brief overview of some key functions.
It should also be noted that `mpfit.py` was written by Mark Rivers and Sergey Kopsov and 
is direct Python port of the `mpfit.pro` IDL procedure written by Craig Markwardt. As such, 
much of the documentation online for the IDL version of the code is still applicable to the 
Python version (please see the Python doc for more information).

## TODO list
Here, in no particular order, is a list of some things that still need work.
* Expand documentation
* Add more unit and integration tests
* Consider adding a subclass of `NDCubeSequence` which can hold multiple spectral windows
* Consider storing the output fit parameters in another `NDCube`
  (might be one too many subclasses to maintain)
* Restructure project to use the Sunpy affiliated package template?
