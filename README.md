# EISPAC - EIS Python Analysis Code

[![eispac CI status](https://github.com/USNavalResearchLaboratory/eispac/workflows/Tests/badge.svg)
](https://github.com/USNavalResearchLaboratory/eispac/actions/workflows/tests.yml)
[![Documentation Status](https://readthedocs.org/projects/eispac/badge/?version=latest)](https://eispac.readthedocs.io/en/latest/?badge=latest)


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

### Using PIP

EISPAC is now available on PyPI. To install, just use the following command,
```
	> python -m pip install eispac
```

To upgrade the package, please use:
```
	> python -m pip install --upgrade eispac
```

pip should automatically install all package dependencies. If it does not, please
see the list of required packages below. Note: if you are using conda to manage your
Python packages, you may wish to install or update the dependencies manually first,
before installing eispac using pip (this is by no means required, but it can help
simplify updating packages).

### Manual Install

1.  Download or clone "eispac" to a convenient location on your computer (it does not matter where).
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
environment. You can now import the package using `import eispac`.

### Required Packages

* python >= 3.8
* numpy >= 1.18
* scipy >= 1.4
* matplotlib >= 3.1
* h5py >= 2.9
* astropy >= 4.2.11
* sunpy >= 4.0
* ndcube >= 2.0.0
* parfive >= 1.5
* python-dateutil>=2.8

### Getting Started

* **Online user's guide**: <https://eispac.readthedocs.io/en/latest/index.html>:

* `QUICK_GUIDE-cli.md`: A very brief description of some command line tools for searching,
  downloading, and fitting the EIS observations

* `QUICK_GUIDE.md`: A very brief description of EISPAC functions and objects.

* `examples`: Tutorials using Juypter notebooks. In particular, `eispac_tutorial.ipynb` contains
  complete overview and introduction to using EISPAC

## Code Organization

There are currently three core directories:

1. **eispac**: main python code directory containing all of the programs required to
   read level 1 HDF5 files and fit templates and fit spectra using mpfit.

   Notable subdirectories:
   * `../eispac/core/`:  Main code directory. All functions here are loaded into the
     top-level namespace (i.e. eispac.{function name})
   * `../eispac/data/`: Contains fitting templates for specific spectral lines. These HDF5
     files are direct conversions of the ".genx" files used by some IDL users. Also included
     is an example EIS raster from 2021-03-06 at 06:44:44.

2. **scripts**: GUI and command line tools

2. **docs**: Source reStructuredText files used to build the online documentation

The `QUICK_GUIDE.md` text document also give a very brief overview of some key functions.
It should also be noted that `mpfit.py` was written by Mark Rivers and Sergey Kopsov and
is direct Python port of the `mpfit.pro` IDL procedure written by Craig Markwardt. As such,
much of the documentation online for the IDL version of the code is still applicable to the
Python version (please see the Python doc for more information).

## TODO list
Here, in no particular order, is a list of some things that may be added in future releases.
* Expanded documentation
* More unit and integration tests
* More detailed logging (with option to send all log information to a file)
* Scripts for quickly viewing data and spectra fits
* Scripts and routines for creating new fit templates
* Consider adding a subclass of `NDCubeSequence` which can hold multiple spectral windows
* Consider storing the output fit parameters in another `NDCube`
* Restructure project to use the Sunpy affiliated package template?
