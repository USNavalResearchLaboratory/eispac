Introduction
============

The EUV Imaging Spectrometer — EIS — was designed to study the solar atmosphere
and answer fundamental questions on the heating of the solar corona, the origin
of the solar wind, and the release of energy in solar flares [#]_. EIS observes
two wavelength ranges in the extreme ultraviolet, 171–-212 Å and
245–-291 Å with a spectral resolution of about 22 mÅ and a plate scale
of 1″ per pixel. Solar images can be made by stepping the slit over a
region of the Sun and taking an exposure at each position. A detailed
description of EIS is given in the instrument [#]_.

This document describes the basic elements of EIS data analysis using
new HDF5 level-1 files and the EIS Python Analysis Code (EISPAC)
package. At the beginning of the Hinode mission the strategy was to
release unprocessed level-0 FITS files and software routines written in
IDL for processing these files into a format that could be used for data
analysis. Additionally, all of the routines for computing ancillary
information, such as the offsets of the detectors or the magnitude of
the instrumental broadening, were all written in IDL. Unfortunately, IDL
is an expensive, proprietary language, little used outside of solar
physics. Python, in contrast, is a free, open source language that has
grown dramatically in popularity since the launch of Hinode, making it
an obvious choice for future software development.

EIS Level-1 HDF5 Files
----------------------

To accelerate the transition to Python we have created a new level-1
product that contains both the processed level-1 data and the ancillary
information needed for data analysis. The alternative approach, to port
all of the existing IDL software to Python, would be time consuming and
create confusion about which routines are being actively supported
during the transition. Distributing level-1 files removes this problem,
but does make the user dependent on the team for reformatting all of the
files as bugs are discovered. Since the mission has been going on for
some time now, the number of bugs is likely to be small.

There are several other design decisions that merit some explanation

-  The data and header information are stored in separate files. Since
   the data is large and unlikely to change, the time-consuming download
   of these files should only need to be done once. The header file is
   very small and can be updated easily.

-  HDF5 is used to store the data. This is a very widely used,
   high-performance file format that is well supported by both IDL and
   Python. The most attractive feature for this application is that data
   is stored in a self-documenting, directory-like tree structure
   instead of binary table extensions.

-  The data is processed from raw “data numbers” to “photon events” or
   “counts”. The default behavior of ``eis_prep`` is to convert to
   calibrated units. With the HDF5 files conversion to absolute units is
   done using a calibration curve in the header file, and several
   different calibration curves can be considered.

The processed level-1 HDF5 files can be downloaded either directly from
the NRL Hinode/EIS website at https://eis.nrl.navy.mil/ or using the
tools included in EISPAC. The chapter on :ref:`sec-prep` describes how
the files were processed.

EIS Python Analysis Code (EISPAC)
---------------------------------

EISPAC provides Python classes and functions that can read the new HDF5
files, perform all of the necessary calibration and pointing
adjustments, and create user-friendly Python objects that can be
manipulated as needed. Also included are functions for fitting the
intensity profiles with multi- Gaussian functions using template files
and a Python port of the venerable ``MPFIT`` library [#]_. For convenience,
command line and GUI tools are provided to help users quickly browse and
download data, copy template files, and fit multiple files at once using
parallel processing.

Requirements
~~~~~~~~~~~~

EISPAC depends on a number of Python packages that are commonly used in
scientific and solar research. Normally, the installation process should
automatically check and install missing dependencies, assuming your
environment is configured appropriately. If it does not, you may wish to
try installing the required packages individually first.

-  python >= 3.7

-  numpy >= 1.18.1

-  scipy >= 1.4.1

-  matplotlib >= 3.1

-  h5py >= 2.9

-  astropy >= 3.1

-  sunpy >= 1.0.3

-  ndcube >= 1.2.1

-  pyqt >= 5.9

Additionally, some of the command line tools in EISPAC depend on two
non-Python software packages - ``wget`` and ``cURL``. Both of these
packages should come preinstalled on most modern operating systems. If
your system does not, please refer to the respective project websites
and/or contact your system administrator.

Installation
~~~~~~~~~~~~

The current release of EISPAC is not yet available on the usual Python
web repositories (PyPi or Conda). As such, the installation process is a
little bit more involved than other packages.

#. Download the entire ``eispac`` repository and extract it to
   a convenient directory on your computer (it does not matter where).

#. Open a terminal and navigate to the directory chosen above

#. Run the install script using your preferred package manager,

   #. **PIP**: type ``python -m pip install .``

   #. **conda**: DETAILS FORTHCOMING

If you later wish to update EISPAC you will need to repeat steps 1 & 2
above and then issue the command ``python -m pip install --upgrade .``

You should be all ready to go now!

.. rubric:: Footnotes and Citations

.. [#] EIS is part of the Hinode mission and was sponsored by the Japan
   Aerospace Exploration Agency (JAXA), the United Kingdom Space Agency (UKSA),
   and National Aeronautics and Space Administration (NASA) with contributions
   from ESA and Norway. Hinode was launched on September 22, 2006 at 21:36 UTC
   from the Uchinoura Space Center in Japan and continues to operate.

.. [#] Culhane, J. L., Harra, L. K., James, A. M., et al. 2007, *Sol.* *Phys.*, **243**, 19

.. [#] Markwardt, C. B. 2009, in Astronomical Society of the Pacific Conference
   Series, Vol. 411, Astronomical Data Analysis Software and
   Systems XVIII, ed. D. A. Bohlender, D. Durand, & P. Dowler, 251
