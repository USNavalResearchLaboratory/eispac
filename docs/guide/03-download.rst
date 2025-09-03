Downloading EIS Data
====================

.. _sec-download:

eis_catalog GUI or Web Browser
------------------------------

The easiest way to search for and download the processed HDF5 files is to
use the  ``eis_catalog`` GUI tool included with EISPAC. This tool provides
an easy interface to the official EIS as-run database and can search based on
a wide range of criteria (see the :ref:`sec-catalog` section for details).
After searching, you can also view the metadata and wavelength lists for each
observation as well as a context image showing the EIS field-of-view on a SDO/AIA 
or SOHO/EIT image (depnding on which was avialable at the time).

Alternatively, you can browse and download files directly from the online archive
(https://eis.nrl.navy.mil/) or use the `~eispac.db.download_hdf5_data()`
function (assuming you know the exact name of the file you want).

Using SunPy's Fido Interface
----------------------------

As of November of 2022, there is now a handy client for Sunpy's ``Fido`` data
interface! `~sumpy.net.Fido` is a unified tool that can query a variety of
solar and heliophysics data repositories and return a standardized format of
results for downloading. A general guide for using `~sumpy.net.Fido`` can be found
in the `Acquiring Data <https://docs.sunpy.org/en/stable/tutorial/acquiring_data/index.html>`_
section of the SunPy User's Guide.

By default, the client for EIS data only returns entries for the data.h5 files
and fetching a data.h5 file will automatically download the matching head.h5
file. The original level-1 FITS files are also available, should a user wish to 
analyze the data using older IDL routines. Fetching a FITS entry will download 
both the "l1" (data) and "er" (errors) files produced by the "EIS_PREP" 
routine in SolarSoft/IDL (Note: these FITS files contain exactly the same data 
as the HDF5 files, but are missing some of the metadata which are computed by 
SSW routines).

European users may get better download speeds from the data mirror hosted by
the Mullard Space Science Laboratory (MSSL) in the UK. To use it, set 
``a.Provider(MSSL)``

Here is a simple example for initializing and using the EISPAC Fido interface.
**Important note:** At this time, only time-based searches of EIS data are available
through ``Fido``. For more complex searches, please use the ``eis_catalog`` GUI.

.. code:: python

   >>> from sunpy.net import Fido, attrs as a
   >>> import eispac.net
   >>> results = Fido.search(a.Time('2020-11-09 00:00:00','2020-11-09 01:00:00'),
   ...                       a.Instrument('EIS'),
   ...                       a.Source('Hinode'),
   ...                       a.Provider('NRL'))  #doctest: +REMOTE_DATA
   >>> results  #doctest: +REMOTE_DATA
   <sunpy.net.fido_factory.UnifiedResponse object at ...>
   Results from 1 Provider:
   <BLANKLINE>
   3 Results from the EISClient:
   Source: https://eis.nrl.navy.mil/
   <BLANKLINE>
          Start Time               End Time        ... Level   FileType
   ----------------------- ----------------------- ... ----- -----------
   2020-11-09 00:10:12.000 2020-11-09 00:10:12.999 ...     1   HDF5 data
   <BLANKLINE>
   <BLANKLINE>
   >>> results = Fido.search(a.Time('2020-11-09 00:00:00','2020-11-09 01:00:00'),
   ...                       a.Instrument('EIS'),
   ...                       a.Source('Hinode'),
   ...                       a.Provider('NRL'),
   ...                       a.eispac.FileType('Any'))  #doctest: +REMOTE_DATA
   >>> results  #doctest: +REMOTE_DATA
   <sunpy.net.fido_factory.UnifiedResponse object at ...>
   Results from 1 Provider:
   <BLANKLINE>
   1 Results from the EISClient:
   Source: https://eis.nrl.navy.mil/
   <BLANKLINE>
          Start Time               End Time        ... Level   FileType
   ----------------------- ----------------------- ... ----- -----------
   2020-11-09 00:10:12.000 2020-11-09 00:10:12.999 ...     1   HDF5 data
   2020-11-09 00:10:12.000 2020-11-09 00:10:12.999 ...     1 HDF5 header
   2020-11-09 00:10:12.000 2020-11-09 00:10:12.999 ...     1        FITS

.. Attention::
   Some laboratory and institution networks (including VPNs) are known to cause
   issues with using certain ``Fido`` data clients. In such cases, your system
   administrator may be able to provide a solution.
