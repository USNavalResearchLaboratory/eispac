Downloading EIS Data
====================

.. _sec-download:

eis_catalog GUI or Web Browser
------------------------------

The easiest way to search for and download the processed HDF5 files is to
use the  ``eis_catalog`` GUI tool included with EISPAC. This tool provides
an easy interface to the official EIS as-run database and can search based on
a wide range of technical criteria (see the :ref:`sec-catalog` section for details).
Alternatively, you can browse and download files directly from the online archive
(https://eis.nrl.navy.mil/) or use the `~eispac.download.download_hdf5_data()`
function (assuming you know the exact name of the file you want).

Using SunPy's Fido Interface
----------------------------

As of November of 2022, there is now a handy client for Sunpy's ``Fido`` data
interface! `~sumpy.net.Fido` is a unified tool that can query a variety of
solar and heliophysics data repositories and return a standardized format of
results for downloading. A general guide for using `~sumpy.net.Fido can be found
in the `Finding and Downloading Data <https://docs.sunpy.org/en/latest/guide/acquiring_data/fido.html#fido-guide>`_
section of the SunPy User's Guide.

Here is a simple example for initializing and using the EISPAC Fido interface.
**Please note:** At this time, only time-based searches of EIS data are available
through ``Fido``.

.. code:: python

   >>> from sunpy.net import Fido, attrs as a
   >>> import eispac.net
   >>> from eispac.net.attrs import FileType
   >>> results = Fido.search(a.Time('2020-11-09 00:00:00','2020-11-09 01:00:00'),
   ...                       a.Instrument('EIS'),
   ...                       a.Physobs.intensity,
   ...                       a.Source('Hinode'),
   ...                       a.Provider('NRL'),
   ...                       a.Level('1'))  #doctest: +REMOTE_DATA
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
   2020-11-09 00:10:12.000 2020-11-09 00:10:12.999 ...     1 HDF5 header
   2020-11-09 00:10:12.000 2020-11-09 00:10:12.999 ...     1        FITS
   <BLANKLINE>
   <BLANKLINE>
   >>> results = Fido.search(a.Time('2020-11-09 00:00:00','2020-11-09 01:00:00'),
   ...                       a.Instrument('EIS'),
   ...                       a.Physobs.intensity,
   ...                       a.Source('Hinode'),
   ...                       a.Provider('NRL'),
   ...                       a.Level('1'),
   ...                       FileType('HDF5 header'))  #doctest: +REMOTE_DATA
   >>> results  #doctest: +REMOTE_DATA
   <sunpy.net.fido_factory.UnifiedResponse object at ...>
   Results from 1 Provider:
   <BLANKLINE>
   1 Results from the EISClient:
   Source: https://eis.nrl.navy.mil/
   <BLANKLINE>
          Start Time               End Time        ... Level   FileType
   ----------------------- ----------------------- ... ----- -----------
   2020-11-09 00:10:12.000 2020-11-09 00:10:12.999 ...     1 HDF5 header

.. Attention::
   Some laboratory and institution networks (including VPNs) are known to cause
   issues with using certain ``Fido`` data clients. In such cases, your system
   administrator may be able to provide a solution.
