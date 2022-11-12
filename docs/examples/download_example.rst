.. _ex-fido:

Downloading EIS Data with Fido
==============================

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
