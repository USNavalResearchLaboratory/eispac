.. _sec-prep:

Level-1 HDF5 Files
==================

This chapter describes in more detail how the EIS level-1 HDF5 files
were processed and saved. These HDF5 files can be downloaded from
https://eis.nrl.navy.mil/ or by using the search & download
functionality of EISPAC (see the section on :ref:`sec-download`).

Prepping the data in IDL
------------------------

The level-0 fits files were prepped using the IDL routine ``eis_prep``
available via SolarSoft [#]_ with the following options:

::

     default = 1
     save = 1
     quiet = 1
     retain = 1
     photons = 1
     refill = 0

There are 400,000+ EIS level-0 files at present, but on a multi-core
machine using the IDL bridge all of the files can be prepped in under 24
hours. We have prepped all of the available EIS files and saved them to
standard fits files in the usual way. Some important points:

* *units* - As mentioned previously, the units for the output in these level-1
  files is “photon events” or “counts.” This means that the statistical
  uncertainty can usually be estimated as :math:`\sqrt{N}`. When EISPAC
  loads the HDF5 files, it also calculates an estimate for the read
  noise contribution. It should be noted, however, that the read noise
  becomes significant only at very low flux levels (1–2 counts).

* *retain* - Note that the retain keyword preserves negative values. One of
  the jobs of ``eis_prep`` is to remove the pedestal from the CCD readout
  and any time-dependent dark current. Since the spectral windows are
  generally narrow, the estimate of the background can be too high and
  the subtracted intensities of the continuum can be negative. This
  will be dealt with during the fitting.

* *refill* - The warm pixel problem complicates the fitting of EIS line
  profiles. As discussed in the EIS software note #13 (found in SSW or on the
  eiswiki), interpolating the values of missing pixels appears to best
  reproduce the original data. This option is left off during
  ``eis_prep`` so that the level-1 fits file preserves the information
  on the missing pixels. As discussed below, the interpolation (via the
  refill option) is done during the read and this data is ultimately
  written to the HDF5 file. A mask indicating which pixels have been
  interpolated will be added to the HDF5 files in a future revision.

Here is an IDL code snippet related to reading the data by looping over
the spectral windows.

.. code:: idl

   for iwin=0, nwin-1 do begin
     d = eis_getwindata(eis_level1_filename, iwin, /refill, /quiet)
     eis_level1_data[iwin] = ptr_new(d)
   endfor

Writing the HDF5 files
----------------------

Each processed level-1 fits file was bundled up with the associated
calibration and metadata and saved as a pair of two HDF5 files:

* **eis_YYYYMMDD_HHMMSS.data.h5** - Contains only the corrected photon counts
  within each spectral window. This is, by far, the larger of the two HDF5
  files. However, they should not need to be updated or downloaded very often.

* **eis_YYYYMMDD_HHMMSS.head.h5** - Contains the original fit file index,
  calibration curves for each spectral window (used to convert counts into
  intensity values), and the corrected pointing information.

Internal Structure
------------------

Users will rarely, if ever, need to access the information inside the
HDF5 files directly. EISPAC contains all of the functions needed to read
the data and apply the calibration and pointing corrections.
Nevertheless, the contents and data structure of the HDF5 files are
summarized below.

**.data.h5**

- **level1** (group)

   - **intensity_units** (dataset) - String with the level1 data units.
     This will usually be "counts"

   - **win##** (dataset) - array of floating point values with the photon
     counts in a given spectral window (e.g. ``win00, win01, ... win24``).
     Note well, each set of EIS observations may have a different
     number of spectral windows, up to a maximum of 24 windows. Window
     numbers are numbered sequentially from 00; a given wavelength
     range may be assigned a different window number in each EIS study.

**.head.h5**

- *Details to be added soon*

Additionally, the contents of the HDF5 files can be displayed using
``h5dump`` command line tool, which is provided along with the Anaconda
Python distribution platform or can be installed on its own. Example
usage,

::

   > h5dump -n eis_20190404_131513.data.h5
   FILE_CONTENTS {
   group      /
   group      /level1
   dataset    /level1/intensity_units
   dataset    /level1/win00
   dataset    /level1/win01
   dataset    /level1/win02
   dataset    /level1/win03
   dataset    /level1/win04
   dataset    /level1/win05
   . . .

The actual data associated with each variable can be printed out using
the ``-d`` option. For example,

::

   > h5dump -d exposure_times/duration eis_20190404_131513.head.h5
   HDF5 "eis_20190404_131513.head.h5" {
   DATASET "exposure_times/duration" {
      DATATYPE  H5T_IEEE_F32LE
      DATASPACE  SIMPLE { ( 87 ) / ( 87 ) }
      DATA {
      (0): 40.0005, 40.0002, 40.0004, 40.0004, 39.9994, 40.0002, 39.9995, 40,
      (8): 40.0007, 39.9999, 40.0005, 40.0004, 39.9997, 40.0002, 39.9994,
      . . .
   }

.. rubric:: Citations

.. [#] Freeland, S. L., & Handy, B. N. 1998, *Sol.* *Phys.*, **182**, 497
