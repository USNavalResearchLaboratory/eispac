.. _sec-download:

Command Line Scripts
====================

There are two main components to the EISPAC software: (1) a set of
command line scripts and (2) the Python package itself. In this chapter
we will give a very brief overview of how to use the scripts to download data,
explore line templates, and batch process multiple observations at once.
The next two chapters will cover how to read and manipulate the level-1 HDF
files directly and fit spectra with Gaussian functions.

Overview of scripts
-------------------

The command line scripts should be automatically installed and
registered with the OS as part of installing EISPAC. These scripts are
designed to help users quickly browse, download, and fit Gaussian
functions to the data. To use a script, simply enter its name in the
command line from any directory in which you have read and write
privileges.

.. tip::
   Some scripts will default to saving files to your current working
   directory, therefore we recommend running the scripts from the directory
   in which you intend to do most of your analysis.

There are currently four command line scripts available,

-  ``eis_catalog`` - GUI tool from searching the as-run EIS data catalog and
   downloading the HDF5 files your computer. Can also generate a text list of
   files to download.

-  ``eis_browse_templates`` - GUI tool for browsing the fit templates
   corresponding to each spectral window in a given observation set and
   copying the template files from EISPAC to your current working directory
   (fit templates are explained more in the next chapter)

-  ``eis_download_files`` - Command line tool for downloading a the level-1
   HDF5 files associated with one or more level-0 EIS fits files. Can also
   download an entire list of files using the text output of ``eis_catalog``.
   Example usage,

   ::

      >>> eis_downdload_files eis_l0_20190404_131513.fits

-  ``eis_fit_files`` - Command line tool for fitting all of the HDF5 files
   in a given directory with each fit template found in another directory.
   Example usage,

   ::

      >>> eis_fit_files ./eis_study/ ./eis_study/templates/
