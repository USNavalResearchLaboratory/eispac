---
title: 'EISPAC - The EIS Python Analysis Code'
tags:
  - Python
  - solar physics
  - astronomy
  - spectroscopy
authors:
  - name: Micah J. Weberg
    orcid: 0000-0002-4433-4841
    corresponding: true
    equal-contrib: true
    affiliation: "1, 2"
  - name: Harry P. Warren
    orcid: 0000-0001-6102-6851
    equal-contrib: true
    affiliation: 2
  - name: Nicholas Crump
    equal-contrib: false
    affiliation: 2
  - name: Will Barnes
    orcid: 0000-0001-9642-6089
    equal-contrib: false
    affiliation: "3, 4"
affiliations:
  - name: George Mason University, USA
    index: 1
  - name: US Naval Research Laboratory, USA
    index: 2
  - name: Department of Physics, American University, USA
    index: 3
  - name: NASA Goddard Space Flight Center, USA
    index: 4
date: 19 August 2022
bibliography: eispac_paper.bib
---

# Summary
Spectral observations of the Sun - particularly in the extreme ultraviolet light (EUV) range - provide key information concerning the elemental composition and physical parameters of solar plasmas. Since its launch in 2006, the EUV Imaging Spectrometer (EIS) on board the Hinode spacecraft [@Culhane:2007] has provided high quality data used in more than 600 publications. The EIS Python Analysis Code (`EISPAC`) provides convenient and easy to use functions and command line tools for searching, downloading, and analyzing EIS data within the scientific Python ecosystem. Additionally, EISPAC interfaces with other packages focused on solar and heliophysics such as `SunPy` [@sunpy] and `NDCube` [@ndcube]. This compatibility supports an array of exciting multi-spacecraft studies with past and present missions. `EISPAC` is accompanied by a new level-1 archive of the EIS data that combines both calibrated data and ancillary information in sets of HDF5 files for efficient storage and accessibility.

# Statement of Need
To date, EIS has taken well over 430,000 observations. This data has been essential to advancing our understanding of the solar atmosphere. Unfortunately, however, the data processing and analysis routines are written in a commercial programming language that, increasingly, universities and research institutions have decided to replace with open-source alternatives such as Python. This shift has made accessing and analyzing EIS data difficult for students and scientists without the appropriate program licenses. `EISPAC` fills this gap and ensures that EIS data analysis remains easily accessible to everyone.

# Key Features
`EISPAC` has three main components: (1) an archive of minimally preprocessed EIS data in the HDF5 file format, (2) a set of GUI and command line tools for quickly searching, downloading, and analyzing EIS data, and (3) the python package itself, which provides classes and functions for detailed analysis.

The HDF5 files are prepared using the official processing routines developed and validated by the EIS instrument team and distributed in the SolarSoftWare package [@ssw] for the Interactive Data Language (IDL) [@idl]. These routines remove known instrumental effects such as bad pixels, and orbital effects. The data are saved in pairs of HDF5 files, a header file containing all metadata and ancillary information and a data file with the observations in units of photon counts. This approach is a pragmatic one that minimizes the risks of introducing errors in processing the data while also storing the information in a widely accessible and well-supported file format. The HDF5 files can be downloaded either directly online at https://eis.nrl.navy.mil/ or using the tools included in `EISPAC`.

The GUI and command line tools included with EISPAC run in Python and are built using the `PyQt` package [@pyqt]. These tools are designed for both discovering new observations and batch analyzing large numbers of files. The included tools are:

   * _eis_catalog_ - a GUI tool for searching the as-run observation database using a variety of parameters such as wavelength, study ID, and science objectives and then downloading the HDF5 files to your local system
   * _eis_browse_templates_ - GUI tool for listing the observed spectral lines and viewing the available, pre-initialized fitting templates.
   * _eis_fit_files_ - command line tool for fitting spectral lines and computed line intensites, dopper velocities, and line widths.
   * _eis_plot_fit_ - command line tool for producing quick-look plots of the results.

The core `EISPAC` package contains all the necessary functions and classes to read EIS data, run spectral fits, and produce maps and images for further scientific analysis. Wherever possible, `EISPAC` takes advantage of the existing Python packages used in solar and heliophysics. Data are loaded into coordinate-aware `EISCube` array objects that are a subclass of `NDCube`. The EUV line spectra are fit with multi-component Gaussian functions using a Python version of the venerable `MPFIT` package [@Markwardt:2009; @mpfitpy] and further accelerated with the Python multiprocessing module. \autoref{fig:fit_example} shows an example of one such fit. The measured line properties can be saved to `.fits` files, which are a standard image format used in solar and astrophysics, and loaded by the `SunPy` package into a `Map` object which allows for advanced analysis such as co-aligning and comparing observations from multiple spacecraft located throughout the solar system and performing accurate reprojections to compatible coordinate systems.

![Example data cutout (left) and a double-Gaussian fit profile (right) for the Fe XI 188.216 $\AA$ line observed on 4 April 2019. The red X shows the location of maximum intensity. \label{fig:fit_example}](ex_fit.png)

# Acknowledgements
EISPAC is based on work supported by NASA's Hinode Project. Important testing and feedback were given by Art Poland, Deborah Baker, Andy To, David Standby, and Richard Morton.

# References
