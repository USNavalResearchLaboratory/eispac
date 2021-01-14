__all__ = ['fit_spectra_astropy']

import copy
from datetime import datetime
import numpy as np
from astropy.modeling import fitting
import eispac.core.fitting_functions as fit_fns
from eispac.core.eiscube import EISCube
from eispac.core.read_template import EISFitTemplate
from eispac.core.eisfitresult import EISFitResult
from eispac.core.scale_guess import scale_guess
from eispac.core.mpfit import mpfit
from eispac.core.generate_astropy_model import generate_astropy_model

def fit_spectra_astropy(inten, template, parinfo=None, wave=None, errs=None,
                        min_points=10):
    """Fit one or more EIS line spectra using mpfit.

    Parameters
    ----------
    inten : EISCube object or array_like
        One or more intensity profiles to be fit. The code will loop over the data
        according to its dimensionality. 3D data is assumed to be a full EIS raster
        (or a sub region), 2D data is assumed to be a single EIS slit, and 1D data
        is assumed to be a single profile.
    template : EISFitTemplate object or dict
        Either an EISFitTemplate or the just a 'template' dictionary.
    parinfo : dict, optional
        Dictionary of fit parameters formatted for use with mpfit. Required if
        the 'template' parameter is given just a dict and ignored otherwise.
    wave : array_like, optional
        Associated wavelength values for the spectra. Required if 'inten' is
        given an array and ignored otherwise.
    errs : array_like, optional
        Intensity error values for the spectra. Required if 'inten' is given an
        array and ignored otherwise.
    min_points : int, optional
        Minimum number of good quality data points (i.e. non-zero values & errs)
        to be used in each fit. Spectra with fewer data points will be skipped.
        Default is 10.

    Returns
    -------
    fit_res :EISFitResult class instance
        An EISFitResult object containing the output fit paramaters.
    """
    # If given an EISCube, extract data arrays and zero out masked values.
    # Otherwise, make local copies and just zero out negative values
    if isinstance(inten, EISCube):
        wave_cube = inten.wavelength.copy()
        errs_cube = inten.uncertainty.array.copy()
        inten_cube = inten.data.copy()
        metadata = inten.meta
        loc_masked = np.where(inten.mask == True)
        inten_cube[loc_masked] = 0
        errs_cube[loc_masked] = 0
    else:
        wave_cube = wave.copy()
        errs_cube = errs.copy()
        inten_cube = inten.copy()
        metadata = {'filename_data':'unknown', 'index':{}, 'pointing':{}}
        loc_neg = np.where(inten_cube < 0)
        inten_cube[loc_neg] = 0
        errs_cube[loc_neg] = 0

    # Check contents of template and parinfo
    if parinfo is None and isinstance(template, EISFitTemplate):
        parinfo_copy = copy.deepcopy(template.parinfo)
        template_copy = copy.deepcopy(template.template)
    elif parinfo is None:
        raise TypeError('Please input either a "parinfo" dictionary or a full'
                        ' EISFitTemplate object.')
    else:
        parinfo_copy = copy.deepcopy(parinfo)
        template_copy = copy.deepcopy(template)

    # Check the dimensions of input data.
    # If the the arrays are not 3D, add shallow dimensions of size 1
    num_dims = inten_cube.ndim
    dims_size = inten_cube.shape
    if num_dims == 1:
        n_pxls = 1
        n_steps = 1
        wave_cube = wave_cube[np.newaxis, np.newaxis, :]
        inten_cube = inten_cube[np.newaxis, np.newaxis, :]
        errs_cube = errs_cube[np.newaxis, np.newaxis, :]
    elif num_dims == 2:
        n_pxls = dims_size[0]
        n_steps = 1
        wave_cube = wave_cube[:, np.newaxis, :]
        inten_cube = inten_cube[:, np.newaxis, :]
        errs_cube = errs_cube[:, np.newaxis, :]
    elif num_dims == 3:
        n_pxls = dims_size[0]
        n_steps = dims_size[1]
        wave_cube = wave_cube
        inten_cube = inten_cube
        errs_cube = errs_cube

    # Initalize output object which will contain the fit results
    fit_res = EISFitResult(wave_cube, template_copy, parinfo_copy, func_name='multigaussian')
    fit_res.meta = metadata
    fit_res.fit_module = 'astropy'
    fit_res.fit_method = 'LevMarLSQ'

    # get the indices of the function parameters
    n_gauss = fit_res.n_gauss
    n_poly = fit_res.n_poly
    # peaks = np.arange(n_gauss)*3
    # cents = np.arange(n_gauss)*3+1
    # wdths = np.arange(n_gauss)*3+2
    # backs = n_gauss*3

    # timer
    t1 = datetime.now()

    # oldguess = template_copy['fit']
    fitter = fitting.LevMarLSQFitter()

    tmplt_model = generate_astropy_model(template.filename_temp)
    tmplt_init_param_vals = copy.deepcopy(tmplt_model.parameters)

    # Get indices of parameter that are NOT tied
    tied_mask = np.zeros(tmplt_model.parameters.size, dtype=bool)
    for p in range(tmplt_model.parameters.size):
        if tmplt_model.__dict__[tmplt_model.param_names[p]].tied != False:
            tied_mask[p] = True
    loc_untied = np.where(tied_mask == False)

    # loop over pixel positions and slit steps in the entire raster
    for ii in range(n_pxls):
        for jj in range(n_steps):

            # Extract a single profile from the raster
            wave_ij = wave_cube[ii,jj,::]
            inten_ij = inten_cube[ii,jj,::]
            errs_ij = errs_cube[ii,jj,::]

            # Only use good data
            loc_good = np.where(errs_ij > 0)
            if len(loc_good[0]) < min_points:
                fit_res.fit['status'][ii,jj] = -1
                continue
            wave_ij = wave_ij[loc_good]
            inten_ij = inten_ij[loc_good]
            errs_ij = errs_ij[loc_good]

            # Do the fit using Astropy
            fit_model = fitter(tmplt_model, wave_ij, inten_ij,
                               weights=1.0/errs_ij, maxiter=2000)

            # check convergence status
            if fitter.fit_info['ierr'] in [1, 2, 3, 4]:
                # update initial parameter values
                setattr(tmplt_model, 'parameters', fit_model.parameters)
                # assemble fit structure
                # TO-DO: figure out a more compact way to do this that is still
                #        easy to read and debug
                fit_dof = wave_ij.size - loc_untied[0].size
                chi_sq = sum(((inten_ij - fit_model(wave_ij))/errs_ij)**2)
                if fitter.fit_info['param_cov'] is not None:
                    untied_perror = np.sqrt(np.diag(fitter.fit_info['param_cov']))
                else:
                    untied_perror = np.zeros(loc_untied[0].size)
                fit_res.fit['status'][ii,jj] = fitter.fit_info['ierr']
                fit_res.fit['chi2'][ii,jj] = chi_sq/fit_dof
                fit_res.fit['wavelength'][ii,jj,::] = wave_cube[ii,jj,::]
                # fit_res.fit['peak'][ii,jj,::] = out.params[peaks]
                # fit_res.fit['err_peak'][ii,jj,::] = out.perror[peaks]
                # fit_res.fit['centroid'][ii,jj,::] = out.params[cents]
                # fit_res.fit['err_centroid'][ii,jj,::] = out.perror[cents]
                # fit_res.fit['width'][ii,jj,::] = out.params[wdths]
                # fit_res.fit['err_width'][ii,jj,::] = out.perror[wdths]
                # fit_res.fit['background'][ii,jj,::] = out.params[backs]
                # fit_res.fit['err_background'][ii,jj,::] = out.perror[backs]
                fit_res.fit['params'][ii,jj,::] = fit_model.parameters
                fit_res.fit['perror'][ii,jj,loc_untied] = untied_perror
                # fit_res.fit['int'][ii,jj,::] = int
                # fit_res.fit['err_int'][ii,jj,::] = err_int
            else:
                print(' ! fit did not converge!')
                fit_res.fit['status'][ii,jj] = fitter.fit_info['ierr']

    # print status
    print(' + fit completed!')

    # timer
    t2 = datetime.now()
    print(' + fit runtime : {}'.format(t2-t1))

    return fit_res



if __name__ == '__main__':

    import pathlib
    import matplotlib.pyplot as plt
    import astropy.units as u
    # from read_cube import read_cube
    # from read_template import read_template

    # input data and template files
    file_data = './data/eis_20190404_131513.data.h5'
    file_data = str(pathlib.Path(file_data).resolve())
    file_template = './templates/eis_template_dir/fe_12_195_119.2c.template.h5'
    file_template = str(pathlib.Path(file_template).resolve())

    # read fit template
    Fe_XII_195_119 = read_template(file_template)

    # read spectra window
    raster = read_cube(file_data, Fe_XII_195_119.central_wave)

    # fit profile
    wave_coords = raster.axis_world_coords('em.wl')
    lower_corner = (300*u.arcsec, 50*u.arcsec, wave_coords[0])
    upper_corner = (400*u.arcsec, 150*u.arcsec, wave_coords[-1])
    # sub_raster = raster.crop_by_coords(lower_corner, upper_corner=upper_corner)
    # sub_raster = raster[0:10, 0:10, :]
    sub_raster = raster[:, 40:50, :]

    fit_res = fit_spectra_astropy(sub_raster, Fe_XII_195_119)

    # Quick plot raster
    plot_aspect_ratio = raster.meta['pointing']['y_scale']/raster.meta['pointing']['x_scale']
    sub_raster[:,:,12].plot(aspect=plot_aspect_ratio)

    # Plot example fit
    ex_pxl_coords = [5, 5]
    # ex_pxl_coords = [4, 8]
    fit_x, fit_y = fit_res.get_fit_profile(coords=ex_pxl_coords, num_wavelengths=100)
    c0_fit_x, c0_fit_y = fit_res.get_fit_profile(component=0, coords=ex_pxl_coords,
                                                 num_wavelengths=100)
    c1_fit_x, c1_fit_y = fit_res.get_fit_profile(component=1, coords=ex_pxl_coords,
                                                 num_wavelengths=100)
    c2_fit_x, c2_fit_y = fit_res.get_fit_profile(component=2, coords=ex_pxl_coords,
                                                 num_wavelengths=100)
    sub_data = sub_raster.data[ex_pxl_coords[0], ex_pxl_coords[1], :]
    sub_wave = sub_raster.wavelength[ex_pxl_coords[0], ex_pxl_coords[1], :]
    sub_err = sub_raster.uncertainty.array[ex_pxl_coords[0], ex_pxl_coords[1], :]

    fig = plt.figure()
    profile_subplt = fig.add_subplot(111)
    profile_subplt.errorbar(sub_wave, sub_data, yerr=sub_err,
                            ls='', marker='o', color='k')
    profile_subplt.plot(fit_x, fit_y, color='b')
    profile_subplt.plot(c0_fit_x, c0_fit_y, color='r')
    profile_subplt.plot(c1_fit_x, c1_fit_y, color='r', ls='--')
    profile_subplt.plot(c2_fit_x, c2_fit_y, color='g')
    profile_subplt.set_xlabel('Wavelength [$\AA$]')
    profile_subplt.set_ylabel('Intensity ['+sub_raster.unit.to_string()+']')
    plt.show()
