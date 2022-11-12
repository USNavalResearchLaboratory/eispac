__all__ = ['fit_spectra']

import sys
import copy
import pathlib
import inspect
from datetime import datetime
import multiprocessing as mp
import numpy as np
import eispac.core.fitting_functions as fit_fns
from eispac.core.eiscube import EISCube
from eispac.core.read_cube import read_cube
from eispac.core.read_template import EISFitTemplate
from eispac.core.eisfitresult import EISFitResult
from eispac.core.eisfitresult import create_fit_dict
from eispac.core.scale_guess import scale_guess
from eispac.core.mpfit import mpfit
from eispac.instr import calc_velocity

# i'm sure that there is a better way to do this!
cntr = mp.Value("i", 0) # a counter
nexp = mp.Value("i", 0) # total exposures

def check_for_name_guard(debug=False):
    """Check to see if the top level script has a "name guard" that will
    protect against multiprocessing spawning infinite child processes.
    """
    name_guard = False
    # Walk up the call stack and find the first program or script with
    # __name__ == "__main__". This should be the top-level program.
    if debug:
        print('')
        print('Checking call stack for a name guard...')
    call_stack = inspect.stack()
    for s in range(len(call_stack)):
        frame = call_stack[s][0]
        frame_filename = call_stack[s][1]
        if debug:
            print('Stack frame index:', s)
            print('   filename = ', frame_filename)
        if '__name__' in frame.f_locals:
            if debug:
                print('   __name__ = ', frame.f_locals['__name__'])
                print('   pathlib.Path.is_file() = ',
                      str(pathlib.Path(frame_filename).is_file()))
            if frame.f_locals['__name__'] == '__main__':
                if not frame_filename.endswith('.py'):
                    # Probably running in a Python terminal or interactive shell
                    name_guard = False
                else:
                    if pathlib.Path(frame_filename).is_file() == False:
                        # Probably running in an entry point executable
                        name_guard = True
                        break
                    # Examine the source code of __main__ script
                    with open(frame_filename, 'r') as f_file:
                        for line in f_file.readlines():
                            sl = line.replace(' ','') # Remove whitespace
                            sl = sl.replace('"', "'") # Convert " to '
                            if sl.startswith("if__name__=='__main__':"):
                                name_guard = True
                break

    return name_guard

def fit_with_mpfit(wave_cube, inten_cube, errs_cube, template, parinfo,
                   min_points=7,  chunk=1, data_units='unknown'):
    """Helper function for fit_spectra(). Fits one or more intensity spectra
    using the mpfit module.
    """

    with cntr.get_lock():
        cntr.value += 1
    #print(f' + working on exposure {cntr.value:03d}') #, end='\r')
    print(f' + working on exposure {chunk:03d}', end='\r')

    inten_size = inten_cube.shape
    n_pxls = inten_size[0]
    n_steps = inten_size[1]
    n_wave = inten_size[2]
    line_ids = template['line_ids']
    if isinstance(line_ids, bytes):
        line_ids = [line_ids.decode('utf-8')]

    # Extract fit information from template
    n_gauss = template['n_gauss']
    n_poly = template['n_poly']
    min_wave = template['data_x'][0]
    max_wave = template['data_x'][-1]
    oldguess = template['fit']

    # Create fit dictionary, mask array, and parameter indices
    fit_dict = create_fit_dict(n_pxls, n_steps, n_wave, n_gauss, n_poly,
                               data_units=data_units)
    fit_dict['line_ids'] = line_ids
    fit_dict['wave_range'][0] = min_wave
    fit_dict['wave_range'][1] = max_wave
    mask_ij = np.ones(n_wave)
    loc_peaks = np.arange(n_gauss)*3
    loc_cen = np.arange(n_gauss)*3+1
    loc_wid = np.arange(n_gauss)*3+2
    loc_backs = np.arange(n_poly)+n_gauss*3

    # loop over pixel positions and slit steps in the entire raster
    for ii in range(n_pxls):
        for jj in range(n_steps):

            # Extract a single profile from the raster
            wave_ij = wave_cube[ii,jj,::]
            inten_ij = inten_cube[ii,jj,::]
            errs_ij = errs_cube[ii,jj,::]

            # Cut data to only "good" nonzero errs within the wavelength range
            mask_ij[0:-1] = 1 # reset the mask (remember: 1 = True = masked val)
            loc_good = np.where((errs_ij > 0) & (wave_ij >= min_wave)
                                & (wave_ij <= max_wave))
            num_good_data = len(loc_good[0])
            if num_good_data > 0:
                mask_ij[loc_good] = 0 # unmask good data

            fit_dict['mask'][ii,jj,:] = mask_ij[:] # always save the mask

            if num_good_data < min_points:
                fit_dict['status'][ii,jj] = -1
                continue
            wave_ij = wave_ij[loc_good]
            inten_ij = inten_ij[loc_good]
            errs_ij = errs_ij[loc_good]

            # Scale guess parameters to data (should speed up fitting some)
            newguess = scale_guess(wave_ij, inten_ij, oldguess, n_gauss, n_poly)

            # Plug in new guess values to parinfo
            for i in range(len(newguess)):
                parinfo[i]['value'] = newguess[i]

            # Assemble dict of extra args to pass to mpfit
            fa = {'x': wave_ij, 'y': inten_ij, 'error': errs_ij,
                  'n_gauss': n_gauss, 'n_poly': n_poly}

            # Fit the profile
            out = mpfit(fit_fns.multigaussian_deviates, parinfo=parinfo,
                        functkw=fa, xtol=1.0E-6, ftol=1.0E-6, gtol=1.0E-6,
                        maxiter=2000, quiet=1)

            # check convergence status for a valid result
            if out.status > 0:
                # compute line inten and errors directly (may need to revisit)
                fpeaks = out.params[loc_peaks]
                fwdths = out.params[loc_wid]
                epeaks = out.perror[loc_peaks]
                ewdths = out.perror[loc_wid]
                l_inten = np.sqrt(2*np.pi)*fpeaks*fwdths
                e_inten = np.zeros(n_gauss)
                for n in range(n_gauss):
                    if fpeaks[n] != 0 and fwdths[n] != 0:
                        e_inten[n] = (l_inten[n]*np.sqrt((epeaks[n]/fpeaks[n])**2
                                                         +(ewdths[n]/fwdths[n])**2))
                    else:
                        e_inten[n] = 0.0

                # assemble fit structure
                fit_dict['status'][ii,jj] = out.status
                fit_dict['chi2'][ii,jj] = out.fnorm/out.dof
                fit_dict['wavelength'][ii,jj,:] = wave_cube[ii,jj,:]
                fit_dict['params'][ii,jj,:] = out.params
                fit_dict['perror'][ii,jj,:] = out.perror
                fit_dict['int'][ii,jj,:] = l_inten
                fit_dict['err_int'][ii,jj,:] = e_inten
            else:
                print(' ! fit did not converge!')
                fit_dict['status'][ii,jj] = out.status

    return fit_dict

def fit_spectra(inten, template, parinfo=None, wave=None, errs=None, min_points=7,
                ncpu='max', unsafe_mp=False, ignore_warnings=False,
                skip_fitting=False, debug=False):
    """Fit one or more EIS line spectra using mpfit (with multiprocessing).

    Parameters
    ----------
    inten : EISCube object, array_like, or filepath
        One or more intensity profiles to be fit. The code will loop over the data
        according to its dimensionality. 3D data is assumed to be a full EIS raster
        (or a sub region), 2D data is assumed to be a single EIS slit, and 1D data
        is assumed to be a single profile.
    template : EISFitTemplate object, dict, or filepath
        Either an EISFitTemplate, a 'template' dictionary, or the path to a
        template file.
    parinfo : list, optional
        List of dictionaries with fit parameters formatted for use with mpfit.
        Will supercede any parinfo lists loaded from an EISFitTemplate. Required
        if the 'template' parameter is given as a dictionary.
    wave : array_like, optional
        Associated wavelength values for the spectra. Required if 'inten' is
        given as an array and ignored otherwise.
    errs : array_like, optional
        Intensity error values for the spectra. Required if 'inten' is given as
        an array and ignored otherwise.
    min_points : int, optional
        Minimum number of good quality data points (i.e. non-zero values & errs)
        to be used in each fit. Spectra with fewer data points will be skipped.
        Must be a number >= the total number of fit parameters. Default is 7.
    ncpu : int, optional
        Number of cpu processes to parallelize over. Must be less than or equal
        to the total number of cores the system has. If set to 'max' or None, the
        code will use the maximum number of cores available. Default is 'max'.
        Important: due to the specifics of how the multiprocessing library works,
        any statements that call fit_spectra() using ncpu > 1 MUST be wrapped in
        a "if __name__ == '__main__':" statement in the top-level program. If such
        a "name guard" statement is not detected, this function will fall back to
        using a single process.
    unsafe_mp : bool, optional
        If set to True (and ncpu > 0), will use multiprocessing even if there is
        no "name guard" in use (see above). Used by the console script
        "eit_fit_files". Default is False (name guard enforced). Disabling the
        name guard runs the risk of spawning infinite processes if run incorrectly.
        USE AT YOUR OWN RISK!
    ignore_warnings : bool, optional
        If set to True, will silence the warning about a missing or disabled name
        guard (we are serious at it, be careful). Default is False.
    skip_fitting : bool, optional
        If set to True, will skip the fitting altogether and just return an empty
        EISFitResult instance. Used mainly for testing. Default is False.
    debug : bool, optional
        If set to True, will print some extra information useful for debugging
        development versions of the code. Default is False.

    Returns
    -------
    fit_res : EISFitResult class instance
        An EISFitResult object containing the output fit parameters.
    """
    # Validate template & parinfo and read / copy as needed
    if isinstance(template, (str, pathlib.Path)):
        template_obj = EISFitTemplate.read_template(template)
        if template_obj is None:
            return None
        template_copy = copy.deepcopy(template_obj.template)
        parinfo_copy = copy.deepcopy(template_obj.parinfo)
        tmplt_filename = template_obj.filename_temp
    elif isinstance(template, EISFitTemplate):
        template_copy = copy.deepcopy(template.template)
        parinfo_copy = copy.deepcopy(template.parinfo)
        tmplt_filename = template.filename_temp
    elif isinstance(template, dict):
        template_copy = copy.deepcopy(template)
        parinfo_copy = None
        tmplt_filename = 'unknown'
    else:
        print('Please input either the path to a template file, an'
             +' EISFitTemplate instance, or dictionary.', file=sys.stderr)
        return None

    if isinstance(parinfo, list):
        # Direct user-input always supercedes automatically loaded data
        parinfo_copy = copy.deepcopy(parinfo)
    elif parinfo_copy is None:
        print('Please input a parinfo list or a full EISFitTemplate.', file=sys.stderr)
        return None

    # If given an EISCube, extract data arrays and zero out masked values.
    # Otherwise, make local copies and just zero out negative values
    # TODO: Add more validation for the case of input arrays
    if isinstance(inten, (str, pathlib.Path)):
        central_wave = np.mean([template_copy['wmin'], template_copy['wmax']])
        eis_cube = read_cube(inten, window=float(central_wave))
        if eis_cube is None:
            return None
        wave_cube = eis_cube.wavelength.copy()
        errs_cube = eis_cube.uncertainty.array.copy()
        inten_cube = eis_cube.data.copy()
        metadata = copy.deepcopy(eis_cube.meta)
        data_units = eis_cube.unit.to_string()
        data_radcal = copy.deepcopy(eis_cube.radcal)
        loc_masked = np.where(eis_cube.mask == True)
        inten_cube[loc_masked] = 0
        errs_cube[loc_masked] = 0
        del eis_cube
    elif isinstance(inten, EISCube):
        wave_cube = inten.wavelength.copy()
        errs_cube = inten.uncertainty.array.copy()
        inten_cube = inten.data.copy()
        metadata = copy.deepcopy(inten.meta)
        data_units = inten.unit.to_string()
        data_radcal = copy.deepcopy(inten.radcal)
        loc_masked = np.where(inten.mask == True)
        inten_cube[loc_masked] = 0
        errs_cube[loc_masked] = 0
    elif isinstance(inten, np.ndarray):
        if not isinstance(wave, np.ndarray):
            print('Please input a wavelength array or a full EISCube.', file=sys.stderr)
            return None
        elif not isinstance(errs, np.ndarray):
            print('Please input an error array or a full EISCube.', file=sys.stderr)
            return None
        wave_cube = wave.copy()
        errs_cube = errs.copy()
        inten_cube = inten.copy()
        metadata = {'filename_data':'unknown', 'index':{}, 'pointing':{},
                    'radcal':'unknown', 'wave':'unknown'}
        data_units = 'unknown'
        data_radcal = 'unknown'
        loc_bad = np.where(errs_cube <= 0)
        inten_cube[loc_bad] = 0
        errs_cube[loc_bad] = 0
    else:
        print('Error: missing or invalid data. Please input a filepath, EISCube,'
             +' or complete set of intensity, wavelength, and error arrays.', file=sys.stderr)
        return None

    # Validate input ncpu value
    if str(ncpu).lower() == 'max' or str(ncpu).lower() == 'none':
        ncpu = mp.cpu_count()
    else:
        ncpu = int(ncpu)

    if ncpu <= 0:
        ncpu = 1
    elif ncpu > mp.cpu_count():
        ncpu = mp.cpu_count()

    # Ensure that multiprocessing will not spawn infinite child processes
    if ncpu > 1:
        name_guard = check_for_name_guard(debug=debug)
        if unsafe_mp == True and name_guard == False:
            if ignore_warnings == False:
                print('CRITICAL WARNING: unsafe_mp == True while no name guard'
                     +' was found in the top-level script! Be aware, parallel'
                     +' processes may freeze or behave unexpectedly.')
        elif name_guard == False:
            ncpu = 1
            if ignore_warnings == False:
                print('WARNING: no name guard was found in the top-level script!'
                     +' Falling back to a single process for safety.')

    # Check value of min_points
    n_params = len(parinfo_copy)
    if min_points is None or min_points < n_params:
        print('WARNING: min_points must be >= total number of fit parameters.'
             +' min_points has been set to '+str(n_params))
        min_points = n_params

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
    fit_res = EISFitResult(wave_cube, template_copy, parinfo_copy,
                           func_name='multigaussian', data_units=data_units,
                           radcal=data_radcal)
    fit_res.meta = metadata
    fit_res.meta['filename_template'] = tmplt_filename
    fit_res.fit_module = 'mpfit'
    fit_res.fit_method = 'LevMarLSQ'

    if skip_fitting != True:
        t1 = datetime.now() # start a simple timer
        nexp.value = n_steps
        print(f' + computing fits for {n_steps:d} exposures, each with {n_pxls:d} spectra')

        # Run fitting in either a single process (default) or using multiprocessing
        if ncpu == 1:
            print(' + running mpfit in a single process')
            fit_dict = fit_with_mpfit(wave_cube, inten_cube, errs_cube,
                                      template_copy, parinfo_copy, min_points,
                                      data_units=data_units)
            fit_res.fit = fit_dict
        else:
            if ncpu > n_steps:
                ncpu = n_steps
            print(f' + running mpfit on {ncpu:d} cores (of {mp.cpu_count():d})')

            # initialize pool of workers
            with mp.Pool(processes=ncpu) as pool:

                # Split out the data for each single slit and run the pool
                args = [(wave_cube[:,jj:jj+1,:], inten_cube[:,jj:jj+1,:], errs_cube[:,jj:jj+1,:],
                        template_copy, parinfo_copy, min_points, jj+1, data_units)
                        for jj in range(n_steps)]
                pool_out = pool.starmap(fit_with_mpfit, args)

            # Now, loop over each slit and copy fit results values to full output object
            for jj in range(n_steps):
                fit_res.fit['status'][:,jj] = pool_out[jj]['status'][:,0]
                fit_res.fit['chi2'][:,jj] = pool_out[jj]['chi2'][:,0]
                fit_res.fit['mask'][:,jj,:] = pool_out[jj]['mask'][:,0,:]
                fit_res.fit['wavelength'][:,jj,:] = pool_out[jj]['wavelength'][:,0,:]
                fit_res.fit['params'][:,jj,:] = pool_out[jj]['params'][:,0,:]
                fit_res.fit['perror'][:,jj,:] = pool_out[jj]['perror'][:,0,:]
                fit_res.fit['int'][:,jj,:] = pool_out[jj]['int'][:,0,:]
                fit_res.fit['err_int'][:,jj,:] = pool_out[jj]['err_int'][:,0,:]

        # Calculate the Doppler velocity for each line
        # TODO: revisit error estimation
        for gg in range(fit_res.n_gauss):
            base_wave = parinfo_copy[1+3*gg]['value']
            obs_cent = fit_res.fit['params'][:,:,1+3*gg]
            obs_errs = fit_res.fit['perror'][:,:,1+3*gg]
            velocity = calc_velocity(obs_cent, base_wave)
            fit_res.fit['vel'][:,:,gg] = velocity
            rel_err = obs_errs/obs_cent
            fit_res.fit['err_vel'][:,:,gg] = rel_err*velocity

        # print status
        t2 = datetime.now() # end timer
        num_fit = len(np.where(fit_res.fit['status'] > -0)[0])
        num_too_few = len(np.where((fit_res.fit['status'] == -2) |
                                   (fit_res.fit['status'] == -1))[0])
        num_bad_params = len(np.where((fit_res.fit['status'] == -3) |
                                      (fit_res.fit['status'] == 0))[0])
        print('\n')
        print('Finished computing fits!')
        print(f'   runtime : {t2-t1}')
        print(f'   {num_fit} spectra fit without issues')
        print(f'   {num_too_few} spectra have < {min_points} good data points')
        print(f'   {num_bad_params} spectra have bad or invalid parameters')

        # reset global counters
        cntr.value = 0
        nexp.value = 0

    return fit_res
