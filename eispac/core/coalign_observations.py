__all__ = ['coalign_observations']

import sys
import copy
import numpy as np
from scipy.signal import correlate

import astropy.units as u
import sunpy.map
from ndcube import NDCubeSequence

from eispac.core.eiscube import EISCube

# Optional imports (not required by EISPAC)
try:
    from skimage.feature import match_template
    scikit_available = True
except:
    scikit_available = False

# The helper functions below were copied from sunkit_image v0.6.1
# Note: these functions were made internal in sunkit v0.5 (2023-08-10).
#       The offical recommendation is to copy these functions, since the sunkit
#       API may change. These functions have no use to EISPAC outside this
#       coalignement tool, so they are copied here rather than placing them in
#       the "extern" sub-module

################################################################################
##### SunKit-Image BSD-2-Clause license ########################################
################################################################################
# Copyright (c) 2024, The SunPy Community
#
# Redistribution and use in source and binary forms, with or without modification, 
# are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright notice, this 
#   list of conditions and the following disclaimer.
# * Redistributions in binary form must reproduce the above copyright notice, 
#   this list of conditions and the following disclaimer in the documentation 
#   and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

def _find_best_match_location(corr):
    """
    Calculate an estimate of the location of the peak of the correlation result
    in image pixels.

    Parameters
    ----------
    corr : `numpy.ndarray`
        A 2D correlation array.

    Returns
    -------
    `~astropy.units.Quantity`
        The shift amounts ``(y, x)`` in image pixels. Subpixel values are
        possible.
    """
    # Get the index of the maximum in the correlation function
    ij = np.unravel_index(np.argmax(corr), corr.shape)
    cor_max_x, cor_max_y = ij[::-1]

    # Get the correlation function around the maximum
    if len(corr) > 1:
        array_maximum = corr[
            np.max([0, cor_max_y - 1]) : np.min([cor_max_y + 2, corr.shape[0] - 1]),
            np.max([0, cor_max_x - 1]) : np.min([cor_max_x + 2, corr.shape[1] - 1]),
        ]
    else:
        # For single-valued corr arrays 
        # This can occur in some edge cases (usually just in testing)
        array_maximum = np.zeros([3,3])
        array_maximum[1,1] = corr[0,0]
    
    y_shift_maximum, x_shift_maximum = _get_correlation_shifts(array_maximum)

    # Get shift relative to correlation array
    y_shift_correlation_array = y_shift_maximum + cor_max_y * u.pix
    x_shift_correlation_array = x_shift_maximum + cor_max_x * u.pix

    return y_shift_correlation_array, x_shift_correlation_array

def _get_correlation_shifts(array):
    """
    Estimate the location of the maximum of a fit to the input array. The
    estimation in the "x" and "y" directions are done separately. The location
    estimates can be used to implement subpixel shifts between two different
    images.

    Parameters
    ----------
    array : `numpy.ndarray`
        An array with at least one dimension that has three elements. The
        input array is at most a 3x3 array of correlation values calculated
        by matching a template to an image.

    Returns
    -------
    `~astropy.units.Quantity`
        The ``(y, x)`` location of the peak of a parabolic fit, in image pixels.
    """
    # Check input shape
    ny = array.shape[0]
    nx = array.shape[1]
    if nx > 3 or ny > 3:
        msg = "Input array dimension should not be greater than 3 in any dimension."
        raise ValueError(msg)

    # Find where the maximum of the input array is
    ij = np.unravel_index(np.argmax(array), array.shape)
    x_max_location, y_max_location = ij[::-1]

    # Estimate the location of the parabolic peak if there is enough data.
    # Otherwise, just return the location of the maximum in a particular
    # direction.
    y_location = _parabolic_turning_point(array[:, x_max_location]) if ny == 3 else 1.0 * y_max_location

    x_location = _parabolic_turning_point(array[y_max_location, :]) if nx == 3 else 1.0 * x_max_location

    return y_location * u.pix, x_location * u.pix

def _parabolic_turning_point(y):
    """
    Find the location of the turning point for a parabola ``y(x) = ax^2 + bx +
    c``, given input values ``y(-1), y(0), y(1)``. The maximum is located at
    ``x0 = -b / 2a``. Assumes that the input array represents an equally spaced
    sampling at the locations ``y(-1), y(0) and y(1)``.

    Parameters
    ----------
    y : `numpy.ndarray`
        A one dimensional numpy array of shape "3" with entries that sample the
        parabola at "-1", "0", and "1".

    Returns
    -------
    `float`
        A float, the location of the parabola maximum.
    """
    numerator = -0.5 * y.dot([-1, 0, 1])
    denominator = y.dot([1, -2, 1])
    return numerator / denominator

# End of Sunkit-Image helper functions

def coalign_observations(eis_observations, reference_map, coalign_index=0, 
                         weight_radius=20, method='scipy_correlate',
                         return_details=False):
    """Coalign one or more EIS observations with a reference map (e.g. SDO/AIA)

    Calculate and apply coordinate shifts required to align one or more EISMaps
    or EISCubes with a reference map. A few different methods are available, 
    but they generally all use some variation of 2D cross-correlation in 
    fourier space.
    
    Note: only translational shifts are applied. This function does not compute
    or apply rotational shifts that may be required for higher-quality alignments.

    Parameters
    ----------
    eis_observations : `EISMap`, `EISCube`, list, `~sunpy.map.MapSequence`, or `~ndcube.NDCubeSequence`
        One or more EISMaps or EISCubes to be coaligned. Multiple observations 
        may be input as either a list, a `~sunpy.map.MapSequence`, or `~ndcube.NDCubeSequence`, 
        however only one obs will be compared to the reference map and used for 
        coalignment. 
    reference_map : `~sunpy.map.Map`
        Reference map to coalign the EIS obs to. Should be the same size or 
        larger than all of the EIS observations and be taken at roughly the same 
        time as the observation used for coalignment.
    coalign_index : int, optional
        Index of EIS observation used for computing coalignment. Only used if 
        multiple obs are input as a list, MapSequence, or NDCubeSequence. Must 
        have a value less than the total number of eis_observations (negative 
        indices are fine). Default is "0".
    weight_radius : int or float, optional
        Radial distance (in units of [Arcsec]) around the input EIS pointing used
        for weighting the correlation array. Locations inside this distance are 
        given a full weight of "1" when locating the best coalignment.
        This is particularly important when using the "scipy_correlate" method,
        which has a known bias for bright regions. Setting a radius of "None"
        or another non-positive value will disable weighting. Default is 20.
    method : str, optional
        Method used to compute the coalignment. All methods other than the 
        default require that the optional package "scikit-image" is installed. 
        The following methods are available:

        `scipy_correlate` (default) - 2D cross-correlation using an FFT
        `match_template` - 2D normalized correlation using FFT convolution 
    return_details :  bool, optional
        If set to "True", will also return a dict containing more deatils about
        the coalignment (e.g. shifts, correlation array). Default is "False".
        
    Returns
    -------
    output_obs : `EISMap`, `EISCube`, list, `~sunpy.map.MapSequence`, `~ndcube.NDCubeSequence`
        Copy of input obs with aligned and updated coordinate information
    details_dict : dict, optional
        Dict containing more details about the coalignment processes such as the
        shifts in [arcsec], pixel location of max correlation, the correlation
        and weighting arrays, as well as names and timestamps of the maps used.
        Only returned in "return_details=True". 
    """
    # Validate input maps and make copies (so we do not overwrite the inputs)
    if isinstance(eis_observations, (sunpy.map.GenericMap, EISCube)):
        output_type = 'single'
        output_obs = [copy.deepcopy(eis_observations)] # list for processing only
    elif isinstance(eis_observations, list):
        output_type = 'list'
        output_obs = copy.deepcopy(eis_observations)
        for OBS in output_obs:
            if not isinstance(OBS, (sunpy.map.GenericMap, EISCube)):
                print('ERROR: Invalid list element dtype. All obs in the list'
                     +' must be an EISMap or EISCube.', file=sys.stderr)
                return None
    elif isinstance(eis_observations, sunpy.map.MapSequence):
        output_type = 'map_sequence'
        output_obs = list(copy.deepcopy(eis_observations))
    elif isinstance(eis_observations, NDCubeSequence):
        output_type = 'cube_sequence'
        output_obs = list(copy.deepcopy(eis_observations))
    else:
        print('ERROR: Invalid "eis_observations" dtype. Please input an EISMap,'
             +' EISCube, list, MapSequence, or NDCubeSequence.', file=sys.stderr)
        return None
        
    # Validate reference map
    # TO-DO: maybe add option to automatically download an AIA 193 image?
    if isinstance(reference_map, sunpy.map.GenericMap):
        ref_map = copy.deepcopy(reference_map)
    else:
        print('ERROR: Invalid "reference_map" dtype. Please input a Sunpy Map.',
              file=sys.stderr)
        return None
        
    # Validate index of template map used for coalignment
    n_output_obs = len(output_obs)
    t_ind = coalign_index
    if ((not isinstance(t_ind, int)) 
    or (t_ind < 0 and np.abs(t_ind) > n_output_obs) or (t_ind >= n_output_obs)):
        print(f'ERROR: coalign_index={coalign_index} is invalid. Please input an'
             +' integer with a value less than the number of eis_observations.', 
             file=sys.stderr)
        return None
    
    # Extract coalignment map (or generate it if given an EISCube)
    coalign_map = copy.deepcopy(output_obs[t_ind])
    if isinstance(coalign_map, sunpy.map.GenericMap):
        coalign_map = coalign_map
    elif isinstance(coalign_map, EISCube):
        cube_line_id = coalign_map.meta['mod_index']['line_id']
        cube_line_ion = cube_line_id.split()[0]
        cube_line_wave = float(cube_line_id.split()[-1])
        if cube_line_ion.lower().startswith('ccd'):
            # Sum "full" CCD image around the ref map wavelength (if valid)
            ref_wave = ref_map.wavelength.to('angstrom').value
            if (ref_wave >= coalign_map.meta['wave'][0] 
            and ref_wave <= coalign_map.meta['wave'][-1]): 
                sum_range = [ref_wave, 2.5]*u.Angstrom
            else:
                sum_range = [cube_line_wave, 2.5]*u.Angstrom
            coalign_map = coalign_map.sum_spectra(sum_range, return_map=True)
        else:
            # Sum entire wavelength axis of template EISCube
            coalign_map = coalign_map.sum_spectra(return_map=True)

    # Check map wavelengths and timestamps and issue warnings as needed
    tmplt_wave = coalign_map.wavelength.to('angstrom').value
    ref_wave = ref_map.wavelength.to('angstrom').value
    if tmplt_wave is None or ref_wave is None or abs(ref_wave - tmplt_wave) > 3:
        print('WARNING: the reference and coalignment maps have wavelengths'
             +' that are either unknown or more than 3 Angstroms apart.'
             +' The coalignment results may be inaccurate.')
        
    tmplt_start = coalign_map.date_start 
    tmplt_start = tmplt_start if tmplt_start is not None else coalign_map.date
    tmplt_end = coalign_map.date_end
    tmplt_end = tmplt_end if tmplt_end is not None else coalign_map.date
    ref_date = ref_map.date_average if ref_map.date_average is not None else ref_map.date
    if tmplt_start - 10*u.min > ref_date or tmplt_end + 10*u.min < ref_date:
        print('WARNING: the reference map has a timestamp that is more than'
             +' 10 min outside the the observation time of the coalignment map.'
             +' The coalignment results may be inaccurate.')
        
    # Check for nans or inf values and issue warnings as needed
    if not np.all(np.isfinite(coalign_map.data)):
        print('WARNING: coalignment map has nonfinite values such as NaN or Inf!'
             +' The coalignment results may be inaccurate. Please consider'
             +' replacing nonfinite values with either zeros or a local mean.')
        
    if not np.all(np.isfinite(ref_map.data)):
        print('WARNING: reference map has nonfinite values such as NaN or Inf!'
             +' The coalignment results may be inaccurate. Please consider'
             +' replacing nonfinite values with either zeros or a local mean.')
        
    # Validate method string
    if method.lower() not in ['scipy_correlate', 'match_template']:
        print(f'ERROR: {method} is not a valid coalignment method.', sys.stderr)
        return None
    elif (not scikit_available) and (method.lower() in ['match_template']):
        print(f'WARNING: ths scikit-image package was not found. Cannot use the'
             +f' {method} method; will use "scipy_correlate" instead.')
        method = 'scipy_correlate'

    # Rescale ref map to the same X and Y scaling as the template map
    new_nx = (ref_map.scale.axis1 * ref_map.dimensions.x) / coalign_map.scale.axis1
    new_ny = (ref_map.scale.axis2 * ref_map.dimensions.y) / coalign_map.scale.axis2
    ref_map = ref_map.resample(u.Quantity([new_nx, new_ny]))

    # Compute the cross-correlation array
    tmplt_data = coalign_map.data
    ref_data = ref_map.data
    if method.lower() == 'scipy_correlate':
        # It is best if we first normalize the data arrays
        # TO-DO: consider ignoring values <= 0 for the mean and std
        # tmplt_stat_arr = tmplt_data[(tmplt_data > 0)] # ignore negative values
        # ref_stat_arr = ref_data[(ref_data > 0)]
        tmplt_data = (tmplt_data - np.nanmean(tmplt_data )) / np.nanstd(tmplt_data)
        ref_data = (ref_data - np.nanmean(ref_data)) / np.nanstd(ref_data)
        corr = correlate(ref_data, tmplt_data, method='fft', mode='valid')
    elif method.lower() == 'match_template':
        corr = match_template(ref_data, tmplt_data, pad_input=False)
        
    # Apply weighting to the correlation arrays
    # Note: "dist" below is in units of [arcsec] from the input EIS pointing
    if isinstance(weight_radius, (int, float)) and weight_radius > 0:
        w_radius = float(weight_radius)
        mesh_x, mesh_y = np.meshgrid(np.arange(corr.shape[1]), 
                                    np.arange(corr.shape[0]), indexing='xy')
        x_scale = ref_map.scale.axis1.value # units of [arcsec/pix]
        y_scale = ref_map.scale.axis2.value
        guess_pix = ref_map.world_to_pixel(coalign_map.bottom_left_coord)
        dist = np.sqrt((x_scale*(mesh_x-guess_pix.x.value))**2 
                    + (y_scale*(mesh_y-guess_pix.y.value))**2)
        weights = np.zeros_like(dist)
        weights[(dist <= w_radius)] = 1
        weights[(dist > w_radius)] = w_radius/dist[(dist > w_radius)]
    else:
        weights = np.ones_like(corr)

    # Compute the sub-pixel location of best correlation
    best_ypix, best_xpix = _find_best_match_location(corr*weights)

    # Compute the relative shifts in units of [arcsec]
    best_coord = ref_map.pixel_to_world(best_xpix, best_ypix)
    Tx_shift = best_coord.Tx - coalign_map.bottom_left_coord.Tx
    Ty_shift = best_coord.Ty - coalign_map.bottom_left_coord.Ty

    # Apply the shifts to all input obs
    # Note: we also adjust the "date_obs" of each object since coaligning also
    #       implicitly shifts the timestamp of the coordinate frame
    for i in range(n_output_obs):
        output_obs[i] = output_obs[i].shift_reference_coord(Tx_shift, Ty_shift)
        output_obs[i]._set_reference_date(ref_date)
    
    # Repack or extract coaligned objects or sequences
    if output_type == 'single':
        output_obs = output_obs[0]
    elif output_type == 'map_sequence':
        output_obs = sunpy.map.MapSequence(output_obs, sortby=None)
    elif output_type == 'cube_sequence':
        output_obs = NDCubeSequence(output_obs, meta=eis_observations.meta, 
                                    common_axis=eis_observations._common_axis)
    
    if return_details:
        details_dict = {'xy_shift':[Tx_shift.to('arcsec').value, 
                                    Ty_shift.to('arcsec').value]*u.arcsec, 
                        'xy_best_pix':[best_xpix.value, best_ypix.value]*u.pix,
                        'coalign_name':coalign_map.name,
                        'coalign_date_avg':coalign_map.date_average,
                        'coalign_index':coalign_index,
                        'ref_name':ref_map.name,
                        'ref_date_avg':ref_date,
                        'corr':corr, 
                        'weights':weights}
            
        return output_obs, details_dict
    else:
        return output_obs