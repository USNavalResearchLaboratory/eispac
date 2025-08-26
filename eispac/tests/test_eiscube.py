import pathlib
import numpy as np
import pytest
from ndcube import NDCube
import astropy.units as u
from astropy.nddata import StdDevUncertainty

import eispac

@pytest.fixture
def ex_mod_index():
    output_hdr = dict()
    # Time information
    output_hdr['naxis'] = 3
    output_hdr['date_obs'] = '2021-03-06T06:44:44.000'
    output_hdr['date_beg'] = '2021-03-06T06:44:44.000'
    output_hdr['date_avg'] = '2021-03-06T06:47:09.000'
    output_hdr['date_end'] = '2021-03-06T06:49:34.000'
    output_hdr['timesys'] = 'UTC'
    # Observation details
    output_hdr['telescop'] = 'Hinode'
    output_hdr['instrume'] = 'EIS'
    output_hdr['line_id'] = 'Fe XII 192.410'
    output_hdr['measrmnt'] = 'intensity'
    output_hdr['bunit'] = 'erg / (cm2 s sr)'
    output_hdr['slit_id'] = '2"'
    output_hdr['slit_ind'] = 2
    output_hdr['obs_type'] = 'scan'
    output_hdr['nraster'] = 25 # 1 for sit-and-stare
    output_hdr['nexp'] = 25 # total num exposures
    output_hdr['nexp_prp'] = 1 # num exposures per slit position
    output_hdr['tr_mode'] = 'FIX' # tracking mode. "FIX" for no tracking
    output_hdr['saa'] = 'IN' # IN / OUT South Atlantic Anomaly (SAA)
    output_hdr['hlz'] = 'OUT' # IN / OUT High-Latitude Zone (HLZ)
    output_hdr['exptime'] = 10.000087
    output_hdr['cadence'] = 11.660708333333332
    output_hdr['timeunit'] = 's'
    # IDs and Study information
    output_hdr['tl_id'] = 76255
    output_hdr['jop_id'] = 396
    output_hdr['study_id'] = 589
    output_hdr['stud_acr'] = 'cam_flare_hot1_study'
    output_hdr['rast_id'] = 568
    output_hdr['rast_acr'] = 'cam_flare_hot1'
    output_hdr['obstitle'] = 'HOP396 Collaboration w/ Chandarayan-2 XSM AR12807'
    output_hdr['obs_dec'] = 'cam_flare_hot1_study'
    output_hdr['target'] = 'Flare'
    output_hdr['sci_obj'] = 'AR'
    output_hdr['noaa_num'] = 0
    # Coordinate information
    output_hdr['naxis1'] = 25
    output_hdr['cname1'] = 'Solar-X'
    output_hdr['crval1'] = -46.731937885284424
    output_hdr['crpix1'] = 1
    output_hdr['cdelt1'] = 3.9936
    output_hdr['ctype1'] = 'HPLN-TAN'
    output_hdr['cunit1'] = 'arcsec'
    output_hdr['naxis2'] = 120
    output_hdr['cname2'] = 'Solar-Y'
    output_hdr['crval2'] = -257.1441192626953
    output_hdr['crpix2'] = 1
    output_hdr['cdelt2'] = 1.0
    output_hdr['ctype2'] = 'HPLT-TAN'
    output_hdr['cunit2'] = 'arcsec'
    output_hdr['naxis3'] = 24
    output_hdr['cname3'] = 'Wavelength'
    output_hdr['crval3'] = 192.14506798083278
    output_hdr['crpix3'] = 1
    output_hdr['cdelt3'] = 0.02228723732869753
    output_hdr['ctype3'] = 'WAVE'
    output_hdr['cunit3'] = 'Angstrom'
    output_hdr['fovx'] = 99.83999729156494
    output_hdr['fovy'] = 120.0
    output_hdr['xcen'] = 3.188060760498047
    output_hdr['ycen'] = -197.1441192626953
    output_hdr['hgln_obs'] = 8.887297412477977e-15
    output_hdr['hglt_obs'] = -7.252204140746373
    output_hdr['dsun_obs'] = 148415597571.13638

    return output_hdr

@pytest.fixture
def ex_cube_meta(ex_mod_index):
    meta_dict = dict()
    meta_dict['wave'] = np.linspace(192.14012857, 192.65273503, num=ex_mod_index['naxis3'])
    meta_dict['radcal'] = 2*np.ones(ex_mod_index['naxis3'], dtype=np.float32)
    meta_dict['mod_index'] = ex_mod_index
    meta_dict['date_obs'] = np.array(['2021-03-06T06:49:23.857', '2021-03-06T06:49:12.256',
                                      '2021-03-06T06:49:00.609', '2021-03-06T06:48:48.826',
                                      '2021-03-06T06:48:37.240', '2021-03-06T06:48:25.561',
                                      '2021-03-06T06:48:13.951', '2021-03-06T06:48:02.344',
                                      '2021-03-06T06:47:50.664', '2021-03-06T06:47:38.889',
                                      '2021-03-06T06:47:27.207', '2021-03-06T06:47:15.496',
                                      '2021-03-06T06:47:03.840', '2021-03-06T06:46:52.244',
                                      '2021-03-06T06:46:40.477', '2021-03-06T06:46:28.830',
                                      '2021-03-06T06:46:17.170', '2021-03-06T06:46:05.473',
                                      '2021-03-06T06:45:53.834', '2021-03-06T06:45:42.254',
                                      '2021-03-06T06:45:30.471', '2021-03-06T06:45:18.813',
                                      '2021-03-06T06:45:07.221', '2021-03-06T06:44:55.631',
                                      '2021-03-06T06:44:44.000'], dtype='<U24')
    meta_dict['duration'] = np.array([9.999931,  9.999666,  9.999325, 10.000283, 
                                      10.000707, 10.000062, 9.999954, 10.000371,  
                                      9.999854,  9.999475, 10.000266, 10.00075 ,
                                      9.999315, 10.000467, 10.000361, 10.000407,  
                                      9.999714, 10.000589, 9.999717,  9.99959 , 
                                      10.00086 , 10.000018,  9.999368, 10.000409,
                                      10.000702], dtype=np.float32)
    meta_dict['notes'] = []

    return meta_dict

@pytest.fixture
def ex_eiscube(ex_mod_index, ex_cube_meta):
    wave_arr = np.tile(ex_cube_meta['wave'][np.newaxis, np.newaxis, :], 
                       (ex_mod_index['naxis2'], ex_mod_index['naxis1'], 1))
    data_arr = np.ones_like(wave_arr)
    err_arr = StdDevUncertainty(0.2*np.ones_like(wave_arr))
    ex_wcs = eispac.eiscube._make_wcs_from_cube_index(ex_mod_index)
    data_mask = data_arr <= 0 # nothing will be masked
    
    ex_cube = eispac.EISCube(data_arr, wcs=ex_wcs, uncertainty=err_arr,
                              wavelength=wave_arr, radcal=ex_cube_meta['radcal'],
                              meta=ex_cube_meta, unit='erg / (cm2 s sr)', 
                              mask=data_mask)
    return ex_cube

def test_slice_box(ex_eiscube):
    slice_cube = ex_eiscube[0:10,0:10,:]
    assert isinstance(slice_cube, eispac.EISCube)

def test_slice_step(ex_eiscube):
    slice_cube = ex_eiscube[:,1,:]
    assert isinstance(slice_cube, eispac.EISCube)

def test_slice_point(ex_eiscube):
    slice_cube = ex_eiscube[1,1,:]
    assert isinstance(slice_cube, eispac.EISCube)

def test_remove_radcal(ex_eiscube):
    # Note: the ex radcal is filled with 2.0 and the ex data is filled with 1.0
    #       so removing the radcal should give data values of 0.5
    photon_cube = ex_eiscube.remove_radcal()
    assert isinstance(photon_cube, eispac.EISCube)
    assert photon_cube.data[0,0,0] == 0.5

def test_apply_radcal(ex_eiscube):
    new_cal_cube = ex_eiscube.apply_radcal()
    assert isinstance(new_cal_cube, eispac.EISCube)

def test_shift_ref_coord(ex_eiscube):
    shift_cube = ex_eiscube.shift_reference_coord(10*u.arcsec, 10*u.arcsec)
    assert isinstance(shift_cube, eispac.EISCube)
    assert shift_cube.meta['mod_index']['crval1'] == ex_eiscube.meta['mod_index']['crval1'] + 10.0
    assert shift_cube.meta['mod_index']['crval2'] == ex_eiscube.meta['mod_index']['crval2'] + 10.0

def test_set_ref_date(ex_eiscube):
    ex_eiscube._set_reference_date('2010-10-10T10:10:10.100')
    assert ex_eiscube.meta['mod_index']['date_obs'] == '2010-10-10T10:10:10.100'

def test_total_intensity(ex_eiscube):
    # Note: the ex data is filled with 1.0, so summing over wavelength should
    #       equal the number of wavelength values in the ex data array
    sum_cube = ex_eiscube.sum_spectra()
    assert isinstance(sum_cube, NDCube)
    assert sum_cube.data[0,0] == float(sum_cube.meta['mod_index']['naxis3'])

def test_sum_map(ex_eiscube):
    sum_map = ex_eiscube.sum_spectra(return_map=True)
    assert isinstance(sum_map, eispac.EISMap)

def test_smooth_cube_unitless(ex_eiscube):
    # Note: the ex data is filled with 1.0, so smoothing will give data values
    #       equal to 1.0 (sum value/num_box_pixels for all pixels in the box)
    sm_cube = ex_eiscube.smooth_cube(width=3)
    assert isinstance(sm_cube, eispac.EISCube)
    assert sm_cube.data[10,10,10] == 1.0
    
    sm_cube = ex_eiscube.smooth_cube(width=[5,3])
    assert isinstance(sm_cube, eispac.EISCube)
    assert sm_cube.data[10,10,10] == 1.0

def test_smooth_cube_arcsec(ex_eiscube):
    # Note: no need to test the math, since that is covered in the above test
    sm_cube = ex_eiscube.smooth_cube(width=3*u.arcsec)
    assert isinstance(sm_cube, eispac.EISCube)
    
    sm_cube = ex_eiscube.smooth_cube(width=[5*u.arcsec,3*u.arcsec])
    assert isinstance(sm_cube, eispac.EISCube)

def test_extract_pix(ex_eiscube):
    point_cube = ex_eiscube.extract_points([[0,0], [1,1]], units='pixel')
    assert isinstance(point_cube, eispac.EISCube)

def test_extract_arcsec(ex_eiscube):
    point_cube = ex_eiscube.extract_points([[0,-200], [20,-160]], units='arcsec')
    assert isinstance(point_cube, eispac.EISCube)