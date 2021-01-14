import pytest
import eispac

test_data_filepath = '../data/test/test_dataset_eis_20190404_131513.data.h5'
test_tmplt_filepath = '../templates/eis_template_dir/fe_12_195_119.2c.template.h5'

def test_empty():
    fit_res = eispac.EISFitResult(empty=True)
    assert isinstance(fit_res, eispac.EISFitResult)

def test_EISFitResult():
    tmplt = eispac.read_template(test_tmplt_filepath)
    eis_cube = eispac.read_cube(test_data_filepath, 195.119)
    wave_cube = eis_cube.wavelength
    fit_res = eispac.EISFitResult(wave_cube, tmplt.template, tmplt.parinfo,
                                  func_name='multigaussian')
    assert isinstance(fit_res, eispac.EISFitResult)

def test_validate_component_num_none():
    fit_res = eispac.EISFitResult(empty=True)
    test_comp = fit_res._validate_component_num(None)
    assert test_comp == None

def test_validate_component_num_neg_index():
    fit_res = eispac.EISFitResult(empty=True)
    test_comp = fit_res._validate_component_num(-2)
    assert test_comp == 0

def test_validate_component_num_invalid():
    fit_res = eispac.EISFitResult(empty=True)
    test_comp = fit_res._validate_component_num(2)
    assert test_comp == 'error'

def test_validate_component_num_multiple():
    fit_res = eispac.EISFitResult(empty=True)
    test_comp = fit_res._validate_component_num([0,1])
    assert all(test_comp == [0,1]) == True

def test_validate_coords_none():
    fit_res = eispac.EISFitResult(empty=True)
    test_coords = fit_res._validate_coords(None)
    assert test_coords is None

def test_validate_coords_neg_index():
    fit_res = eispac.EISFitResult(empty=True)
    test_coords = fit_res._validate_coords([10,-18])
    assert test_coords == [10, 4]

def test_validate_coords_invalid():
    fit_res = eispac.EISFitResult(empty=True)
    test_coords = fit_res._validate_coords([9001,4])
    assert test_coords == 'error'

def test_validate_coords_wrong_shape():
    fit_res = eispac.EISFitResult(empty=True)
    test_coords = fit_res._validate_coords([1])
    assert test_coords == 'error'
