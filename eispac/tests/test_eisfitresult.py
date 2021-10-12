import pytest
import numpy as np

import eispac


@pytest.fixture
def empty_fit_result():
    return eispac.EISFitResult(empty=True)


def test_EISFitResult(test_data_filepath, test_template_filepath):
    tmplt = eispac.EISFitTemplate.read_template(test_template_filepath)
    eis_cube = eispac.read_cube(test_data_filepath, 192.394)
    wave_cube = eis_cube.wavelength
    fit_res = eispac.EISFitResult(wave_cube, tmplt.template, tmplt.parinfo,
                                  func_name='multigaussian')
    assert isinstance(fit_res, eispac.EISFitResult)


def test_empty_is_fit_result(empty_fit_result):
    assert isinstance(empty_fit_result, eispac.EISFitResult)


@pytest.mark.parametrize('components,answer', [
    (None, None),
    (-2, 0),
    (2, 'error'),
    ([0,1], [0,1]),
])
def test_validate_component_num(empty_fit_result, components, answer):
    compare = empty_fit_result._validate_component_num(components) == answer
    # wrapp in np.array().all() to enable array comparison
    assert np.array(compare).all()


@pytest.mark.parametrize('coords,answer', [
    (None, None),
    ([10, -18], [10, 4]),
    ([9001, 4], 'error'),
    ([1], 'error'),
])
def test_validate_coords(empty_fit_result, coords, answer):
    assert empty_fit_result._validate_coords(coords) == answer
