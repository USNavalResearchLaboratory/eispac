import pathlib

import numpy as np
import pytest
import astropy.units as u
from sunpy.time import parse_time

import eispac

@pytest.fixture
def ex_eiscube(test_data_filepath):
    return eispac.read_cube(test_data_filepath, 192.41)

@pytest.fixture
def ex_ref_map(ex_eiscube):
    # Just the base cube shifted by [10,10] arcsec and summed over intensity
    shifted_cube = ex_eiscube.shift_reference_coord(10*u.arcsec, 10*u.arcsec)
    shifted_cube._set_reference_date(shifted_cube.meta['mod_index']['date_avg'])
    shifted_ref_map = shifted_cube.sum_spectra(return_map=True)
    return shifted_ref_map

def test_invalid_input_maps(ex_ref_map):
    assert eispac.coalign_observations(3.14, ex_ref_map) is None
    assert eispac.coalign_observations([ex_ref_map, 3.14], ex_ref_map) is None
    assert eispac.coalign_observations(ex_ref_map, 3.14) is None

def test_invalid_coalign_index(ex_ref_map):
    assert eispac.coalign_observations(ex_ref_map, ex_ref_map, '0') is None
    assert eispac.coalign_observations(ex_ref_map, ex_ref_map, 99) is None

def test_unknown_method(ex_ref_map):
    assert eispac.coalign_observations(ex_ref_map, ex_ref_map, method='unknown') is None

def test_coalign(ex_eiscube, ex_ref_map):
    new_cube = eispac.coalign_observations(ex_eiscube, ex_ref_map)
    assert isinstance(new_cube, eispac.EISCube)
    assert new_cube.meta['mod_index']['crval1'] == ex_eiscube.meta['mod_index']['crval1'] + 10
    assert new_cube.meta['mod_index']['date_obs'] == ex_eiscube.meta['mod_index']['date_avg']

# def test_details(ex_eiscube, ex_ref_map):
#     new_cube, details = eispac.coalign_observations(ex_eiscube, ex_ref_map, return_details=True)
#     assert isinstance(details, dict)