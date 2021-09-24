import pathlib

import numpy as np
import pytest

import eispac

empty_inten = np.zeros((5,3,10))
empty_wave = np.zeros((5,3,10))
empty_errs = np.zeros((5,3,10))


@pytest.fixture
def test_template(test_template_filepath):
    return eispac.EISFitTemplate.read_template(test_template_filepath)


def test_invalid_template():
    fit_res = eispac.fit_spectra(empty_inten, 0)
    assert fit_res is None


def test_dict_template_and_invalid_parinfo(test_template):
    fit_res = eispac.fit_spectra(empty_inten, test_template.template, parinfo='invalid')
    assert fit_res is None


def test_template_str_filepath(test_template_filepath):
    fit_res = eispac.fit_spectra(empty_inten, test_template_filepath, wave=empty_wave,
                                 errs=empty_errs, skip_fitting=True)
    assert isinstance(fit_res, eispac.EISFitResult)


def test_template_pathlib_filepath(test_data_filepath, test_template_filepath):
    tmplt_path_obj = pathlib.Path(test_template_filepath)
    eis_cube = eispac.read_cube(test_data_filepath, 192.394)
    sub_raster = eis_cube[0:3, 0:3, :]
    fit_res = eispac.fit_spectra(sub_raster, tmplt_path_obj, skip_fitting=True)
    assert isinstance(fit_res, eispac.EISFitResult)


def test_missing_wavelengths(test_template_filepath):
    fit_res = eispac.fit_spectra(empty_inten, test_template_filepath)
    assert fit_res is None


def test_missing_errs(test_template_filepath):
    fit_res = eispac.fit_spectra(empty_inten, test_template_filepath, wave=empty_wave)
    assert fit_res is None


def test_fit_spectra_single_thread(test_data_filepath, test_template):
    eis_cube = eispac.read_cube(test_data_filepath, 192.394)
    sub_raster = eis_cube[0:3, 0:3, :]
    fit_res = eispac.fit_spectra(sub_raster, test_template)
    assert isinstance(fit_res, eispac.EISFitResult)


def test_fit_spectra_parallel(test_template, test_data_filepath):
    if __name__ == '__main__':
        eis_cube = eispac.read_cube(test_data_filepath, 192.394)
        sub_raster = eis_cube[0:3, 0:3, :]
        fit_res = eispac.fit_spectra(sub_raster, test_template, ncpu=4)
        assert isinstance(fit_res, eispac.EISFitResult)


def test_fit_spectra_parallel_without_guard(test_template, test_data_filepath):
    eis_cube = eispac.read_cube(test_data_filepath, 192.394)
    sub_raster = eis_cube[0:3, 0:3, :]
    fit_res = eispac.fit_spectra(sub_raster, test_template, ncpu=4)
    assert isinstance(fit_res, eispac.EISFitResult)
