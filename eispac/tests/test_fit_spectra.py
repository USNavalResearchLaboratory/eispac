import numpy as np
import pytest
import eispac

# test_data_filepath = '../data/test/eis_20210306_064444.data.h5'
# test_tmplt_filepath = '../templates/eis_template_dir/fe_12_192_394.1c.template.h5'
test_data_filepath = eispac.data.test_data
test_tmplt_filepath = eispac.templates.template_dir+'/fe_12_192_394.1c.template.h5'

empty_inten = np.zeros((5,3,10))
empty_wave = np.zeros((5,3,10))
empty_errs = np.zeros((5,3,10))

def test_invalid_template():
    fit_res = eispac.fit_spectra(empty_inten, 0)
    assert fit_res is None

def test_dict_template_and_invalid_parinfo():
    tmplt = eispac.read_template(test_tmplt_filepath)
    tmplt_dict = tmplt.template
    fit_res = eispac.fit_spectra(empty_inten, tmplt_dict, parinfo='invalid')
    assert fit_res is None

def test_template_str_filepath():
    fit_res = eispac.fit_spectra(empty_inten, test_tmplt_filepath, wave=empty_wave,
                          errs=empty_errs, skip_fitting=True)
    assert isinstance(fit_res, eispac.EISFitResult)

def test_template_pathlib_filepath():
    import pathlib
    tmplt_path_obj = pathlib.Path(test_tmplt_filepath)
    eis_cube = eispac.read_cube(test_data_filepath, 192.394)
    sub_raster = eis_cube[0:3, 0:3, :]
    fit_res = eispac.fit_spectra(sub_raster, test_tmplt_filepath, skip_fitting=True)
    assert isinstance(fit_res, eispac.EISFitResult)

def test_missing_wavelengths():
    fit_res = eispac.fit_spectra(empty_inten, test_tmplt_filepath)
    assert fit_res is None

def test_missing_errs():
    fit_res = eispac.fit_spectra(empty_inten, test_tmplt_filepath, wave=empty_wave)
    assert fit_res is None

def test_fit_spectra_single_thread():
    tmplt = eispac.read_template(test_tmplt_filepath)
    eis_cube = eispac.read_cube(test_data_filepath, 192.394)
    sub_raster = eis_cube[0:3, 0:3, :]
    fit_res = eispac.fit_spectra(sub_raster, tmplt)
    assert isinstance(fit_res, eispac.EISFitResult)

def test_fit_spectra_parallel():
    if __name__ == '__main__':
        tmplt = eispac.read_template(test_tmplt_filepath)
        eis_cube = eispac.read_cube(test_data_filepath, 192.394)
        sub_raster = eis_cube[0:3, 0:3, :]
        fit_res = eispac.fit_spectra(sub_raster, tmplt, ncpu=4)
        assert isinstance(fit_res, eispac.EISFitResult)

def test_fit_spectra_parallel_without_guard():
    tmplt = eispac.read_template(test_tmplt_filepath)
    eis_cube = eispac.read_cube(test_data_filepath, 192.394)
    sub_raster = eis_cube[0:3, 0:3, :]
    fit_res = eispac.fit_spectra(sub_raster, tmplt, ncpu=4)
    assert isinstance(fit_res, eispac.EISFitResult)
