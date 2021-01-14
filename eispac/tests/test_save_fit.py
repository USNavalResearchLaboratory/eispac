import pytest
import eispac

test_data_filepath = '../data/test/test_dataset_eis_20190404_131513.data.h5'
test_fit_filepath = '../data/test/eis_20190404_131513_fe_12_195_119_c0.fit.h5'
empty_fit = eispac.EISFitResult(empty=True)
empty_fit.meta['filename_head'] = test_data_filepath

def test_invalid_save_dir():
    saved_files = eispac.save_fit('no result')
    assert saved_files is None

def test_invalid_save_dir():
    saved_files = eispac.save_fit(empty_fit, save_dir=42)
    assert saved_files is None
