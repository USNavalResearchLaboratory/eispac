import pytest
import eispac

# test_data_filepath = '../data/test/eis_20210306_064444.data.h5'
# test_fit_filepath = '../data/test/eis_20210306_064444.fe_12_192_394.1c-0.fit.h5'
test_data_filepath = eispac.data.test_data
test_fit_filepath = eispac.data.test_fit
empty_fit = eispac.EISFitResult(empty=True)
empty_fit.meta['filename_head'] = eispac.data.test_head

def test_invalid_save_dir():
    saved_files = eispac.save_fit('no result')
    assert saved_files is None

def test_invalid_save_dir():
    saved_files = eispac.save_fit(empty_fit, save_dir=42)
    assert saved_files is None
