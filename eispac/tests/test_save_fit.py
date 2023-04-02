import eispac

"""
 A trivial change
"""

def test_invalid_save_dir():
    saved_files = eispac.save_fit('no result')
    assert saved_files is None

def test_invalid_save_dir(test_head_filepath):
    empty_fit = eispac.EISFitResult(empty=True)
    empty_fit.meta['filename_head'] = test_head_filepath
    saved_files = eispac.save_fit(empty_fit, save_dir=42)
    assert saved_files is None
