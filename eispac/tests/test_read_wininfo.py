import pathlib
import numpy as np
import eispac


def test_read_wininfo_invalid_type():
    assert eispac.read_wininfo(123) is None
    assert eispac.read_wininfo(None) is None
    assert eispac.read_wininfo(['sample.head.h5']) is None

def test_read_wininfo_nonexistent_file():
    assert eispac.read_wininfo('non-existent_file.head.h5') is None

def test_read_wininfo_wrong_filetype(test_toml_template_filepath):
    wininfo = eispac.read_wininfo(test_toml_template_filepath)
    assert wininfo is None

def test_read_wininfo_from_head(test_head_filepath):
    wininfo = eispac.read_wininfo(test_head_filepath)
    assert isinstance(wininfo, np.recarray)
    assert len(wininfo) > 0

def test_read_wininfo_from_data(test_data_filepath):
    wininfo = eispac.read_wininfo(test_data_filepath)
    assert isinstance(wininfo, np.recarray)
    assert len(wininfo) > 0

def test_read_wininfo_pathlib(test_head_filepath):
    wininfo = eispac.read_wininfo(pathlib.Path(test_head_filepath))
    assert isinstance(wininfo, np.recarray)
    assert len(wininfo) > 0

def test_read_wininfo_fields_and_values(test_head_filepath):
    wininfo = eispac.read_wininfo(test_head_filepath)
    expected_fields = {'iwin', 'line_id', 'wvl_min', 'wvl_max', 'nl', 'xs'}

    # Check field names
    assert expected_fields.issubset(set(wininfo.dtype.names))

    # Check first window contents and types
    first_win = wininfo[0]
    assert isinstance(first_win.iwin, (int, np.integer))
    assert isinstance(first_win.line_id, str)
    assert isinstance(first_win.wvl_min, (float, np.floating))
    assert isinstance(first_win.wvl_max, (float, np.floating))
    assert first_win.wvl_max > first_win.wvl_min
    assert first_win.nl > 0