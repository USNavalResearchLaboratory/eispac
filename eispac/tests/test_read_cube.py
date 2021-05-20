import pytest
import eispac

# test_data_filepath = '../data/test/eis_20210306_064444.data.h5'
# test_head_filepath = '../data/test/eis_20210306_064444.head.h5'
test_data_filepath = eispac.data.test_data
test_head_filepath = eispac.data.test_head

def test_invalid_filepath_type():
    eis_cube = eispac.read_cube(42)
    assert eis_cube is None

def test_missing_files():
    eis_cube = eispac.read_cube('non-existent_file.head.h5', 0)
    assert eis_cube is None

def test_invalid_window():
    eis_cube = eispac.read_cube(test_data_filepath, -7)
    assert eis_cube is None

def test_read_cube_data_str_filepath():
    eis_cube = eispac.read_cube(test_data_filepath, 192.394)
    assert isinstance(eis_cube, eispac.EISCube)

def test_read_cube_head_str_filepath():
    eis_cube = eispac.read_cube(test_head_filepath, 0)
    assert isinstance(eis_cube, eispac.EISCube)

def test_read_fit_pathlib_filepath():
    import pathlib
    path_obj = pathlib.Path(test_data_filepath)
    eis_cube = eispac.read_cube(path_obj)
    assert isinstance(eis_cube, eispac.EISCube)

def test_total_intensity():
    from ndcube import NDCube
    eis_cube = eispac.read_cube(test_data_filepath, 192.394)
    sum_inten = eis_cube.sum_spectra()
    assert isinstance(sum_inten, NDCube)
