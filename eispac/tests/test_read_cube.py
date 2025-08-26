import pathlib

from ndcube import NDCube

import eispac


def test_invalid_filepath_type():
    eis_cube = eispac.read_cube(42)
    assert eis_cube is None

def test_missing_files():
    eis_cube = eispac.read_cube('non-existent_file.head.h5', 0)
    assert eis_cube is None

def test_invalid_window(test_data_filepath):
    eis_cube = eispac.read_cube(test_data_filepath, -7)
    assert eis_cube is None

def test_read_cube_data_str_filepath(test_data_filepath):
    eis_cube = eispac.read_cube(test_data_filepath, 192.394)
    assert isinstance(eis_cube, eispac.EISCube)

def test_read_cube_head_str_filepath(test_head_filepath):
    eis_cube = eispac.read_cube(test_head_filepath, 0)
    assert isinstance(eis_cube, eispac.EISCube)

def test_read_fit_pathlib_filepath(test_data_filepath):
    path_obj = pathlib.Path(test_data_filepath)
    eis_cube = eispac.read_cube(path_obj)
    assert isinstance(eis_cube, eispac.EISCube)