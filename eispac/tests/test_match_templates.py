import pathlib
import eispac


# Note: string and pathlib inputs are passed to eispac.read_wininfo
#       As such, we don't need to test all of the input types again here.

def test_invalid_input_type():
    assert eispac.match_templates(42) is None
    assert eispac.match_templates(None) is None
    assert eispac.match_templates(['sample.head.h5']) is None


def test_match_from_hdf5_file(test_head_filepath):
    matches = eispac.match_templates(test_head_filepath)
    # Note: output should be a list of lists
    assert isinstance(matches, list)
    assert len(matches) > 0
    assert isinstance(matches[0], list)
    num_matches = [len(m) for m in matches]
    assert all(num_matches) # at least 1 match for each window
    assert isinstance(matches[0][0], pathlib.Path)
    assert matches[0][0].is_file()


def test_match_from_eiscube(test_data_filepath):
    eis_cube = eispac.read_cube(test_data_filepath, window=192.394)
    matches = eispac.match_templates(eis_cube)
    # Note: Output should be a single list
    assert isinstance(matches, list)
    assert len(matches) > 0
    assert isinstance(matches[0], pathlib.Path)
    assert matches[0].is_file()