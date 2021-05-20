import pathlib as _pathlib # should be hidden

test_data = str(_pathlib.Path(__file__).parent.joinpath('test/eis_20210306_064444.data.h5'))
test_head = str(_pathlib.Path(__file__).parent.joinpath('test/eis_20210306_064444.head.h5'))
test_fit = str(_pathlib.Path(__file__).parent.joinpath('test/eis_20210306_064444.fe_12_192_394.1c-0.fit.h5'))
