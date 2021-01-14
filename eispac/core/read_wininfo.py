__all__ = ['read_wininfo']

import sys
import pathlib
import numpy as np
import h5py

def read_wininfo(filename):

    # Input type validation (value checks are implemented later)
    if not isinstance(filename, (str, pathlib.Path)):
        print('Error: Please input a valid filepath as '
             +'either a string or pathlib.Path object', file=sys.stderr)
        return None

    # Parse filename and determine the directory and filename
    abs_filepath = pathlib.Path(filename).resolve()
    input_name = str(abs_filepath.name)
    input_dir = abs_filepath.parent
    if str(input_dir) == '.':
        input_dir = pathlib.Path().cwd()

    # Determine header filename, even if given a data filename by mistake
    head_filepath = input_dir.joinpath(input_name.replace('.data.h5', '.head.h5'))

    # Check that the file exists and has the correct file extension
    if not head_filepath.is_file():
        print('Error: Header file does not exist, ' + str(head_filepath),
              file=sys.stderr)
        return None

    if not str(head_filepath).endswith('.head.h5'):
        print('Error: Invalid file extension. Please input a filename that ends'
             +' with .head.h5', file=sys.stderr)
        return None

    wininfo = None
    with h5py.File(head_filepath, 'r') as f_head:
        num_win = f_head['/wininfo/nwin'][0]
        dt = np.dtype([('iwin', 'i4'), ('line_id', 'U64'),
                       ('wvl_min', 'f'), ('wvl_max','f'),
                       ('nl', 'i4'), ('xs','i4')])
        wininfo = np.recarray((num_win,), dtype=dt)
        for iwin in range(num_win):
            line_id = f_head[f'/wininfo/win{iwin:02d}/line_id'][0]
            wvl_min = f_head[f'/wininfo/win{iwin:02d}/wvl_min'][0]
            wvl_max = f_head[f'/wininfo/win{iwin:02d}/wvl_max'][0]
            win_nl = f_head[f'/wininfo/win{iwin:02d}/nl'][0]
            win_xs = f_head[f'/wininfo/win{iwin:02d}/xs'][0]
            wininfo[iwin].iwin = iwin
            wininfo[iwin].line_id = line_id.decode('utf-8')
            wininfo[iwin].wvl_min = wvl_min
            wininfo[iwin].wvl_max = wvl_max
            wininfo[iwin].nl = win_nl
            wininfo[iwin].xs = win_xs

    return wininfo
