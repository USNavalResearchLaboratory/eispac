#!/usr/bin/env python

import os
import numpy as np
import h5py

def eis_read_wininfo(filename_head):

    if not os.path.isfile(filename_head):
        print(' ! file does not exist ')
        return None

    if not filename_head.endswith('.head.h5'):
        print(' ! file extension not understood ')
        return None        

    wininfo = None
    with h5py.File(filename_head, 'r') as f_head:
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
            wininfo[iwin].nl = win_nl # NEW
            wininfo[iwin].xs = win_xs # NEW
            
    return wininfo

