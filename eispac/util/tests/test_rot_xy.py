import pytest

import numpy as np

from eispac.util import rot_xy


@pytest.mark.parametrize('cen,cen_rot_expected,start,end',[
    ([-131.43,-1.64], [-122.28026191424578,-1.7654771021022289], '2012-09-24T10:50:26.000', '2012-09-24T11:50:26.000'),
    ([0,0], [9.471882229186884, 0.08095646740947912], '2021-JAN-01 00:00', '2021-JAN-01 01:00'),
])
def test_rot_xy(cen, cen_rot_expected, start, end):
    cen_rot = rot_xy(*cen, start, end)
    cen_rot = [cen_rot.Tx.value, cen_rot.Ty.value]
    assert np.allclose(cen_rot, cen_rot_expected)
