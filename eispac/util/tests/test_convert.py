import pytest
from eispac.util.convert import tai2utc, utc2tai

def test_tai2utc_invalid_input():
    with pytest.raises(ValueError):
        for invalid_tai in ["not_a_number", None, [123456], {"tai": 1234}]:
            tai2utc(invalid_tai)


def test_utc2tai_invalid_input():
    with pytest.raises(ValueError):
        for invalid_utc in [123456, None, [2018, 9, 17], 1915991204.0]:
            utc2tai(invalid_utc)

    with pytest.raises(Exception):
        utc2tai("not a valid date string")


def test_convert_round_trip_utc():
    utc_str = '2002-09-22T21:36'
    tai = utc2tai(utc_str)
    assert isinstance(tai, float)
    recovered_utc = tai2utc(tai)
    # Compare with startswith or ISO format (ignoring ms)
    assert recovered_utc.startswith(utc_str.split(".")[0])


def test_convert_round_trip_tai():
    initial_tai = 1411421792.0  # 2002-09-22T21:36 (Hinode launch date)
    utc = tai2utc(initial_tai)
    assert isinstance(utc, str)
    recovered_tai = utc2tai(utc)
    # Roundtrip should be accurate within 1 ms
    assert abs(recovered_tai - initial_tai) < 1e-3


@pytest.mark.parametrize(
    "date_format",
    [
        "2018-09-17T19:46:25.000",
        "2018-09-17 19:46:25",
        "2018/09/17 19:46:25.000",
        "2018-sep-17 19:46:25.000",
        "2018-sept-17 19:46:25.000",
        "17-Sep-2018 19:46:25.000",
        "20180917 19:46:25.000",
    ],
)
def test_utc2tai_various_string_formats(date_format):
    tai = utc2tai(date_format)
    ref_tai = utc2tai("2018-09-17T19:46:25.000")
    assert abs(tai - ref_tai) < 1e-3


def test_utc2tai_short_format_fixes():
    tai_expected = utc2tai("2018-09-01T00:00:00.000")
    # '201809' should be normalized to '20180901'
    tai_yyyymm = utc2tai("201809")
    assert abs(tai_yyyymm - tai_expected) < 1e-3

    # '2018.09' should be normalized to '2018-09'
    tai_dot = utc2tai("2018.09")
    assert abs(tai_dot - tai_expected) < 1e-3


