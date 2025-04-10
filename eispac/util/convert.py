"""Useful time-related functions."""

__all__ = ['tai2utc', 'utc2tai']

import sys
from datetime import datetime
from astropy.time import Time
from dateutil.parser import parse

epoch_diff = 694656000.0

def tai2utc(tai):
    """Convert TAI (in sec from midnight 1958 Jan 1) to UTC.

    Parameters
    ----------
    tai : int or float
        Time in International Atomic Time (TAI), which is defined as the number
        of seconds after 1958-01-01 00:00 (NOT including leap seconds)

    Returns
    -------
    utc : str
        Date string of the form YYYY-MM-DDTHH:MM:SS
    """
    # Validate input
    if not isinstance(tai, (int,float)):
        print(f'ERROR: {tai} is not a valid TAI int or float!', sys.stderr)
        raise ValueError
    
    date_tai = Time(tai-epoch_diff-19.0, format='gps', scale='tai')
    d = Time(date_tai, format='isot', scale='utc')
    return d.isot

def utc2tai(utc):
    """Convert UTC to TAI (in sec from midnight 1958 Jan 1).

    Parameters
    ----------
    utc : str
        Date string of the form YYYY-MM-DDTHH:MM:SS. Other formats
        are also accepted, since the date parser is extremely flexible,
        however some odd inputs may have unexpected results (e.g. "May 4"
        will be parsed as {CURRENT YEAR}-05-04)

    Returns
    -------
    tai : float
        Time in International Atomic Time (TAI), which is defined as the number
        of seconds after 1958-01-01 00:00 (NOT including leap seconds)
    """
    # Validate input and remove extra whitespace
    if not isinstance(utc, str):
        print(f'ERROR: {utc} is not a valid time string!', sys.stderr)
        raise ValueError
    else:
        input_utc = utc.strip()
    
    # Set the default datetime used for filling in missing data
    default_datetime = datetime(datetime.now().year, 1, 1, 0, 0, 0)

    # Fix short strings not normally caught by the parser
    if any(c.isalpha() for c in input_utc):
        # Skip any string containing letters (e.g. 2015-May-04, feb29)
        pass
    elif len(input_utc) == 6 and input_utc.isdecimal():
        # Fix YYYYMM
        input_utc = input_utc + '01' # now YYYYMM01
    elif len(input_utc) <= 7 and '.' in input_utc:
        # Fix YYYY.MM or MM.YYYY
        input_utc = input_utc.replace('.', '-')

    date = parse(input_utc, default=default_datetime)
    date_string = date.strftime('%Y-%m-%dT%H:%M:%S.%f')
    date_utc = Time(date_string, format='isot', scale='utc')
    tai = date_utc.gps + epoch_diff + 19.0
    return tai

if __name__ == "__main__":
    date = ['2018-09-17T19:46:25.000', '2018-09-17 19:46:25.000',
            '2018-sep-17T19:46:25.000','2018-sep-17 19:46:25.000',
            '20180917 19:46:25.000', '17-Sep-2018 19:46:25.000',
            '17-Sep-2018 19:46:25', '2018/09/17 19:46:25.000']
    for day in date:
        tai = utc2tai(day)
        utc = tai2utc(tai)
        print(f'{day:>30} {tai:20.0f} {utc:>30}')
