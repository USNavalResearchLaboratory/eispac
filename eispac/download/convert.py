"""Useful time-related functions."""

import datetime
from astropy.time import Time
from dateutil.parser import parse

epoch_diff = 694656000.0

def tai2utc(tai):
    """Convert TAI (in sec from midnight 1958 Jan 1) to UTC."""
    date_tai = Time(tai-epoch_diff-19.0, format='gps', scale='tai')
    d = Time(date_tai, format='isot', scale='utc')
    return d.isot

def utc2tai(utc):
    """Convert UTC to TAI (in sec from midnight 1958 Jan 1).

    UTC is a string of the form YYYY-MMM-DDTHH:MM:SS. Other forms
    should also work as the date parser is pretty good.
    """
    date = parse(utc)
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
