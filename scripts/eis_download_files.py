#!/usr/bin/env python

"""
eis_download_files.py - a command line script for downloading eis files from eis.nrl.navy.mil

Usage:

* Input a list of eis filenames. Only the stem is important.

  > eis_downdload_files.py eis_l0_20201124_071409.fits eis_l0_20201124_072858

* Input a filename that contains a list of files

  > eis_download_files.py eis_files.txt

* Use '-datetree' to write the outfiles to a date tree, otherwise all files are written to
  'data_eis/' in the current directory.

"""
__all__ = ['eis_download_files']

import sys
from pathlib import Path
import eispac

def parse_filename(input_filename, extension=False):
    """
    Parse an input filename into name + '.' + ext, return these as strings
    """
    f = Path(input_filename)
    out = f.name.split('.')
    if len(out) == 1:
        name = out[0]
        ext = ''
    else:
        name = out[0]
        ext = '.'.join(out[1:])
    if extension:
        return name, ext
    else:
        return name

def parse_file(input_file):
    """
    Parse a text file containing a list of files, and perhaps some comments, return a list
    """
    out = []
    with open(input_file, 'r') as f:
        files = f.readlines()
        out = []
        for f in files:
            f = f.split('#')[0] # remove comments
            f = f.strip('  \t\n\r') # remove tabs and newlines
            if f == '': continue
            out.append(parse_filename(f))
    return out

def eis_download_files():
    # nothing is input, exit
    if len(sys.argv) == 1:
        print(' ! enter some filenames . . .')
        exit()

    # everything left is input
    arglist = sys.argv[1:]

    # look for datetree option
    datetree = False
    if '-datetree' in arglist:
        datetree = True
        arglist.remove('-datetree')

    # look for remote source options
    data_source = 'nrl'
    if '-nasa' in arglist:
        data_source = 'nasa'
        arglist.remove('-nasa')
    if '-mssl' in arglist:
        data_source = 'mssl'
        arglist.remove('-mssl')
    if '-nrl' in arglist:
        data_source = 'nrl'
        arglist.remove('-nrl')

    # loop over inputs
    for file in arglist:
        # parse the input
        name, ext = parse_filename(file, extension=True)

        if ext == 'txt':
            # a file with a list of files in it
            names = parse_file(file)
        else:
            # a scalar filename
            names = [name]

        # download the files
        for name in names:
            print(f'+ processing {name}')
            o = eispac.download.download_hdf5_data(name, source=data_source,
                                                   datetree=datetree)

if __name__ == '__main__':
    eis_download_files()
