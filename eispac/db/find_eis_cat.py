__all__ = ['find_eis_cat']

import os
import sys
import pathlib

def find_eis_cat(dir=None, quiet=False):
    """Search for EIS as-run database (eis_cat.sqlite)

    Attempts to find the EIS as-run database by searching probable locations.

    Search priority:
    1) User input dir
    2) SolarSoft (SSW) IDL installation, if available and EIS is configured
    3) User home dir
    4) Current working dir 

    Parameters
    ----------
    dir : str or `pathlib.Path` object, optional
        Directory to search FIRST. If the the catalog is not found, will then
        begin searching using environment variables and standard SSW install
        locations. Default is "None"
    quiet : bool, optional
        If set to "True", will not print status updates to the command line.
        Default is "False"


    Returns
    -------
    cat_path : `pathlib.Path` object
        Path to the eis_cat.sqlite database
    """

    # Validate user input dir
    if dir is None:
        input_dir = None
    elif not isinstance(dir, (str, pathlib.Path)):
        input_dir = None
        if not quiet:
            print('WARNING: Input dir is not a valid string or pathlib.Path'
                +' object. Reverting to default search pattern.', 
                file=sys.stderr)
    else:
        # Ensure pathlib.Path points to a valid dir (not a file)
        input_dir = pathlib.Path(dir).resolve()
        if input_dir.is_file():
            input_dir = input_dir.parent
        elif not input_dir.is_dir():
            input_dir = None
            if not quiet:
                print('WARNING: Input dir not found! Reverting to default'
                    +' search pattern.', file=sys.stderr)

    # Assemble list of search locations
    test_path_list = []

    # (1) Check user input dir
    if input_dir is not None:
        test_path_list.append(input_dir / 'eis_cat.sqlite')

    # (2) Check for SSW instalation
    ssw_dir = os.environ.get('SSW') # should exist if SSW is fully configured
    if ssw_dir is None:
        # Search for SSW in common installation directories
        for test_ssw_dir in ['/usr/local/ssw', os.path.expanduser('~')+'/ssw',
                         'C:\\ssw', 'D:\\ssw']:
            if os.path.isdir(test_ssw_dir):
                ssw_dir = test_ssw_dir
                break
    if ssw_dir is not None:
        test_path_list.append(pathlib.Path(ssw_dir) / 'hinode' / 'eis' / 
                              'database' / 'catalog' / 'eis_cat.sqlite')

    # (3) Check in user's home dir
    test_path_list.append(pathlib.Path.home() / 'eis_cat.sqlite')

    # (4) Check in current working directory
    test_path_list.append(pathlib.Path.cwd() / 'eis_cat.sqlite')

    # Finally, search all paths and return the first eis_cat found
    for PATH in test_path_list:
        if PATH.is_file():
            if not quiet:
                print(f'Found eis_cat.sqlite database!')
                print(f'   filepath: {PATH}')
            return PATH
    
    # If no catalog found, return none
    if not quiet:
        print(f'WARNING: No eis_cat.sqlite database found!', file=sys.stderr)
    return None