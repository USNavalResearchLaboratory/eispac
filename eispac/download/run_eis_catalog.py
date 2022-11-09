__all__ = ['run_eis_catalog']

import sys
import pathlib
import subprocess

def run_eis_catalog(db_dir=None):
    """Launch the "eis_catalog" GUI tool""

    Parameters
    ----------
    db_dir : str or pathlib.Path object, optional
        Directory containing the "eis_cat.sqlite" file. If None, EISPAC will
        for the catalog in an existing SSW installation or ask if you want to
        download it. Default is None. Note: any filename included with the path
        will be ignored (please do not rename your local catalog).

    Returns
    -------
    None
    """

    if db_dir is not None:
        if not isinstance(db_dir, (str, pathlib.Path)):
            print('Error: Please input a valid database dir as '
                 +'either a string or pathlib.Path object', file=sys.stderr)
            return None

        # Parse dir and determine the directory and filename (if any)
        abs_filepath = pathlib.Path(db_dir).resolve()
        input_dir = abs_filepath.parent
        if str(input_dir) == '.':
            input_dir = pathlib.Path().cwd()

        if not input_dir.is_dir():
            print("Error: db_dir is invalid or missing!", file=sys.stderr)
            return None

        db_filepath = input_dir.joinpath('eis_cat.sqlite')

        subprocess.call(['eis_catalog', db_filepath])
    else:
        subprocess.call('eis_catalog')
