__all__ = ['run_eis_browse_templates']

import sys
import pathlib
import subprocess

def run_eis_browse_templates(filepath=None):
    """Launch the "eis_browse_templates" GUI tool""

    Parameters
    ----------
    filepath : str or pathlib.Path object, optional
        Filepath to an EIS level-1 HDF5 data or header file. If None, you may
        broswe for and select the observation file from within the GUI. Default
        is None.

    Returns
    -------
    None
    """

    if filepath is not None:
        if not isinstance(filepath, (str, pathlib.Path)):
            print('Error: Please input a valid database dir as '
                 +'either a string or pathlib.Path object', file=sys.stderr)
            return None

        # Parse dir and determine the directory and filename (if any)
        abs_filepath = pathlib.Path(filepath).resolve()
        input_name = str(abs_filepath.name)
        input_dir = abs_filepath.parent
        if str(input_dir) == '.':
            input_dir = pathlib.Path().cwd()

        obs_filepath = input_dir.joinpath(input_name)

        subprocess.call(['eis_browse_templates', obs_filepath])
    else:
        subprocess.call('eis_browse_templates')
