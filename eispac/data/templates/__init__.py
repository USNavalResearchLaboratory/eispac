"""
Template files for fitting EIS spectra
"""
import os
import pathlib

from astropy.utils.data import get_pkg_data_filename

import eispac

__all__ = ['get_fit_template_filepath', 'fit_template_filenames']


def get_fit_template_filepath(filename, **kwargs):
    """
    Return the full path to a template file in the ``templates/eis_template_dir`` directory.

    Parameters
    ----------
    filename : `str`
        The name of the file inside the ``templates/eis_template_dir`` directory.

    Return
    ------
    filepath : `str`
        The full path to the file.

    Notes
    -----
    This is a wrapper around `astropy.utils.data.get_pkg_data_filename` which
    sets the ``package`` kwarg to be 'eispac.templates.eis_template_dir'.
    """
    return get_pkg_data_filename(filename, package="eispac.data.templates", **kwargs)


def fit_template_filenames():
    """
    Return a list of all available fitting template files.
    """
    rootdir = pathlib.Path(os.path.dirname(eispac.__file__)) / "data" / "templates"
    return sorted(list(rootdir.glob('*.template.h5')))
