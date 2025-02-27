from ..db.download_db import download_db as new_download_db
from ..db.download_hdf5_data import download_hdf5_data as new_download_hdf5_data
from .run_eis_catalog import run_eis_catalog

from warnings import warn

def download_db(*args, **kwargs):
    warn('eispac.download.download_db() has moved! Please use'
         +' eispac.db.download_db() instead', DeprecationWarning, 2)
    return new_download_db(*args, **kwargs)
    
def download_hdf5_data(*args, **kwargs):
    warn('eispac.download.download_hdf5_data() has moved! Please use'
         +' eispac.db.download_hdf5_data() instead', DeprecationWarning, 2)
    return new_download_hdf5_data(*args, **kwargs)