__all__ = ['get_remote_db_modtime']

import sys
import time
import urllib

def get_remote_db_modtime(source='nrl'):
    """Fetch the modification date of the EIS as-run catalog on a remote server

    Parameters
    ----------
    source : str, optional
        String indicating which remote server to check. Choose from 'NRL' or
        'NASA'. Default is 'NRL'

    Returns
    -------
    mod_date : str
        ISO 8601 format timestamp of the remote file modification date. If
        the remote server is unavailable, then this fuction will return None. 
    """

    if source.lower().startswith('nrl'):
        remote_filepath = 'https://eis.nrl.navy.mil/level1/db/eis_cat.sqlite'
    else:
        remote_filepath = 'https://hesperia.gsfc.nasa.gov/ssw/hinode/eis/database/catalog/eis_cat.sqlite'

    try:
        with urllib.request.urlopen(remote_filepath, timeout=30) as conn:
            last_modified = conn.headers['last-modified'] # has a stupid format
            t_struct = time.strptime(last_modified, '%a, %d %b %Y %H:%M:%S %Z')
            mod_date = time.strftime('%Y-%m-%dT%H:%M:%S', t_struct)
    except Exception as ex_text:
        print('ERROR: Unable to check the remote server for the latest EIS'
             +' as-run catalog! Please check your connection and try again.')
        print(f'Remote URL: {remote_filepath}')
        print(f'Exception: {ex_text}', sys.stderr)
        mod_date = None
    
    return mod_date
