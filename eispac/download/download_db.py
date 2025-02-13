
import os
from pathlib import Path
import parfive

def download_db(download_dir=None, source='nrl'):
    """Download the offical EIS as-run SQLite database

    Parameters
    ----------
    download_dir : str or `pathlib.Path` object, optional
        Local directory where the database should be saved. If there is already
        an EIS as-run database in the target directory, it will be overwritten.
        Defaults to the user's home directory.
    source : str, optional
        Short name of remote server to download the catalog from. Select from 
        'NRL' or 'NASA'. Default is 'NRL', which is usually more up-to-date.

    Returns
    -------
    local_filepath : str
        Full filepath to the downloaded database. If the download failed, for
        any reason, the filename will end with .part (any exisiting database
        will NOT be overwirtten with a partial or failed download).
    """

    # Set default download dir to the user's home dir
    if download_dir is None:
        local_dir = str(Path.home()) # Path(__file__).parent.absolute()
    else:
        local_dir = str(download_dir)

    # Ensure local dir is a string path to a valid dir (and not a file)
    if os.path.isfile(local_dir):
        local_dir = os.path.dirname(local_dir)
    elif not os.path.isdir(local_dir):
        raise OSError(f'{local_dir} is not a valid dir! Unable to download EIS as-run catalog.')
    
    db_name = 'eis_cat.sqlite'
    local_filepath = os.path.join(local_dir, db_name)

    if source.lower().startswith('nrl'):
        remote_filepath = 'https://eis.nrl.navy.mil/level1/db/eis_cat.sqlite'
    else:
        remote_filepath = 'https://hesperia.gsfc.nasa.gov/ssw/hinode/eis/database/catalog/eis_cat.sqlite'

    print(f' + downloading eis_cat.sqlite')
    print(f'   remote: {remote_filepath}')
    print(f'   local: {local_filepath}')
    dl = parfive.Downloader(max_conn=2, overwrite=True)
    dl.enqueue_file(remote_filepath, path=local_dir,
                    filename=db_name+'.part', verify_ssl=False)
    dl_result = dl.download()

    # If there were no errors, delete the old db (if any) and rename the new
    # NB: there is an edge case where a file fully downloads (after one or more
    #     failed tries), but is not properly closed by parfive. Hence the "try"
    if len(dl_result.errors) == 0:
        try:
            os.replace(local_filepath+'.part', local_filepath)
        except Exception as err:
            print(err)
    else:
        print(dl_result.errors[0][2])
        print(f'Incomplete file download. Please check your connection and try again.')
        if not os.path.isfile(local_filepath):
            local_filepath = local_filepath+'.part'

    return local_filepath
