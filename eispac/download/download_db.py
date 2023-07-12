
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
        Defaults to the same directory as this function (e.g. ../eispac/download/).

    Returns
    -------
    local_filepath : str
        Full filepath to the downloaded database. If the download failed, for
        any reason, the filename will end with .part (any exisiting database
        will NOT be overwirtten with a partial or failed download).
    """

    if download_dir is None:
        local_dir = Path(__file__).parent.absolute()
    else:
        local_dir = str(download_dir)
    db_name = 'eis_cat.sqlite'
    local_filepath = os.path.join(local_dir, db_name)

    # if os.path.isfile(local_filepath):
    #     os.remove(local_filepath)

    if source.lower().startswith('nrl'):
        remote_filepath = 'https://eis.nrl.navy.mil/level1/db/eis_cat.sqlite'
    else:
        remote_filepath = 'https://hesperia.gsfc.nasa.gov/ssw/hinode/eis/database/catalog/eis_cat.sqlite'

    print(f' + downloading eis_cat.sqlite')
    print(f'   remote: {remote_filepath}')
    print(f'   local: {local_filepath}')
    # wget.download(remote_file, local_file)
    # print(f'   done')
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
