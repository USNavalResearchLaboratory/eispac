
import os
from pathlib import Path
import parfive

def download_db(download_dir=None):

    if download_dir is None:
        local_dir = Path(__file__).parent.absolute()
    else:
        local_dir = str(download_dir)
    db_name = 'eis_cat.sqlite'
    local_filepath = os.path.join(local_dir, db_name)

    # if os.path.isfile(local_filepath):
    #     os.remove(local_filepath)

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
