
import wget
import os
from pathlib import Path

def download_db(download_dir=None):

    if download_dir is None:
        local_dir = Path(__file__).parent.absolute()
    else:
        local_dir = str(download_dir)
    local_file = 'eis_cat.sqlite'
    local_file = os.path.join(local_dir, local_file)

    if os.path.isfile(local_file):
        os.remove(local_file)

    remote_file = 'https://hesperia.gsfc.nasa.gov/ssw/hinode/eis/database/catalog/eis_cat.sqlite'

    print(f' + downloading eis_cat.sqlite')
    print(f'   remote: {remote_file}')
    print(f'   local: {local_file}')
    wget.download(remote_file, local_file)
    print(f'   done')

    return local_file
