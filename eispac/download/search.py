
import os
from pathlib import Path
import subprocess
from .download_db import download_db

def search():

    path = Path(__file__).parent.absolute()
    exe = os.path.join(path, 'eis_catalog.py')
    db = os.path.join(path, 'eis_cat.sqlite')

    if not os.path.isfile(db):
        print(' ! no local db file exists, downloading from the web')
        download_db()
    
    subprocess.call([exe, db])
