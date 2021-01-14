import os
from pathlib import Path
import subprocess

def search(eis_filename=None):

    path = Path(__file__).parent.absolute()
    exe = os.path.join(path, 'eis_copy_templates.py')

    com = [exe]
    if eis_filename is not None:
        com.append(eis_filename)
    
    subprocess.call(com)
