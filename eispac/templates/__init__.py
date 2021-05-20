from .search import search
from .eis_find_templates import eis_find_templates
from .eis_read_wininfo import eis_read_wininfo
from .read_wininfo import read_wininfo

import pathlib as _pathlib # should be hidden
template_dir = str(_pathlib.Path(__file__).parent.joinpath('eis_template_dir')) 
