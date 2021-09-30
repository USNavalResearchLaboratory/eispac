from .search import search
from .eis_find_templates import eis_find_templates

import pathlib as _pathlib # should be hidden
template_dir = str(_pathlib.Path(__file__).parent.joinpath('eis_template_dir'))
