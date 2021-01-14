__all__ = ['EISFitTemplate', 'create_funcinfo', 'read_template']

import os
import sys
import pathlib
import numpy as np
import h5py

class EISFitTemplate:

    def __init__(self):
        self.filename_temp = None
        self.central_wave = None
        self.template = None
        self.parinfo = None
        self.funcinfo = None

    def print_parinfo(self):
        if self.parinfo is None:
            print('parinfo is currently empty! Please load a template file.')
        else:
            print('--- FIT TEMPLATE PARAMETER CONSTRAINTS ---')
            print(f"{'*':>4} {'Value':>16} {'Fixed':>10} "
                 +f"{'Limited':>18} {'Limits':>22} {'Tied':>18}")
            for i in range(len(self.parinfo)):
                print(f"{'p['+str(i)+']':>6} "
                     +f"{self.parinfo[i]['value']:14.4f} "
                     +f"{self.parinfo[i]['fixed']:10d} "
                     +f"{self.parinfo[i]['limited'][0]:10d} "
                     +f"{self.parinfo[i]['limited'][1]:10d} "
                     +f"{self.parinfo[i]['limits'][0]:12.4f} "
                     +f"{self.parinfo[i]['limits'][1]:12.4f} "
                     +f"{self.parinfo[i]['tied']:>18}")

def create_funcinfo(template):
    funcinfo = []
    for g in range(template['n_gauss']):
        funcinfo.append({'func':'Gaussian1D',
                         'name':template['line_ids'][g],
                         'n_params':3})
    if template['n_poly'] > 0:
        funcinfo.append({'func':'Polynomial1D',
                         'name':'Background',
                         'n_params':template['n_poly']})
    return funcinfo

def read_template(filename=None, quiet=False):

    # NOTE: added 'central_wave' value
    # Quick input validation (value checks are implemented later)
    if not isinstance(filename, (str, pathlib.Path)):
        print('Error: Please input a valid template filepath as '
             +'either a string or pathlib.Path object', file=sys.stderr)
        return None

    # Initialize output object
    output = EISFitTemplate()

    # Parse filename and determine the directory and filename
    abs_filepath = pathlib.Path(filename).resolve()
    input_name = str(abs_filepath.name)
    input_dir = abs_filepath.parent
    if str(input_dir) == '.':
        input_dir = pathlib.Path().cwd()

    tmplt_filepath = input_dir.joinpath(input_name)

    # Check for the template file. Exit if it does not exist.
    if not tmplt_filepath.is_file():
        print('Error: Template file does not exist, ' + str(tmplt_filepath), file=sys.stderr)
        return None
    else:
        output.filename_temp = str(tmplt_filepath)
        print('Template file,\n   ' + str(tmplt_filepath))

    # read template data
    with h5py.File(tmplt_filepath, 'r') as f_temp:
        template = {}
        for key in f_temp['template']:
            val = f_temp['template/'+key]
            if key == 'line_ids':
                val = np.char.decode(val).flatten() # convert bytes to unicode
            elif len(val) > 1:
                val = np.array(val)
            else:
                val = val[0]
            template[key] = val

        output.template = template

        # read parinfo
        nstr = len(f_temp['parinfo/value'])
        parinfo = []
        for istr in range(nstr):
            parameter = {}
            for key in f_temp['parinfo']:
                val = f_temp['parinfo/'+key][istr]
                if key == 'tied':
                    val = np.char.decode(val) # convert bytes to unicode
                parameter[key] = val
            parinfo.append(parameter)

        output.parinfo = parinfo

    output.funcinfo = create_funcinfo(output.template)

    # Calculate central wavelength (saves the user a little time and energy)
    output.central_wave = template['wmin'] + (template['wmax']-template['wmin'])*0.5

    if not quiet:
        output.print_parinfo()

    return output

if __name__ == '__main__':

    import pathlib

    filename = './templates/eis_template_dir/fe_12_195_119.2c.template.h5'
    filename = str(pathlib.Path(filename).resolve())
    Fe_XII_195_119 = read_template(filename)

    print('\ntemplate')
    print('---------')
    t = Fe_XII_195_119.template
    for key in t.keys():
        print('{:12} {:12} {:4d}'.format(key,str(t[key].dtype),np.size(t[key])))

    print('\nparinfo['+str(len(Fe_XII_195_119.parinfo))+']')
    print('---------')
    t = Fe_XII_195_119.parinfo[0]
    for key in t.keys():
        print('{:12} {:12} {:4d}'.format(key,str(t[key].dtype),np.size(t[key])))
