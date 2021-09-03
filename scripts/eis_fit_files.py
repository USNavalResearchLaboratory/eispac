#!/usr/bin/env python
__all__ = ['eis_fit_files']

import os
import sys
import warnings
import pathlib
# from astropy.utils.exceptions import AstropyWarning, AstropyDeprecationWarning
import eispac

def eis_fit_files():

    if len(sys.argv) !=3:
        print('ERROR: Invalid number of input arguments.')
        print('   Please input either lists or directories containing')
        print('   the data template files you wish to fit. For example,')
        print('   >> eis_fit_files {DATA_DIR} {TEMPLATE_DIR}')
        sys.exit()

    input_data = str(sys.argv[1])
    if isinstance(input_data, (str, pathlib.Path)):
        abs_data_path = pathlib.Path(input_data).resolve()
        if abs_data_path.is_file():
            data_list = [abs_data_path]
        elif abs_data_path.is_dir():
            data_list = list(abs_data_path.glob('*data.h5'))
        else:
            print(f'ERROR: Data file or directory does not exist: {input_data}')
            sys.exit()
    elif isinstance(input_data, list):
        data_list = input_data
    else:
        print('ERROR: Invalid type for input data file or directory')
        print('Please input either a string, list, or pathlib.Path.')
        sys.exit()

    input_template = str(sys.argv[2])
    if isinstance(input_template, (str, pathlib.Path)):
        abs_template_path = pathlib.Path(input_template).resolve()
        if abs_template_path.is_file():
            template_list = [abs_template_path]
        elif abs_template_path.is_dir():
            template_list = list(abs_template_path.glob('*template.h5'))
        else:
            print(f'ERROR: Template file or directory does not exist: {input_template}')
            sys.exit()
    elif isinstance(input_template, list):
        template_list = input_template
    else:
        print('ERROR: Invalid type for input template file or directory')
        print('Please input either a string, list, or pathlib.Path.')
        sys.exit()

    n_files = len(data_list)
    n_templates = len(template_list)
    print(f'Found {n_files} data files to be fit with {n_templates} templates.')

    # Note: we should look into the astropy warnings raised by the code
    #       For now, we are just silencing them to avoid spamming the console
    # warnings.simplefilter('ignore', AstropyWarning)
    # warnings.simplefilter('ignore', AstropyDeprecationWarning)
    if not sys.warnoptions:
        warnings.simplefilter('ignore')
        os.environ['PYTHONWARNINGS'] = 'ignore'
    for t in range(n_templates):
        print('')
        print('-----')
        print('')
        tmplt = eispac.read_template(template_list[t])

        for d in range(n_files):
            print('')
            eis_cube = eispac.read_cube(data_list[d], tmplt.central_wave)
            fit_res = eispac.fit_spectra(eis_cube, tmplt, ncpu='max',
                                         unsafe_mp=True, ignore_warnings=True)

            temp_data_path = pathlib.Path(data_list[d]).resolve()
            output_dir = temp_data_path.parent
            saved_files = eispac.save_fit(fit_res, save_dir=output_dir)
            saved_fits_files = eispac.export_fits(fit_res, save_dir=output_dir)

    warnings.resetwarnings()
    os.environ['PYTHONWARNINGS'] = 'default'

if __name__ == '__main__':
    eis_fit_files()
