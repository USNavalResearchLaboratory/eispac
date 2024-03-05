__all__ = ['EISFitTemplate']

import sys
import copy
import pathlib
import warnings

import numpy as np
import h5py

if sys.version_info >= (3, 11):
    import tomllib
else:
    import tomli as tomllib

class EISFitTemplate:
    """Multigaussian fitting template for use with MPFIT and `~eispac.core.fit_spectra`

    Parameters
    ----------
    filename : str or `pathlib.Path`, optional
        Path to the template file. Default is "unknown_file". Note: passing a
        filename directly to EISFitTemplate will NOT automatically load the file;
        please use `~eispac.core.read_template` instead.
    template : dict, optional
        Dictionary of template parameters. Valid keys include:

        * ``n_gauss`` (int) - number of Gaussian components
        * ``n_poly`` (int) - Number of background polynomial terms. 
            Common values are: 0 (no background), 1 (constant), and 2 (linear).
        * ``line_ids`` (array_like) - Strings giving the line identification 
            for each Gaussian component. For example, "Fe XII 195.119". 
            If not specified, placeholder values of 
            "unknown I {INITAL CENTROID VALUE}" will be used.
        * ``wmin`` (float) - min wavelength value of data to use for fitting 
        * ``wmax`` (float) - max wavelength value of data to use for fitting 

    parinfo : list or dict, optional
        Either a list of dicts or a dict of lists giving the initial fitting 
        parameters and constraints. If given a `list`, each entry must be a dict
        with the correct keys for a single fit parameter. If given a `dict`, each
        key must contain an array of list for all parameters. The order of 
        parameters is assumed to be sets of [PEAK, CENTROID, WIDTH] for each 
        Gaussian component followed by any coefficients for the background 
        polynomial, starting with the LOWEST (constant) order term first. 
        Valid keys include:

        * ``value`` (float) - initial parameter guess
        * ``fixed`` (0 or 1) - If set to "1", will not fit and just use initial value
        * ``limited`` (two-element array_like) - If set to "1" in the first/second
            value, will apply and limit to the parameter on the lower/upper side
        * ``limits`` (two-element array_like) - Values of the limits on the
            lower/upper side. Both "limited" and "limits" must be give together.
        * ``tied`` (str) - String defining a fixed relationship between this
            parameters one or more others. For example "p[0] / 2" would define
            a parameter to ALWAYS be exactly half the value of the first parameter.

        Additional keys are available, please see the MPFIT documentation in the
        EISPAC user's guide for more details.
    kwargs : value or array_like, optional
        For convenience, users can also pass in any "template" or "parinfo" value
        or array as a seperate keyword. Any values defined in this way will
        take priority and overwrite any previous values for that key.

    Attributes
    ----------
    template : dict
        Dictoary of template paramaters
    parinfo : list of dicts
        List of parameter constraint dicts
    """

    def __init__(self, filename=None, template=None, parinfo=None, **kwargs):
        if filename is None:
            self.filename_temp = 'unknown_file'
        else:
            self.filename_temp = str(filename)
        
        self.template = dict()
        if isinstance(template, dict):
            self.template = copy.deepcopy(template)
            # Ensure that all keys are lower case
            self.template = {KEY.lower():VAL for KEY, VAL in self.template.items()}
        elif template is not None:
            print('Error: Invalid datatype for "template". Please input a '
                 +'dictionary. Initializing an empty dict using the other input '
                 +'parameters.', file=sys.stderr)
            
        self.parinfo = dict()
        if isinstance(parinfo, (list, tuple, dict)):
            self.parinfo = copy.deepcopy(parinfo)
        elif parinfo is not None:
            print('Error: Invalid datatype for "parinfo". Please input a '
                 +'list of dicts or a dict of lists with the correct keys. '
                 +'Initializing an empty list of dicts using the other input '
                 +'parameters.', file=sys.stderr)
            
        # Copy over kwargs (if any) to the correct dict or list
        if len(kwargs) > 0:
            # Ensure that all keys are lower case
            kwargs = {KEY.lower():VAL for KEY, VAL in kwargs.items()}
            for KEY in kwargs.keys():
                if KEY in ['component', 'n_gauss', 'n_poly', 'line_ids',
                           'wmin', 'wmax', 'data_e', 'data_x', 'data_y',
                           'fit', 'fit_back', 'fit_gauss', 'order']:
                    # .template keys
                    self.template[KEY] = copy.deepcopy(kwargs[KEY])
                elif KEY in ['value', 'fixed', 'limited', 'limits', 
                             'tied', 'parname', 'step', 
                             'mpmaxstep', 'mpside', 'mpprint']: 
                    # .parinfo keys
                    if isinstance(self.parinfo, dict):
                        self.parinfo[KEY] = copy.deepcopy(kwargs[KEY])
                    elif isinstance(self.parinfo, list):
                        n_new = len(kwargs[KEY])
                        n_par = len(self.parinfo)
                        for p in range(min(n_new, n_par)):
                            self.parinfo[p][KEY] = copy.deepcopy(kwargs[KEY][p])
                else: 
                    print(f'Warning: "{KEY}" is not a valid template or parinfo '
                          +f'keyword and will be ignored. Please see the docs '
                          +f'for a list of supported inputs.', file=sys.stderr)

        # Initialize values of n_gauss and n_poly (if not given)
        if ('n_gauss' not in self.template.keys()) or ('n_poly' not in self.template.keys()):
            if len(self.parinfo) > 0:
                # If parinfo is defined, guess from number of values
                if isinstance(self.parinfo, (list, tuple)):
                    # Standard list of dicts
                    n_params = len(self.parinfo)
                elif isinstance(self.parinfo, dict):
                    # Alternative dict of arrays
                    n_params = 1
                    for KEY in self.parinfo.keys():
                        n_params = max(n_params, len(self.parinfo[KEY]))
                # Add the missing value to the template
                if 'n_gauss' not in self.template.keys():
                    self.template['n_gauss'] = int(n_params / 3)
                if 'n_poly' not in self.template.keys():
                    self.template['n_poly'] = int(n_params % 3)
            else:
                # If neither template nor parinfo are defined, use defaults
                if 'n_gauss' not in self.template.keys():
                    self.template['n_gauss'] = 1
                if 'n_poly' not in self.template.keys():
                    self.template['n_poly'] = 1

        # If parinfo is not given, initialize minimal parinfo list of dicts
        if len(self.parinfo) == 0:
            n_gauss = max(1, self.template.get('n_gauss', 1))
            n_poly = max(0, self.template.get('n_poly', 1))

            par_list = []
            for g in range(3*n_gauss):
                par_list.append({'value': 1.0,
                                 'limited': np.array([1, 0], dtype='int16')})
            for p in range(n_poly):
                par_list.append({'value': 1.0,
                                 'limited': np.array([0, 0], dtype='int16')})
            self.parinfo = par_list

        # Finally, validate the actual values and fix format issues
        # Note: we need to validate parinfo first, so we can convert the format
        #       and more easily check the template values of n_gauss and n_poly
        self._validate_parinfo_list()
        self._validate_template_dict()

    def _validate_parinfo_list(self):
        """Helper function for validating the length and keys of .parinfo
        """

        # Default dict for a single parameter
        default_parinfo = {'fixed': 0,
                           'limited': np.array([1, 0], dtype='int16'),
                           'limits': np.zeros(2),
                           'tied': np.array(' ', dtype='<U16'),
                           'value': 1.0}
        
        # If given a dict of arrays (or lists), convert to a list of dicts
        if isinstance(self.parinfo, dict):
            # First, ensure all keys are lower case
            input_parinfo = {KEY.lower():VAL for KEY, VAL in self.parinfo.items()}
            if 'value' in input_parinfo.keys():
                n_params = len(self.parinfo['value'])
            else:
                max_len = 0
                for KEY in input_parinfo.keys():
                    this_len = len(input_parinfo[KEY])
                    max_len = max(this_len, max_len)
                n_params = max_len
                print('Warning: No initial values given in parinfo! Please '
                    +'define a "value" key for each parameter', file=sys.stderr)

            # Extract the values and assemble the correct list of dicts
            par_list = []
            for p in range(n_params):
                info_dict = {}
                for KEY in input_parinfo.keys():
                    try:
                        if KEY in ['tied', 'parname']: # strings
                            info_dict[KEY] = np.array(str(input_parinfo[KEY][p]), 
                                                      dtype='<U16')
                        elif KEY in ['limited', 'limits']: # two-element
                            info_dict[KEY] = copy.deepcopy(default_parinfo[KEY])
                            info_dict[KEY][0] = input_parinfo[KEY][p][0]
                            info_dict[KEY][1] = input_parinfo[KEY][p][1]
                        elif KEY in ['fixed', 'mpside', 'mpprint']: # integers
                            info_dict[KEY] = int(input_parinfo[KEY][p])
                        elif KEY in ['value', 'step', 'mpmaxstep']: # floats
                            info_dict[KEY] = float(input_parinfo[KEY][p])
                        else:
                            print(f'Warning: {KEY} is not a valid parinfo key '
                                 +'and will be ignored.', file=sys.stderr)
                    except:
                        # If there is an issue, give a warning and continue
                        print(f'Warning: There was an issue loading the {KEY} '
                             +f'key for parameter index {p}. Skipping for now; '
                             +f'please check the format and datatype.', 
                             file=sys.stderr)
                        continue
                par_list.append(info_dict)
            
            self.parinfo = par_list
        
        # Check for missing keys expected by eispac and fill with default values
        print_limit_warning = False
        n_params = len(self.parinfo)
        for p in range(n_params):
            # First, ensure all keys are lower case
            self.parinfo[p] = {KEY.lower():VAL for KEY, VAL in self.parinfo[p].items()}
            for KEY, VALUE in default_parinfo.items():
                # Add missing keys that are required for printing details
                if KEY not in self.parinfo[p].keys():
                    self.parinfo[p][KEY] = copy.deepcopy(VALUE)
                else: 
                    # Check shape of two-element keys that DO exist
                    if KEY in ['limited', 'limits']:
                        if len(self.parinfo[p][KEY]) != 2:
                            # Overwrite if given wrong shape
                            self.parinfo[p][KEY] = copy.deepcopy(VALUE)
                            print_limit_warning = True

        if print_limit_warning:
            print(f'Warning: incorrect length of "limited" and/or "limits" '
                  +f'in found parinfo. Both "limited" and "limits" should be '
                  +f'two-element lists or arrays. Invalid inputs have been '
                  +f'replaced with default values.', file=sys.stderr)
    
    def _validate_template_dict(self):
        """Helper function for validating the keys and values of .template
        """
        # Ensure that all keys are lower case
        self.template = {KEY.lower():VAL for KEY, VAL in self.template.items()}

        # Validate keys actaully used by eispac
        n_gauss = self.template.get('n_gauss', -1)
        if isinstance(n_gauss, (int, float)) and int(n_gauss) > 0:
            self.template['n_gauss'] = int(n_gauss)
        else:
            self.template['n_gauss'] = -1
            print('Error: Invalid value for n_gauss. '
                 +'Please input a non-zero integer', file=sys.stderr)

        n_poly = self.template.get('n_poly', -1)
        if isinstance(n_poly, (int, float)) and int(n_poly) >= 0:
            self.template['n_poly'] = int(n_poly)
        else:
            self.template['n_poly'] = -1
            print('Error: Invalid value for n_poly. '
                  +'Please input an integer >= 0.', file=sys.stderr)
        
        wmin = self.template.get('wmin', 170.0)
        if isinstance(wmin, (int, float)):
            self.template['wmin'] = float(wmin)
        else:
            self.template['wmin'] = 170.0
            print('Error: Invalid datatype for wmin. Please input a float. '
                 +'Using a default value of 170.', file=sys.stderr)

        wmax = self.template.get('wmax', 292.0)
        if isinstance(wmax, (int, float)):
            self.template['wmax'] = float(wmax)
        else:
            self.template['wmax'] = 292.0
            print('Error: Invalid datatype for wmax. Please input a float. '
                 +'Using a default value of 292.', file=sys.stderr)

        comp_num = self.template.get('component', 1)
        if isinstance(comp_num, (int, float)):
            self.template['component'] = int(comp_num)
        else:
            self.template['component'] = 1
            print('Error: Invalid datatype for component. Please input an '
                 +'integer. Using a default value of 1.', file=sys.stderr)
                
        # Check number of parameters and components
        n_params = len(self.parinfo)
        if n_params != (3*self.template['n_gauss'] + self.template['n_poly']):
            # When mismatched, assume input parinfo is correct
            # Then, update/overwrite parts of the template as needed
            n_gauss = int(n_params / 3)
            n_poly = int(n_params % 3)

            self.template['n_gauss'] = n_gauss
            self.template['n_poly'] = n_poly
            self.template['line_ids'] = np.zeros(n_gauss, dtype='<U32')
            for g in range(n_gauss):
                self.template['line_ids'][g] = 'unknown I 001.000'

            print(f'Warning: the values of n_gauss and n_poly do not match the '
                 +f'length of parinfo. Defaulting to the best guess values of '
                 +f'n_gauss = {n_gauss} and n_poly = {n_poly}. Template dict '
                 +f' values and line_ids have been overwritten to match.', 
                  file=sys.stderr)

        # Check for line_id and convert as needed
        n_gauss = self.template['n_gauss']
        n_poly = self.template['n_poly']
        if 'line_ids' not in self.template.keys():
            self.template['line_ids'] = np.zeros(n_gauss, dtype='<U32')
            for g in range(n_gauss):
                self.template['line_ids'][g] = 'unknown I 001.000'
        elif isinstance(self.template['line_ids'], (str, int, float)):
            # Convert bare values to an array (actual values checked later)
            bare_id = self.template['line_ids']
            self.template['line_ids'] = np.array([bare_id], dtype='<U24')
        else: 
            # Force everything else to be an array of strings
            current_ids = self.template['line_ids']
            self.template['line_ids'] = np.array(current_ids, dtype='<U24')
        
        # Trim or expand line_ids array to have the correct length
        if len(self.template['line_ids']) != n_gauss:
            old_id_arr = self.template['line_ids']
            self.template['line_ids'] = np.zeros(n_gauss, dtype='<U32')
            for g in range(n_gauss):
                if g < len(old_id_arr):
                    self.template['line_ids'][g] = old_id_arr[g]
                else:
                    self.template['line_ids'][g] = 'unknown I 001.000'

        # Validate format and values of line_id strings
        for g in range(n_gauss):
            parinfo_wave_str = f"{self.parinfo[3*g+1]['value']:07.3f}"

            this_id = self.template['line_ids'][g].lower()
            if this_id.startswith(('unknown', 'no line')) or len(this_id) < 1:
                # If a line is unknown, just update the wavelength substring
                self.template['line_ids'][g] = f"unknown I {parinfo_wave_str}"
            else:
                # Parse the ID string and check values
                warn_bad_id = False
                input_id = self.template['line_ids'][g]
                split_id = input_id.split()

                # Element - can be any string, really
                elem_id = split_id[0]

                # Ionization state - only allow roman numerals less than 100
                ion_id = 'I'
                if len(split_id) >= 2:
                    test_ion = split_id[1].upper()
                    num_roman = sum(map(test_ion.count, ['I','V','X','L']))
                    if len(test_ion) == num_roman:
                        ion_id = test_ion
                    else:
                        warn_bad_id = True
                else:
                    warn_bad_id = True

                # Wavelength - needs to have a decimal point
                wave_id = parinfo_wave_str
                if len(split_id) >= 3:
                    test_wave = split_id[2]
                    if test_wave.isdigit():
                        wave_id = test_wave+'.000'
                    elif '.' in test_wave:
                        wave_id = test_wave
                    else:
                        warn_bad_id = True
                else:
                    warn_bad_id = True

                # Copy over the clean ID (all other string parts are ignored)
                clean_id = ' '.join([elem_id, ion_id, wave_id])
                self.template['line_ids'][g] = clean_id

                if warn_bad_id:
                    print(f'Warning: "{input_id}" is not a valid or complete '
                         +f'line ID string. Using a replacement ID of '
                         +f'"{clean_id}" instead.', file=sys.stderr)

        # Extract initial values and update the 'fit' array
        # Note: this is important for using the scale_guess() function during
        #       the fitting process. 
        self.template['fit'] = np.zeros(3*n_gauss + n_poly)
        for p in range(3*n_gauss + n_poly):
            self.template['fit'][p] = self.parinfo[p]['value']

    def __repr__(self):
        rows = []
        rows.append('--- EISFitTemplate SUMMARY ---')
        rows.append(f"filename_temp: {self.filename_temp}")
        rows.append(f"n_gauss: {self.template['n_gauss']}")
        rows.append(f"n_poly: {self.template['n_poly']}")
        rows.append(f"line_ids: {self.template['line_ids']}")
        rows.append(f"wmin, wmax: {self.template['wmin']}, {self.template['wmax']}")
        rows.append('')
        rows.append('--- PARAMETER CONSTRAINTS ---')
        rows.append(f"{'*':>4} {'Value':>16} {'Fixed':>7} "
                   +f"{'Limited':>11} {'Limits':>19} {'Tied':>19}")
        for i, p in enumerate(self.parinfo):
            rows.append(f"{f'p[{i}]':>6} {p['value']:14.4f} {p['fixed']:7d} "
                       +f"{p['limited'][0]:5d} {p['limited'][1]:5d} "
                       +f"{p['limits'][0]:12.4f} {p['limits'][1]:12.4f} {p['tied']:>18}")
        return '\n'.join(rows)

    @property
    def central_wave(self):
        """Wavelength value in the center of the template wavelength range
        
        Note: this is calculated using the current values of 'wmin' and 'wmax'
        contained in the template dict.
        """
        return self.template['wmin'] + (self.template['wmax'] - self.template['wmin']) * 0.5
    
    @property
    def funcinfo(self):
        """List of dicts specifying each subcomponent function used in the template
        
        Note: this is generated on-demand using the current template dict and
        the ``get_funcinfo`` class method.
        """
        return self.get_funcinfo(self.template)

    @staticmethod
    def get_funcinfo(template):
        """Generate the ``funcinfo`` dict

        Returns a list of dictionaries where each entry describes the
        parameters for one of the fitting basis functions

        Parameters
        ----------
        template : dict

        Returns
        -------
        funcinfo : list
        """
        funcinfo = []
        for g in range(template['n_gauss']):
            funcinfo.append({'func': 'Gaussian1D',
                             'name': template['line_ids'][g],
                             'n_params': 3})
        if template['n_poly'] > 0:
            funcinfo.append({'func': 'Polynomial1D',
                             'name': 'Background',
                             'n_params': template['n_poly']})
        return funcinfo

    @classmethod
    def read_template(cls, filename):
        """Load an `EISFitTemplate` from an HDF5 or TOML template file

        Parameters
        ----------
        filename : str or `pathlib.Path`
            Path to a HDF5 template file provided with eispac or a 
            user-made TOML template file.

        Returns
        -------
        cls : `EISFitTemplate` class instance
            Object containing the fit template
        """
        # NOTE: return None here rather than allow h5py to handle
        # exception so that spectral fitting pipeline can error
        # more gracefully
        # FIXME: replace with proper exception handling and logging
        if not isinstance(filename, (str, pathlib.Path)):
            warnings.warn('Error: Template filepath must be either a string or pathlib.Path')
            return None
        
        filename = pathlib.Path(filename)
        if not filename.is_file():
            warnings.warn(f'Error: Template filepath {filename} does not exist')
            return None

        file_type = filename.suffix

        if file_type.lower() in ['.h5', '.hdf5']:
            # Load a standard format HDF5 file (probably packaged with eispac)
            with h5py.File(filename, 'r') as f_temp:
                # Template
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
                # Parinfo
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

                # Fix the datatypes of important keys
                template['n_gauss'] = int(template.get('n_gauss', 1))
                template['n_poly'] = int(template.get('n_poly', 1))
                template['wmin'] = float(template.get('wmin', 100))
                template['wmax'] = float(template.get('wmax', 1000))
                template['component'] = int(template.get('component', 1))
                template['line_ids'] = template['line_ids'].astype('<U24')
        elif file_type.lower() == '.toml':
            # Load a custom template stored in a TOML file
            with open(filename, 'rb') as f_temp:
                toml_dict = tomllib.load(f_temp)
                # Ensure top-level keys are all lower case
                toml_dict = {KEY.lower(): VALUE for KEY, VALUE in toml_dict.items()}

                template = toml_dict.get('template', None)
                parinfo = toml_dict.get('parinfo', None)
                # note: this parinfo will be an alternative dict of lists

        return cls(filename, template, parinfo)
