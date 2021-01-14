__all__ = ['generate_astropy_model']

from astropy.modeling import models
from eispac.core.read_template import EISFitTemplate, read_template

# if __name__ == '__main__':
#     # Import local versions of submodules
#     print('Notice: Loading local version of eispac.read_template')
#     from read_template import EISFitTemplate, read_template
# else:
#     # Import from installed package
#     from eispac.read_template import EISFitTemplate, read_template

def generate_astropy_model(template):
    """Create an Astropy fitting model using a template intended for MPFIT.

    Parameters
    ----------
    template : EISFitTemplate object or string
        Either a single MPFIT-style template object or a string containing the
        full path to a template file. The template MUST contain both 'template'
        and 'parinfo' attributes.

    Returns
    -------
    combined_model : CompoundModel object from Astropy.modeling
        A fully initialized CompoundModel object with the same parameter values
        and constraints (including tied parameter functions) as given in the
        'parinfo' attribute of the input template.
    """
    # Define list of forbidden strings in "tied" parameter expressions.
    # This should improve security and reduce the scope of possible abuse.
    black_list = ['import', 'def', 'class', '@', 'sys.', 'os.', 'subprocess.']

    # Validate input and, if needed, load the template file
    if isinstance(template, str):
        use_template = read_template(template, quiet=True)
    elif isinstance(template, EISFitTemplate):
        use_template = template
    else:
        raise TypeError('Please input either a valid template filename'
                        ' or an EISFitTemplate object.')

    # Extract basic model information from the template
    parinfo = use_template.parinfo
    num_gauss = use_template.template['n_gauss']
    poly_degree = use_template.template['n_poly'] - 1 # n_poly is just number of terms
    num_params = len(parinfo)
    num_subcomp = num_gauss + 1 if poly_degree > -1 else num_gauss

    # Create each of the component models and add them together
    m_set = [] # list of individual astropy fitting models
    p = 0 # Counter for current parinfo parameter number
    for m in range(num_subcomp):
        if m < num_gauss:
            # Gaussian component with values of (amplitude, mean, stddev)
            m_set.append(models.Gaussian1D(parinfo[p]['value'],
                                           parinfo[p+1]['value'],
                                           parinfo[p+2]['value']))
        elif m == num_gauss and poly_degree > -1:
            # Polynomial background profile
            # TODO: double check the order of the template polynomial coeffs
            poly_coeffs = {'c'+str(poly_degree-i) : parinfo[p+i]['value']
                           for i in range(poly_degree+1)}
            m_set.append(models.Polynomial1D(degree=poly_degree,
                                             domain=None, window=None,
                                             **poly_coeffs))

        # loop over all model parameters and set basic constraints
        for param in m_set[m].param_names:
            if parinfo[p]['fixed'] == 1:
                m_set[m].__dict__[param].fixed = True
            if parinfo[p]['limited'][0] == 1:
                m_set[m].__dict__[param].min = parinfo[p]['limits'][0]
            if parinfo[p]['limited'][1] == 1:
                m_set[m].__dict__[param].max = parinfo[p]['limits'][1]
            p += 1

        # Add the new model to the combined model
        if m == 0:
            combined_model = m_set[m]
        else:
            combined_model = combined_model + m_set[m]

    # Check for tied parameters and convert MPFIT-style strings into functions
    generic_names = ['p['+str(i)+']' for i in range(num_params)]
    param_names = combined_model.param_names
    for p in range(num_params):
        tied_str = str(parinfo[p]['tied'])
        if any(gname in tied_str for gname in generic_names):
            if any(forbidden in tied_str for forbidden in black_list):
                raise ValueError('Unsupported code found in tied parameter'
                                 ' string. \nReminder: only simple mathematical'
                                 ' expressions and generic parameter names are'
                                 ' permitted.')
            tied_str = tied_str.replace('^', '**') # Translate pow notation
            for g in range(num_params):
                tied_str = tied_str.replace(generic_names[g], 'model.'+param_names[g])
            tied_func = eval('lambda model : '+tied_str)
            combined_model.__dict__[param_names[p]].tied = tied_func

    return combined_model

if __name__ == '__main__':

    import pathlib

    filename = './templates/eis_template_dir/fe_12_195_119.2c.template.h5'
    filename = str(pathlib.Path(filename).resolve())
    Fe_XII_195_119 = read_template(filename)

    model = generate_astropy_model(Fe_XII_195_119)
