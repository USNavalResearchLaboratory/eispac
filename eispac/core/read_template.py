import pathlib
import warnings

import numpy as np
import h5py

__all__ = ['EISFitTemplate', 'read_template']


class EISFitTemplate:
    """
    Representation of fitting parameters for a particular line or lines

    Parameters
    ----------
    filename : str or `pathlib.Path`
        Path to fitting template file
    template : dict
        Dictionary of template parameters
    parinfo : list
        List of fitting parameters where each
        entry is a `dict`
    """

    def __init__(self, filename, template, parinfo):
        self.filename_temp = str(filename)
        self.template = template
        self.parinfo = parinfo

    def __repr__(self):
        rows = []
        rows.append('--- FIT TEMPLATE PARAMETER CONSTRAINTS ---')
        rows.append(
            f"{'*':>4} {'Value':>16} {'Fixed':>10} {'Limited':>18} {'Limits':>22} {'Tied':>18}")
        for i, p in enumerate(self.parinfo):
            rows.append(
                f"{f'p[{i}]':>6} {p['value']:14.4f} {p['fixed']:10d} {p['limited'][0]:10d} "
                + f"{p['limited'][1]:10d} {p['limits'][0]:12.4f} {p['limits'][1]:12.4f} {p['tied']:>18}")
        return '\n'.join(rows)

    @property
    def funcinfo(self):
        """List of dicts specifying each subcomponent function used in the template"""
        return self.get_funcinfo(self.template)

    @property
    def central_wave(self):
        """Wavelength value in the center of the template wavelength range"""
        return self.template['wmin'] + (self.template['wmax'] - self.template['wmin']) * 0.5

    @staticmethod
    def get_funcinfo(template):
        """
        Return a list of dictionaries where each entry describes the
        parameters for one of the fitting basis functions

        Parameters
        ----------
        template : list

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
        """
        Create `EISFitTemplate` from template file

        Parameters
        ----------
        filename : str or `pathlib.Path`
            Path to template file

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

        return cls(filename, template, parinfo)


# This is an alias to make the interface to reading template files a bit more
# user-friendly.
read_template = EISFitTemplate.read_template


if __name__ == '__main__':

    import pathlib

    filename = './templates/eis_template_dir/fe_12_195_119.2c.template.h5'
    filename = pathlib.Path(filename).resolve()
    Fe_XII_195_119 = EISFitTemplate.read_template(filename)

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
