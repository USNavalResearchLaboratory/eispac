import pathlib
import numpy as np
import h5py

__all__ = ['EISFitTemplate']


class EISFitTemplate:
    """
    Representation of fitting parameters for a particular line or lines

    Parameters
    ----------
    filename
    template
    parinfo
    """

    def __init__(self, filename, template, parinfo):
        self.filename_temp = str(filename)
        self.template = template
        self.parinfo = parinfo

    def __repr__(self):
        rows = []
        rows.append('--- FIT TEMPLATE PARAMETER CONSTRAINTS ---')
        rows.append(f"{'*':>4} {'Value':>16} {'Fixed':>10} "
                    +f"{'Limited':>18} {'Limits':>22} {'Tied':>18}")
        for i in range(len(self.parinfo)):
            rows.append(
                f"{'p['+str(i)+']':>6} "
                    +f"{self.parinfo[i]['value']:14.4f} "
                    +f"{self.parinfo[i]['fixed']:10d} "
                    +f"{self.parinfo[i]['limited'][0]:10d} "
                    +f"{self.parinfo[i]['limited'][1]:10d} "
                    +f"{self.parinfo[i]['limits'][0]:12.4f} "
                    +f"{self.parinfo[i]['limits'][1]:12.4f} "
                    +f"{self.parinfo[i]['tied']:>18}")
        return '\n'.join(rows)

    @property
    def funcinfo(self):
        return self.get_funcinfo(self.template)

    @property
    def central_wave(self):
        return self.template['wmin'] + (self.template['wmax'] - self.template['wmin']) * 0.5


    @staticmethod
    def get_funcinfo(template):
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
        filename : `str`, `pathlib.Path`
            Path to template file
        """
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
