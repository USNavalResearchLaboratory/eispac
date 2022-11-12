.. _ex-basic:

Basic Example
================

Minimal example program that just loads an observation, fits all the spectra,
and saves the results.

.. code:: python

   import eispac

   if __name__ == '__main__':
       # input data and template files
       data_filepath = './eis_20190404_131513.data.h5'
       template_filepath = './fe_12_195_119.2c.template.h5'

       # read fit template
       tmplt = eispac.read_template(template_filepath)

       # Read spectral window into an EISCube
       data_cube = eispac.read_cube(data_filepath, tmplt.central_wave)

       # Fit the data, then save it to disk and test loading it back in
       fit_res = eispac.fit_spectra(data_cube, tmplt, ncpu='max')
       save_filepaths = eispac.save_fit(fit_res, save_dir='cwd')
       FITS_file = eispac.export_fits(fit_res, save_dir='cwd')
       load_fit = eispac.read_fit(save_filepaths[0])
