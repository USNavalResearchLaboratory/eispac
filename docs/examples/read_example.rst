.. _ex-read:

Reading and Exploring EIS Data
==============================

Reading in a level-1 HDF5 data file and printing some metadata

.. code:: python

   >>> import eispac
   >>> data_filename = 'eis_20190404_131513.data.h5'
   >>> data_cube = eispac.read_cube(data_filename, 195.119)
   >>> data_cube.meta.keys()
   dict_keys(['filename_data', 'filename_head', 'wininfo', 'iwin', 'iwin_str',
              'index', 'pointing', 'wave', 'radcal', 'slit_width',
              'slit_width_units', 'ccd_offset', 'wave_corr', 'wave_corr_t',
              'wave_corr_tilt', 'date_obs', 'date_obs_format', 'duration',
              'duration_units', 'mod_index', 'aspect', 'aspect_ratio', 'notes'])

   >>> data_cube.meta['pointing']['x_scale']
   2.9952

   >>> data_cube.meta['radcal']
   array([8.06751  , 8.060929 , 8.054517 , 8.048271 , 8.042198 , 8.036295 ,
          8.030562 , 8.024157 , 8.017491 , 8.010971 , 8.0046015, 7.998385 ,
          7.9923196, 7.9864078, 7.980654 , 7.975055 , 7.969617 , 7.9643393,
          7.959224 , 7.9542727, 7.949487 , 7.9448686, 7.9404206, 7.9361415],
         dtype=float32)

Viewing information about the spectral windows available in a file.

.. code:: python

   >>> header_filename = 'eis_20190404_131513.head.h5'
   >>> wininfo = eispac.read_wininfo(header_filename)
   >>> wininfo.dtype.names
   ('iwin', 'line_id', 'wvl_min', 'wvl_max', 'nl', 'xs')

   >>> print(wininfo)
   [(0, 'Fe XI 180.400', 180.03426, 180.72559, 32, 661),
    (1, 'Ca XV 182.100', 181.75139, 182.44266, 32, 738),
    (2, 'Fe X 184.720', 183.82512, 185.5865 , 80, 831),
    (3, 'Fe XII 186.750', 186.3891 , 187.0802 , 32, 946)],
    ... ... ...
    (23, 'Mg VII 280.390', 279.7766 , 280.9996 , 56, 3720),
    (24, 'Fe XV 284.160', 283.89   , 284.40134, 24, 3905)]
