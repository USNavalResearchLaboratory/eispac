import matplotlib.pyplot as plt
import astropy.units as u
import eispac

if __name__ == '__main__':
    # input data and template files
    data_filepath = '../data/eis_20190404_131513.data.h5'
    template_filepath = '../templates/eis_template_dir/fe_12_195_119.2c.template.h5'

    # read fit template
    tmplt = eispac.EISFitTemplate.read_template(template_filepath)

    # Read spectral window into an EISCube
    data_cube = eispac.read_cube(data_filepath, tmplt.central_wave)

    # Select a cutout of the raster (note the order of array & plotting indices!)
    # Note: we want the full wavelength axis, so we set the range far outside the bounds of EIS
    cutout_extent = [48, 165, 254, 378] # units of [arcsec]
    # w_coords = data_cube.axis_world_coords('em.wl')
    lower_left = (cutout_extent[2]*u.arcsec, cutout_extent[0]*u.arcsec, 0*u.angstrom)
    upper_right = (cutout_extent[3]*u.arcsec, cutout_extent[1]*u.arcsec, 1000*u.angstrom)
    raster_cutout = data_cube.crop_by_coords(lower_left, upper_corner=upper_right)

    # Fit the data, then save it to disk and test loading it back in
    fit_res = eispac.fit_spectra(raster_cutout, tmplt, ncpu='max')
    save_filepaths = eispac.save_fit(fit_res, save_dir='cwd')
    load_fit = eispac.read_fit(save_filepaths[0])

    # Extract array of total data and fit intensites
    sum_data_inten = raster_cutout.sum_spectra().data
    fit_wave_cube, fit_inten_cube = fit_res.get_fit_profile(component=[0,1])
    sum_fit_inten = fit_inten_cube.sum(axis=2)

    # Extract example fit profiles at a higher spectral resolution than the data
    ex_coords = [43, 28] # [Y,X] array coords in units of [pixels]
    fit_x, fit_y = fit_res.get_fit_profile(coords=ex_coords, num_wavelengths=100)
    c0_fit_x, c0_fit_y = fit_res.get_fit_profile(component=0, coords=ex_coords,
                                                 num_wavelengths=100)
    c1_fit_x, c1_fit_y = fit_res.get_fit_profile(component=1, coords=ex_coords,
                                                 num_wavelengths=100)
    c2_fit_x, c2_fit_y = fit_res.get_fit_profile(component=2, coords=ex_coords,
                                                 num_wavelengths=100)
    sub_data = raster_cutout.data[ex_coords[0], ex_coords[1], :]
    sub_wave = raster_cutout.wavelength[ex_coords[0], ex_coords[1], :]
    sub_err = raster_cutout.uncertainty.array[ex_coords[0], ex_coords[1], :]

    # Make a multi-panel figure with the cutout and example
    fig = plt.figure()
    plot_grid = fig.add_gridspec(nrows=2, ncols=2, hspace=0.5, wspace=0.3)

    data_subplt = fig.add_subplot(plot_grid[0,0])
    data_subplt.imshow(sum_data_inten, origin='lower', extent=cutout_extent, cmap='gray')
    data_subplt.set_title('Data Cutout')
    data_subplt.set_xlabel('Solar-X [arcsec]')
    data_subplt.set_ylabel('Solar-Y [arcsec]')

    fit_subplt = fig.add_subplot(plot_grid[0,1])
    fit_subplt.imshow(sum_fit_inten, origin='lower', extent=cutout_extent, cmap='gray')
    fit_subplt.set_title('Total Fit Intensity')
    fit_subplt.set_xlabel('Solar-X [arcsec]')
    fit_subplt.set_ylabel('Solar-Y [arcsec]')

    profile_subplt = fig.add_subplot(plot_grid[1,:])
    profile_subplt.errorbar(sub_wave, sub_data, yerr=sub_err,
                           ls='', marker='o', color='k')
    profile_subplt.plot(fit_x, fit_y, color='b', label='Combined profile')
    profile_subplt.plot(c0_fit_x, c0_fit_y, color='r', label='Gaussian 1')
    profile_subplt.plot(c1_fit_x, c1_fit_y, color='r', ls='--', label='Gaussian 2')
    profile_subplt.plot(c2_fit_x, c2_fit_y, color='g', label='Background')
    profile_subplt.set_title(f'Cutout indices iy = {ex_coords[0]}, ix = {ex_coords[1]}')
    profile_subplt.set_xlabel('Wavelength [$\AA$]')
    profile_subplt.set_ylabel('Intensity ['+raster_cutout.unit.to_string()+']')
    profile_subplt.legend(loc='upper left')
    plt.show()
