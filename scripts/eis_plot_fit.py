#!/usr/bin/env python

import matplotlib.pyplot as plt
import sys
import eispac
import numpy as np
from pathlib import Path

speed_of_light_km = 299792.458

vmin = -40. # velocity min scaling in km/s
vmax = 40.  # velocity max scaling in km/s
wmin = 25   # width min scaling in mA
wmax = 45   # width max scaling in mA

class eis_plot_fit:

    def __init__(self, eis_fit_file=None):
        if eis_fit_file is not None:
            self.eis_fit_file = Path(eis_fit_file)
            self.read_fit()
            self.read_template()
            self.read_cube()
            self.display_fit()

    def check_file(self, filename):
        if not filename.is_file():
            print(f' ! {filename} not found, exiting . . . ')
            exit()

    def read_fit(self):
        self.check_file(self.eis_fit_file)
        self.fit = eispac.read_fit(self.eis_fit_file)

    def read_template(self):
        self.eis_template_file = Path(self.fit.meta['filename_template'])
        self.check_file(self.eis_template_file)
        self.template = eispac.read_template(self.eis_template_file)

    def read_cube(self):
        self.eis_data_file = Path(self.fit.meta['filename_data'])
        self.check_file(self.eis_data_file)
        self.cube = eispac.read_cube(self.eis_data_file, self.template.central_wave)

    def scale_intensity(self, intensity):
        range = np.percentile(intensity, (1,99))
        range = [range[1]*1.0E-2, range[1]]
        if range[0] < 10: range[0] = 10.0
        scaled = np.clip(intensity, range[0], range[1])
        scaled = np.log10(scaled)
        imin = np.log10(range[0])
        imax = np.log10(range[1])
        scaled = (scaled-imin)/(imax-imin)
        return scaled

    def scale_velocity(self, centroid, line_id):
        wave0 = float(line_id.split()[2])
        velocity = speed_of_light_km*(centroid-wave0)/wave0
        med = np.median(velocity)
        velocity = velocity - med
       
        scaled = np.clip(velocity, vmin, vmax)
        scaled = (scaled-vmin)/(vmax-vmin)
        return scaled

    def scale_width(self, width):
        width *= 1000.
        scaled = np.clip(width, wmin, wmax)
        scaled = (scaled-wmin)/(wmax-wmin)
        return scaled        

    def display_fit(self):
        component = self.fit.fit['main_component'] # the line of interest, could be multiple
        line_id = self.template.template['line_ids'][component] # the line id for this line
        
        intensity = self.fit.fit['int'][:,:,component] # the integrated line intensity
        intensity_error = self.fit.fit['err_int'][:,:,component]
        
        centroid, error_centroid = self.fit.get_params(param_name='centroid') # fit centroid
        centroid = centroid[:,:,component]
        
        width, error_width = self.fit.get_params(param_name='width') # fit gaussian width
        width = width[:,:,component]

        x_scale = self.cube.meta['pointing']['x_scale'] # arcsec per x steps
        y_scale = self.cube.meta['pointing']['y_scale'] # arcsec per y, always 1!
        aspect = y_scale/x_scale

        fig, (ax1, ax2, ax3) = plt.subplots(1,3, figsize=(10,5))
        
        scaled_int = self.scale_intensity(intensity)
        scaled_vel = self.scale_velocity(centroid, line_id)
        scaled_wdh = self.scale_width(width)

        ax1.imshow(scaled_int, aspect=aspect, cmap='inferno', origin='lower', clim=(0,1))
        ax1.set_title('Intensity')
        ax1.set_xlabel('Solar X (pixels)')
        ax1.set_ylabel('Solar Y (pixels)')
        
        ax2.imshow(scaled_vel, aspect=aspect, cmap='RdBu_r', origin='lower', clim=(0,1))
        ax2.set_title('Velocity')
        ax2.set_xlabel('Solar X (pixels)')
        
        ax3.imshow(scaled_wdh, aspect=aspect, cmap='viridis', origin='lower', clim=(0,1))
        ax3.set_title('Line Width')
        ax3.set_xlabel('Solar X (pixels)')        

        t = self.eis_fit_file.stem + ' | ' + line_id
        fig.text(0.5, 0.96, t, ha='center', va='center')

        opf = self.eis_fit_file.with_suffix('.png')
        plt.savefig(opf, dpi=200)
        print(f' + saved {opf}')
        plt.show()

def main():
    if len(sys.argv) != 2:
        print(' ! input the name of an EIS fit file')
        print(' > eis_plot_fit data_eis/eis_20190404_131513.fe_12_195_119.2c-0.fit.h5')
        exit()

    o = eis_plot_fit(sys.argv[1])
        
if __name__ == '__main__':
    main()
