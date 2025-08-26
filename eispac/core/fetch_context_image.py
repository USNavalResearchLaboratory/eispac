__all__ = ['fetch_context_image']

import os
import sys
import pathlib
import numpy as np
import asyncio
import astropy.units as u
import sunpy.map
from sunpy import __version__ as sunpy_ver
from sunpy.time import parse_time
from sunpy.net import Fido, attrs as a
from astropy.visualization import ImageNormalize, AsinhStretch
from eispac.core.eiscube import EISCube
from eispac.core.eisfitresult import EISFitResult

def fetch_context_image(eis_obs, wavelength=193, units=u.Angstrom, 
                        near='avg', instrument='aia', source='vso', 
                        jsoc_email=None, save_dir=None, overwrite=False, 
                        return_map=True, quiet=False):
    """Download an SDO/AIA or SOHO/EIT context image for a given EIS observation

    Parameters
    ----------
    eis_obs : `EISCube`, `EISMap`, or `EISFitResult` object
        EIS observation you would like to download a context image for. Only the
        time information is used, therefore the exact spectral window of the data 
        is unimportant. Can accept either an EISCube, EISMap (SunPy map subclass), 
        or EISFitResult object.
    wavelength : str, int_float, or Astropy `Quantity`, optional
        Context wavelength channel to fetch. Default is 193 Angstroms for SDO/AIA 
        and 195 Angstroms for SOHO/EIT.
    units : str or Astropy Unit, optional
        Unit for "wavelength". Ignored if "wavelength" is an Astropy Quanitity.
        Default is "Angstrom". 
    near : str, optional
        Nearest timestamp to use for selecting the context image. Choose from
        "start", "average", or "end" ("beginning" and "stop" are also recognized). 
        Default is "average".
    instrument : str, optional
        Context image instrument. Choose from "eit" or "aia". Default is "eit"
        before 2010-05-13 and "aia" after.
    source : str, optional
        Remote data server to search and download from. Choose from "VSO" 
        (virtual solar observatory) or "jsoc" (joint science operations center).
        The VSO has both AIA & EIT data while JSOC only has AIA.
        NB: before using "jsoc", users MUST register an email with the data
        export system, see <http://jsoc.stanford.edu/ajax/register_email.html>
        for more details. Default is "vso".
    jsoc_email : str, optional
        Registered email to notify if source="jsoc". If set to None, will attempt
        to get the address from a user environment variable named "JSOC_EMAIL" 
        (if available). Default is None.
    save_dir : str or `pathlib.Path`, optional
        Directory to save the downloaded map in. Can also set to "cwd" to select
        the current working directory or "none" to select the same dir as the
        input eis_obs (if an EISMap is input, will use "cwd" instead).
        Default is "cwd"
    overwrite : bool, optional
        If set to "True", will overwrite any previously downloaded context image
        Default is "False".
    return_map : bool, optional
        If set to "True", will return the context image as a SunPy Map instead
        of returning the path to the downloaded file. Default is "True"
    quiet : bool, optional
        If set to "True", will not print final summary of download results.
        Download progress from SunPy Fido will also be hidden. Default is "False" 
        
    Returns
    -------
    context : sunpy.map.Map or str
        Sunpy Map containing the context image or, if "return_map=False", the
        path to the downloaded context image
    """

    # Set reference info for each instrument
    eit_channels = np.array([171, 195, 284, 304])
    aia_channels = np.array([94, 131, 171, 193, 211, 304, 335, 
                             1600, 1700, 4500])
    aia_start_date = parse_time('2010-05-13T00:00')
    
    # Validate input eis_obs and load timestamps
    if isinstance(eis_obs, sunpy.map.GenericMap):
        eis_date_beg = eis_obs.date_start
        eis_date_avg = eis_obs.date_average
        eis_date_end = eis_obs.date_end
        eis_dir = None
    elif isinstance(eis_obs, (EISCube, EISFitResult)):
        eis_date_beg = parse_time(eis_obs.meta['mod_index']['date_beg'])
        eis_date_avg = parse_time(eis_obs.meta['mod_index']['date_avg'])
        eis_date_end = parse_time(eis_obs.meta['mod_index']['date_end'])
        eis_dir = pathlib.Path(eis_obs.meta['filename_head']).resolve().parent
    else:
        print("ERROR: Please input a valid EISCube, EISMap,"
             +" or EISFitResult object", file=sys.stderr)
        return None

    # Validate save dir
    if str(save_dir).casefold() == 'cwd':
        output_dir = pathlib.Path().cwd()
    elif str(save_dir).casefold() == 'none':
        if eis_dir is not None:
            output_dir = eis_dir
        else:
            output_dir = pathlib.Path().cwd()
    elif isinstance(save_dir, (str, pathlib.Path)):
        output_dir = pathlib.Path(save_dir).resolve()
    else:
        output_dir = pathlib.Path('**invalid_save_dir_datatype**')
    
    if not output_dir.is_dir():
        print("ERROR: invalid or missing save directory. Please input a string"
             +"or pathlib.Path to where you would like to save the context image" 
             +"(or set save_dir='cwd' to save to the current working dir).", 
             file=sys.stderr)
        return None
    
    # Parse input units
    if isinstance(units, (str, u.Quantity, u.Unit)):
        wave_unit = u.Unit(units) # default
    else:
        wave_unit = u.second # set to a non-length unit, will raise error later
    
    # Parse input aia wavelength
    if isinstance(wavelength, u.Quantity):
        input_wave = wavelength
        wave_unit = wavelength.unit
    elif isinstance(wavelength, (str, int, float)):
        input_wave = wavelength*wave_unit

    # Check unit type and convert to Angstroms
    if wave_unit.physical_type == 'length':
        input_wave = input_wave.to('Angstrom')
    else:
        print("ERROR: invalid wavelength unit. Please input a length unit such"
             +" as 'Angstrom' or 'nm'.", file=sys.stderr)
        return None
    
    # Validate instrument and correct default wavelength values
    if not isinstance(instrument, str) or instrument.lower() not in ['eit', 'aia']:
        print(f"ERROR: {instrument} is not a recognized instrument."
             +f" Please choose from 'eit' or 'aia.", file=sys.stderr)
        return None
    elif instrument.lower() == 'aia' and eis_date_beg < aia_start_date:
        instrument = 'eit'
        print(f"WARNING: AIA data is unavailable before {aia_start_date.iso}"
             +f" Searching for EIT data instead.", file=sys.stderr)

    # Select correct list of wavelength channels and adjust default wavelength
    if instrument.lower() == 'eit':
        valid_channels = eit_channels
        if input_wave.value == 193:
            input_wave = 195*u.Angstrom
    elif instrument.lower() == 'aia':
        valid_channels = aia_channels
        if input_wave.value == 195:
            input_wave = 193*u.Angstrom # VSO accepts 193-195 for AIA

    # Ensure selected wavelength is a valid wavelength channel
    if input_wave.value not in valid_channels:
        print(f"ERROR: {wavelength} is not a recognized {instrument.upper()} channel."
              +f" Please choose from the following list: {list(valid_channels)}" 
              " Angstroms", file=sys.stderr)
        return None

    # Select timestamp for "near" keyword
    if str(near).lower().startswith(('beg', 'start')):
        eis_date_near = eis_date_beg
    elif str(near).lower().startswith('avg'):
        eis_date_near = eis_date_avg
    elif str(near).lower().startswith(('end', 'stop')):
        eis_date_near = eis_date_end

    # Silence ignored exceptions on Windows (needed due to outdated parfive code)
    if os.name == 'nt':
        asyncio.set_event_loop_policy(asyncio.WindowsSelectorEventLoopPolicy())

    show_progress = True
    if quiet:
        show_progress = False

    # Generate and run the search query
    output_filepath = None # Default return value in case of failed downloads
    output_timestamp = None
    if instrument.lower() == 'eit':
        # SOHO/EIT data from VSO
        time_range = a.Time(eis_date_beg - 12*u.hour,
                            end=eis_date_end + 12*u.hour)
                
        results = Fido.search(time_range & a.Wavelength(input_wave) 
                              & a.Instrument.eit)
        
        # Filter list for nearest file of highest level data available
        if len(results) >= 1:
            if sunpy_ver >= '7.0':
                # Level-1 data (calibrated and pre-processed)
                loc_lv = np.where(results[0,:]['Size'] >= 4*u.Mibyte)
                eit_level = 1
                if len(loc_lv[0]) == 0:
                    # No Level 1 found, fallback to level 0.5
                    print(f"WARNING: No calibrated, level 1 EIT data found. " 
                        +f" Returning level 0.5 data instead.", file=sys.stderr)
                    loc_lv = np.where(results[0,:]['Size'] <= 4*u.Mibyte)
                    eit_level = 0.5
            else:
                # Level 0.5 data (old versions of SunPy cannot load level 1)
                print(f"WARNING: Calibrated, level 1 EIT data is available but" 
                     +f" requires SunPy >= 7.0 (please update if you wish to use"
                     +f" it). Returning Lv 0.5 data instead.", file=sys.stderr)
                loc_lv = np.where(results[0,:]['Size'] <= 3*u.Mibyte)
                eit_level = 0.5

            time_diff = (parse_time(results[0,loc_lv]['Start Time']) 
                         - eis_date_near).value
            nearest_file_index = loc_lv[0][np.argmin(np.abs(time_diff))]

            files = Fido.fetch(results[0,nearest_file_index], path=output_dir, 
                               overwrite=overwrite, progress=show_progress)
            if len(files) > 0:
                output_filepath = files[0]
                output_timestamp = parse_time(results[0,nearest_file_index]['Start Time'])

    elif instrument.lower() == 'aia' and str(source).lower() == 'vso':
        # SDO/AIA data from VSO
        time_range = a.Time(eis_date_beg - 30*u.minute,
                            end=eis_date_end + 30*u.minute,
                            near=eis_date_near)

        results = Fido.search(time_range & a.Wavelength(input_wave) 
                              & a.Physobs.intensity & a.Instrument.aia)
        files = Fido.fetch(results, path=output_dir, overwrite=overwrite, 
                           progress=show_progress)
        if len(files) > 0:
            output_filepath = files[0]
            output_timestamp = parse_time(results[0,0]['Start Time'])

    elif instrument.lower() == 'aia' and str(source).lower() == 'jsoc':
        # SDO/AIA data from JSOC (slower, but should still work when VSO is down)
        if str(jsoc_email).lower() == 'none':
            # Note: since linux is case-sensitive, we must nest the getenv here 
            notify_email = os.getenv('JSOC_EMAIL', os.getenv('jsoc_email', 'unknown'))
        elif isinstance(jsoc_email, str) and '@' in jsoc_email:
            notify_email = jsoc_email
        else:
            print(f"ERROR: {jsoc_email} is not a valid email address."
              +f" Please input an address registered for JSOC export (to register,"
              +f"visit <http://jsoc.stanford.edu/ajax/register_email.html>)",
              file=sys.stderr)
            return None
            
        search_time_range = a.Time(eis_date_near - 30*u.minute,
                                   end=eis_date_near + 30*u.minute)
        results = Fido.search(search_time_range & a.Wavelength(input_wave) 
                              & a.jsoc.Series('aia.lev1_euv_12s')
                              & a.jsoc.Notify(notify_email))
        
        # Note: JSOC does not allow downloading only a subset of a search.
        #       Therefore, we have to filter the files ourselves and run a
        #       smaller search containing ONLY the files we want.
        if len(results) > 0 and len(results[0]) > 0:
            time_diff = (parse_time(results[0,:]['T_REC']) - eis_date_near).value
            nearest_file_index = np.argmin(np.abs(time_diff))
            nearest_timestamp = parse_time(results[0,nearest_file_index]['T_REC'])

            time_range = a.Time(nearest_timestamp - 1*u.second,
                                end=nearest_timestamp + 1*u.second)
            results = Fido.search(time_range & a.Wavelength(input_wave) 
                                & a.jsoc.Series('aia.lev1_euv_12s')
                                & a.jsoc.Notify(notify_email))
        
        # NOW, we finally download the files
        files = Fido.fetch(results, path=output_dir, overwrite=overwrite, 
                           progress=show_progress)
        if len(files) > 1:
            # Delete extra "spikes" file (we do not need or want it)
            for FILEPATH in files:
                if 'spikes' in FILEPATH:
                    os.remove(FILEPATH)
                else:
                    output_filepath = FILEPATH
                    output_timestamp = nearest_timestamp
        elif len(files) == 1:
            output_filepath = files[0]
            output_timestamp = nearest_timestamp
        
    # Print details
    if not quiet:
        if output_filepath is None:
            print(f"Unable to download {instrument.upper()} context image!"
                 +f" Please check your internet connection and input parameters.")
        else:
            print(f"Finished downloading {instrument.upper()} conext image:"
                 +f"\n   EIS date range: {eis_date_beg.iso}  to  {eis_date_end.iso}"
                 +f"\n   Context image: {instrument.upper()} {str(input_wave)}, {output_timestamp.iso}"
                 +f"\n   save_dir: {output_dir}"
                 +f"\n   filename: {pathlib.Path(output_filepath).name}")

    # Output the filepath (default) or Map with the context image
    if return_map and not output_filepath is None:
        if instrument.lower() == 'eit' and eit_level < 1:
            # Loading the EIT level 0.5 data requires fixing the color scaling
            # For a quick fix, we clip the top and bottom 0.5%
            output_map = sunpy.map.Map(output_filepath)
            vmin = np.percentile(output_map.data, 0.5)
            vmax = np.percentile(output_map.data, 99.5)
            output_map.plot_settings['norm'] = ImageNormalize(vmin=vmin, vmax=vmax,
                                                        stretch=AsinhStretch())
            print(f"NOTICE: The color scale for the output EIT level 0.5 map has"
                 +f" been adjusted to clip off the top and bottom 0.5% of data")
            return output_map
        else:
            # Everything else, just load the map
            return sunpy.map.Map(output_filepath)
    else:
        # Will also return "None", if the search yields no results
        return output_filepath