#! /usr/bin/env python

"""Routines to access the EIS asrun sqlite catalog (eis_cat.sqlite).

{NOTE: THIS CLASS IS CURRENTLY UNDERGOING REFACTORING}
The EIS_DB class loads the catalog into memory. Searches on the
catalog then result in instances of an EIS_Struct class, usually
an array of results, from which catalog info can be retrieved.
There is a get_by_date method that mimics some of the
functionality of Peter Young's IDL eis_obs_structure.pro routine.
Other methods are designed to provide results to the
eis_catalog.py program, but are also useful on their own.

Note that searches on the eis experiment database are relatively
fast, but if the eis main database is searched first, the time
goes up, since a main database entry can have several experiment
database entries and there end up being more sqlite calls.

The main program at the end is to demonstrate some of the ways
one might want to access the catalog.

"""

__all__ = ['EISAsRun']

import sys
import os
import re
import time
import pathlib
import sqlite3
import numpy as np
from astropy.table import Table
from eispac.util.convert import tai2utc, utc2tai
from .download_db import download_db
from .find_eis_cat import find_eis_cat
from .get_remote_db_modtime import get_remote_db_modtime

def regexp(expr, item):
    """Returns 1 if item matches the regular expression, 0 otherwise."""
    r = re.compile(expr, re.IGNORECASE)
    return r.match(item) is not None

class EISAsRun():
    """Connect with the EIS as-run sqlite catalog."""

    def __init__(self, filepath=None):
        """Establish connection and load some useful info."""
        
        # Validate input filepath
        input_cat_filepath = None
        if filepath is None:
            print(f'WARNING: Missing path to the EIS as-run catalog!',
                  file=sys.stderr)
        elif not isinstance(filepath, (str, pathlib.Path)): 
            print(f'ERROR: Invalid EIS as-run catalog path! Please input a'
                 +f' valid string or pathlib.Path object.',
                  file=sys.stderr)
        else:
            input_cat_filepath = pathlib.Path(filepath).resolve()
            if (not input_cat_filepath.is_file() 
            or not str(filepath).lower().endswith('.sqlite')):
                input_cat_filepath = None
                print(f'ERROR: Input path does not point to a SQLite file!'
                      +f' Please input a valid string or pathlib.Path object.',
                      file=sys.stderr)

        if input_cat_filepath is None:
            print('Attempting to automatically locate the EIS as-run catalog...')
            input_cat_filepath = find_eis_cat()
            if input_cat_filepath is None:
                # IF still cannot find a valid catalog file, just exit
                self.cat_filepath = None
                return None
        
        # Connect to the catalog file
        self.cat_filepath = input_cat_filepath
        self.conn = sqlite3.connect(input_cat_filepath)
        if sys.version_info[0] < 3:
            self.conn.text_factory = str
        else:
            self.conn.text_factory = lambda x: str(x, 'latin1')
        self.conn.row_factory = sqlite3.Row
        # Allow us to use REGEX
        self.conn.create_function("regexp", 2, regexp)
        self.cur = self.conn.cursor()

        # Get reference timestamps
        # Note: the tl_id filter is needed due to an odd engi raster 
        #       with a June 2025 timestamp (as of Feb of 2025)
        #       Likewise, the ll_id filder is needed for first date_obs
        self.cur.execute('select MAX(date_mod), * from eis_experiment where tl_id > 10000')
        self.last_date_mod = tai2utc(self.cur.fetchall()[-1]['date_mod'])
        self.cur.execute('select MIN(date_obs), * from eis_experiment where date_obs > 0 and ll_id > 0')
        self.first_date_obs = tai2utc(self.cur.fetchall()[0]['date_obs'])
        self.cur.execute('select MAX(date_obs), * from eis_experiment where tl_id > 10000')
        self.last_date_obs = tai2utc(self.cur.fetchall()[-1]['date_obs'])
        self.cur.execute('select MAX(date_end), * from eis_experiment where tl_id > 10000')
        self.last_date_end = tai2utc(self.cur.fetchall()[-1]['date_end'])

        # Create indexes for tl_id (should speed up mk_obs_list)
        self.cur.execute('CREATE INDEX IF NOT EXISTS main_tl_idx ON eis_main (tl_id);')
        self.cur.execute('CREATE INDEX IF NOT EXISTS exp_tl_idx ON eis_experiment (tl_id);')

        # Create lists of column names (for reference and easy query generation)
        self.cur.execute('PRAGMA table_info(eis_main);')
        self.main_cols = [col['name'] for col in self.cur.fetchall()]
        self.cur.execute('PRAGMA table_info(eis_experiment);')
        self.exp_cols = [col['name'] for col in self.cur.fetchall()]

        self.unique_main_cols = []
        for col in self.main_cols:
            if col not in self.exp_cols:
                self.unique_main_cols.append(col)

        self.unique_exp_cols = []
        for col in self.exp_cols:
            if col not in self.main_cols:
                self.unique_exp_cols.append(col)

        # Load useful info from line linelist db
        ll_string = "id, acronym, title, n_lines, wavelength"
        self.cur.execute("select " + ll_string + " from eis_linelist_db")
        self.ll = self.cur.fetchall()

        # Check integrity of LL ID numbers
        self.n_ll_ids = len(self.ll)
        self.ll_sorted = True # assume LL list index == ID num
        for ID in range(self.n_ll_ids):
            if ID != self.ll[ID]['id']:
                self.ll_sorted = False # LL list index != ID num
                break

        # Load useful info from raster db
        rast_string = """id, acronym, title, rastertype, scan_fm_nsteps,
           scan_fm_stepsize, sns_nexps, sns_duration, ll_id, n_windows,
           data_windows, wind_height, slit_index, nexp, exposures"""
        self.cur.execute("select " + rast_string + " from eis_raster_db")
        self.rast = self.cur.fetchall()

        # Check integrity of raster ID numbers
        self.n_rast_ids = len(self.rast)
        self.rast_sorted = True # assume rast list index == ID num
        for ID in range(self.n_rast_ids):
            if ID != self.rast[ID]['id']:
                self.rast_sorted = False # rast list index != ID num
                break

        # Create placeholders
        self.results = []
        self.skipped_obs = []
        self.unknown_main_row = {'stud_acr':'unknown', 'study_id':-999, 'jop_id':-999, 
                                 'obstitle':'NO STUDY INFO FOUND! ', 
                                 'obs_dec':'No study info with same tl_id found within +-12 hrs ', 
                                 'sci_obj':'?? ', 'target':'unknown'}
        self.unknown_exp_row = {'filename':' ', 'date_obs':1, 'date_end':2, 
                                'xcen':-90000, 'ycen':-90000, 'fovx':0.0, 'fovy':0.0, 
                                'tl_id':0, 'rast_acr':'unknown', 'rast_id':1}

    @property
    def eis_str(self):
        """Reformat result table to be (mostly) backwards compatible"""
        return self.results.as_array().view(np.recarray)
    
    def check_date(self, test_date, quiet=False):
        """Check if a given date is AFTER the last date_end in the catalog"""

        if self.cat_filepath is None or self.cat_filepath == '':
            print('ERROR: No EIS as-run catalog loaded!', sys.stderr)
            return False
        
        # Compare the timestamps
        test_utc = np.datetime64(test_date)
        first_obs_utc = np.datetime64(self.first_date_obs)
        last_end_utc = np.datetime64(self.last_date_end)
        if test_utc < first_obs_utc:
            if not quiet:
                print(f'NOTICE: Input date is before the first date_obs found in'
                    +f' the catalog ({self.first_date_obs}). Please check your'
                    +f' date range and try again.')
            return None
        elif test_utc > last_end_utc:
            if not quiet:
                print(f'NOTICE: Input date is after the last date_end found in'
                    +f' the catalog ({self.last_date_end}). New files are'
                    +f' typically available 7-14 days after observation.'
                    +f' Use .update() to download the latest catalog.')
            return False
        else:
            return True

    def update(self, force=False):
        """Download an updated version of the eis_cat.sqlite catalog"""
        
        # Close current connections
        if self.cat_filepath is not None:
            self.cur.close()
            self.conn.close()
        else:
            # Why are you updating when you don't even have a file loaded?
            return None
        
        # Check current time and compare to last mod time
        # If time diff < 24 hours, there is probably NOT a new file to download
        now_utc = np.datetime64(time.strftime('%Y-%m-%dT%H:%M:%S',time.gmtime()))
        mod_utc = np.datetime64(self.last_date_mod)
        next_mod_utc = mod_utc + np.timedelta64(24, 'h')
        if now_utc < next_mod_utc and not force == True:
            wait_time = (next_mod_utc - now_utc)/np.timedelta64(1, 'h')
            print(f'NOTICE: Database was last modified < 24 hrs ago and should'
                 +f' be up-to-date. Please try again in {wait_time:.1f} hours'
                 +f' or use the "force=True" option.')
            self.__init__(self.cat_filepath) # reload current file
            return None
            
        # Attempt to download new cat (starting with NRL hosted file first)
        # Note: the NASA file is mirrored from NRL ~1 day later
        for source, time_buffer in zip(['nrl', 'nasa'], [12, 36]):
            try:
                source_utc = np.datetime64(get_remote_db_modtime(source))
                if source_utc - mod_utc < np.timedelta64(time_buffer, 'h'):
                    print(f'EIS as-run catalog is already up-to-date'
                          f' (last date_mod: {self.last_date_mod})')
                    new_cat_filepath = str(self.cat_filepath)
                    if not force == True:
                        # Don't re-download the same file unless forced to.
                        break
                new_cat_filepath = download_db(self.cat_filepath.parent, 
                                               source=source)
                if not new_cat_filepath.lower().endswith('.part'):
                    # Download successful!
                    break
            except:
                new_cat_filepath = f'cannot_download_{source}_eis_cat.part'
        
        if new_cat_filepath.endswith('.part'):
            print('ERROR: Database download failed! Please check your internet'
                  +' connection.', sys.stderr)
        else:
            self.cat_filepath = new_cat_filepath

        # Load the new catalog (or just reload the old, if download failed)
        self.__init__(self.cat_filepath)

    def log_skipped_obs(self, exp_row, rast_row=None, log_note=' '):
        """Append details to list of skipped obs"""
        new_entry = {'filename':exp_row['filename'], 
                     'tl_id':exp_row['tl_id'], 
                     'date_obs':tai2utc(exp_row['date_obs']), 
                     'date_end':tai2utc(exp_row['date_end']), 
                     'rast_id':exp_row['rast_id'], 
                     'rast_acr':exp_row['rast_acr'], 
                     'log_note':log_note}
        if rast_row is not None:
            new_entry['ll_id'] = rast_row['ll_id']
            new_entry['rastertype'] = rast_row['rastertype']
            new_entry['slit_index'] = rast_row['slit_index']
            new_entry['n_windows'] = rast_row['n_windows']
            new_entry['nexp'] = rast_row['nexp']
            new_entry['acronym'] = rast_row['acronym']
            new_entry['title'] = rast_row['title']
        
        self.skipped_obs.append(new_entry)

    def mk_trigger_main_row(self, tl_id, rast_acr):
        main_row = {'stud_acr':'triggered_study', 'study_id':0, 'jop_id':0, 
                    'obstitle':' ', 'obs_dec':' ', 
                    'sci_obj':' ', 'target':'Flare Site'}
        if tl_id == 1:
            main_row['obs_dec'] = 'XRT flare triggered raster '
            main_row['sci_obj'] = 'AR, FLR '
            main_row['target'] = 'Flare Site'
        elif tl_id == 3:
            main_row['obs_dec'] = 'EIS flare triggered raster '
            main_row['sci_obj'] = 'AR, FLR '
            main_row['target'] = 'Active Region'
        elif tl_id == 4:
            main_row['obs_dec'] = 'EIS bright point triggered raster '
            main_row['sci_obj'] = 'QS, BP '
            main_row['target'] = 'Quiet Sun'

        main_row['obstitle'] = f"{rast_acr} (rast_acr) trigger response "
        
        return main_row

    def mk_obs_list(self):
        """Make merged list with info for each selected observation.
        Will work regardless of which DB table was searched first. 
        """
        # First, figure out which DB was searched
        unknown_rows = self.cur.fetchall()
        main_rows = []
        exp_rows = []
        if len(unknown_rows) == 0:
            primary_db = 'N/A'
        elif 'study_id' in unknown_rows[0].keys():
            main_rows = unknown_rows
            exp_rows = [0]
            primary_db = 'main'
        elif 'rast_id' in unknown_rows[0].keys():
            main_rows = [0]
            exp_rows = unknown_rows
            primary_db = 'experiment'

        self.results = []
        self.skipped_obs = []
        for loop_m_row in main_rows:
            if primary_db.lower().startswith('main'):
                # Search the EXPERIMENT DB for info
                t_start = loop_m_row['date_obs'] - 43200 # 12 hr BEFORE
                t_end = loop_m_row['date_end'] + 43200 # 12 hr AFTER
                exp_string = ('filename, date_obs, date_end, xcen, ycen,'
                                +' fovx, fovy, tl_id, rast_acr, rast_id')
                self.cur.execute('SELECT '+exp_string+' FROM'
                                +' eis_experiment WHERE tl_id==?'
                                +' AND date_obs BETWEEN ? and ?',
                                (loop_m_row['tl_id'], t_start, t_end))
                exp_rows = self.cur.fetchall()
                if len(exp_rows) <= 0:
                    empty_e_row = self.unknown_exp_row
                    empty_e_row['tl_id'] = loop_m_row['tl_id']
                    empty_e_row['date_obs'] = loop_m_row['date_obs']
                    empty_e_row['date_end'] = loop_m_row['date_end']
                    exp_rows = [empty_e_row]
            
            for e_row in exp_rows:
                # Validate rast ID
                exp_rast_id = e_row['rast_id']
                if exp_rast_id <= 0:
                    # skip engineering studies?
                    self.log_skipped_obs(e_row, log_note='rast_id <= 0')
                    continue
                elif not self.rast_sorted:
                    pass # if rast DB is not clean, don't check vs max ID
                elif exp_rast_id >= self.n_rast_ids:
                    # skip rows with an invalid rast ID
                    self.log_skipped_obs(e_row, log_note='rast_id > MAX(rast_id)')
                    continue

                # Select info from raster DB
                if self.rast_sorted:
                    rast_row = self.rast[exp_rast_id]
                else:
                    # If the rast DB is not clean, loop and find row
                    rast_row = None
                    for item in self.rast:
                        if item['id'] == exp_rast_id:
                            rast_row = item
                    if rast_row is None:
                        # rast ID not found!
                        self.log_skipped_obs(e_row, log_note='rast_id not found')
                        continue
                
                # Validate lineline ID
                ll_id = rast_row['ll_id']
                if ll_id <= 0:
                    self.log_skipped_obs(e_row, rast_row=rast_row,
                                            log_note='ll_id <= 0')
                    continue
                elif not self.ll_sorted:
                    pass # if ll DB is not clean, don't check vs max ID
                elif ll_id >= self.n_ll_ids:
                    # skip rows with an invalid linelist ID
                    self.log_skipped_obs(e_row, rast_row=rast_row,
                                        log_note='ll_id > MAX(ll_id)')
                    continue

                # Select info from linelist DB
                if self.ll_sorted:
                    ll_row = self.ll[ll_id]
                else:
                    ll_row = None
                    for item in self.ll:
                        if item['id'] == ll_id:
                            ll_row = item
                    if ll_row is None:
                        # linelist ID not found!
                        self.log_skipped_obs(e_row, rast_row=rast_row,
                                                log_note='ll_id not found')
                        continue

                # Select the correct main row for this obs 
                exp_tl_id = e_row['tl_id']
                if exp_tl_id in [1, 3, 4]:
                    # Triggered rasters (NOT IN MAIN DATABASE!)
                    # Will overide the incorrect (or missing) main info!
                    m_row = self.mk_trigger_main_row(exp_tl_id, e_row['rast_acr'])
                elif primary_db.lower().startswith('main'):
                    # Use the current main DB row
                    m_row = loop_m_row
                elif primary_db.lower().startswith('exp'):
                    # Search the MAIN DB for information
                    t_start = e_row['date_obs'] - 43200 # 12 hr BEFORE
                    t_end = e_row['date_end'] + 43200 # 12 hr AFTER
                    main_string = ('stud_acr, study_id, jop_id, obstitle,'
                                    +' obs_dec, sci_obj, target')
                    self.cur.execute('SELECT '+main_string+ ' FROM eis_main'
                                    +' WHERE tl_id = ?'
                                    +' AND date_obs BETWEEN ? and ?', 
                                    (e_row['tl_id'], t_start, t_end))
                    study_rows = self.cur.fetchall()
                    if len(study_rows) <= 0:
                        # No matching row found in eis_main!
                        m_row = self.unknown_main_row
                    else:
                        m_row = study_rows[0]

                # Extract, merge, and append the obs info from all rows
                # self.results.append(EIS_Struct(e_row, ll_row, rast_row, m_row))

                # NEW, create a list of dicts and copy to an astropy Table
                self.results.append(_make_obs_dict(e_row, ll_row, rast_row, m_row))
        
        if len(self.results) > 0:
            self.results = Table(rows=self.results)
            self.results.sort('date_obs')

    def search(self, noreturn=False, quiet=False, print_sql=False, 
               auto_update=False, **kwargs):
        """Retrieve info from any eis_cat.sqlite DB using multiple criteria.
        Will automatically identify which DB to search first.
        """

        # Determine which DB to search first based on the first unique key found
        primary_db = 'eis_experiment' # Default search DB
        for key in kwargs.keys():
            if key in self.unique_exp_cols:
                primary_db = 'eis_experiment'
                break
            if key in self.unique_main_cols:
                primary_db = 'eis_main'
                break
        
        if not quiet:
            print(f'Searching {primary_db} table first. Please wait ...')
                
        if primary_db.lower().startswith('eis_exp'):
            select_cols = ('filename, date_obs, date_end, xcen, ycen,'
                          +' fovx, fovy, tl_id, rast_acr, rast_id')
            text_cols = ['rast_acr', 'll_acr']
            valid_cols = self.exp_cols
        elif primary_db.lower().startswith('eis_main'):
            select_cols = ('tl_id, date_obs, date_end, stud_acr, study_id,'
                          +' jop_id, obstitle, obs_dec, sci_obj, target')
            text_cols = ['stud_acr', 'obstitle', 'obs_dec', 'target', 'sci_obj',
                         'join_sb', 'solb_sci', 'eis_sc', 
                         'st_auth', 'observer', 'planner', 'tohbans']
            valid_cols = self.main_cols
        
        query = 'SELECT '+select_cols+' FROM '+primary_db+' WHERE'
        k_vals = tuple()
        next_sep = ' '

        # Parse (and pop out) keys that require special handling (e.g. date)
        # Note: valid dates NOT in a list will be converted as needed
        t0 = None
        t1 = None
        if 'date' in kwargs:
            input_date = kwargs.pop('date')
            if isinstance(input_date, int):
                input_date = list(input_date) # Note: use with care!
            if isinstance(input_date, str):
                if ',' in input_date:
                    input_date = input_date.split(',')
                else:
                    input_date = [input_date]
            if isinstance(input_date, (list, tuple)):
                if input_date[0] is None or input_date[0] == '':
                    # don't bother searching by date
                    pass
                else:
                    t0 = utc2tai(str(input_date[0]))
                    if len(input_date) == 1:
                        # No end time given
                        t1 = t0 + 86400.0 # just add a day
                    elif (len(str(input_date[1])) < 4 
                          or str(input_date[1]).lower() == 'none'):
                        # Invalid end time dtype or length
                        t1 = t0 + 86400.0 # just add a day
                    else:
                        t1 = utc2tai(str(input_date[1]))
                    k_vals = (t0, t1,)
                    query = query+next_sep+'date_obs BETWEEN ? and ?'
                    next_sep = ' AND '
            else:
                print(f'ERROR: {input_date} is not a valid date string!'
                     +f' Please input a date similar to YYYY-MM-DD HH:MM',
                      file=sys.stderr)
                self.mk_obs_list()
                return None # Don't just search the entire database!

        # Loop over all remaining keys, ignoring any that are not in the DB
        for key in kwargs:
            if key in valid_cols:
                if key in text_cols:
                    # Search for case-insensitive sub-string anywhere in text
                    k_vals = k_vals + ('%'+kwargs[key]+'%',)
                    query = query+next_sep+key+' LIKE ?'
                elif ',' in str(kwargs[key]):
                    # Search for values in a list
                    num_tup = ()
                    clean_kwarg_str = str(kwargs[key])
                    for CHAR in ['[', ']', '(', ')', "'", '"']:
                        clean_kwarg_str = clean_kwarg_str.replace(CHAR,'')
                    for NUM in clean_kwarg_str.split(','):
                        if '-' in NUM:
                            # Expand number range
                            start = int(NUM.split('-')[0])
                            end = int(NUM.split('-')[1]) + 1
                            num_tup = num_tup + tuple(range(start, end))
                        else:
                            num_tup = num_tup + (NUM,)
                    q_list = "("+','.join('?' for i in range(len(num_tup)))+")"
                    k_vals = k_vals + num_tup
                    query = query+next_sep+key+' IN '+q_list
                elif '-' in str(kwargs[key]):
                    # Search for numbers in a range
                    start = str(kwargs[key]).split('-')[0]
                    end = str(kwargs[key]).split('-')[1]
                    k_vals = k_vals + (start, end,)
                    query = query+next_sep+key+' BETWEEN ? and ?'
                else:
                    # Search for a single value
                    k_vals = k_vals + (kwargs[key],)
                    query = query+next_sep+key+'=?'
                next_sep = ' AND '
            elif not quiet:
                print(f'WARNING: {key} is not a column in {primary_db}.'
                     +f' Skipping ...')
                
        if print_sql:
            print('\nSQL query:\n  ', query)
            print('\nParameter list:\n  ', k_vals)

        # Check if date range falls outside current catalog time period
        if t1 is not None:
            is_date_good = self.check_date(tai2utc(t1), quiet=quiet)
            
            if is_date_good == False and auto_update == True:
                print('Attempting to automatically update the as-run catalog...')
                self.update()

        # Run the actual search query and then assemble list of observations
        self.cur.execute(query, k_vals)
        self.mk_obs_list()

        if not quiet:
            print(f'Finished!\n'
                 +f'   {len(self.results)} science files found'
                 +f'   {len(self.skipped_obs)} engineering/other files omitted')
            
        if not noreturn:
            return self.results

    def query_main(self, **kwargs):
        """DEPRICATED: Legacy interface for eis_main. Now an alias for search()
        """

        print(f'DEPRICATION WARNING: .query_main() will be removed soon.'
             +f' Please use .search(**kwargs) instead.', file=sys.stderr)
        self.search(**kwargs)

    def get_by_date(self, t0, t1):
        """DEPRICATED: Retrieve info from eis_experiment for date range. 
        Time is expected to be in ISO standard form.
        """

        print(f'DEPRICATION WARNING: .get_by_date() will be removed soon.'
             +f' Please use .search(date=[t0,t1]) instead.', file=sys.stderr)
        self.search(date=[t0,t1])

    def get_by_rast_id(self, rast_id, date=None):
        """DEPRICATED: Retrieve info from eis_experiment for raster id."""

        print(f'DEPRICATION WARNING: .get_by_rast_id() will be removed soon.'
             +f' Please use .search(rast_id, [date]) instead.', file=sys.stderr)
        self.search(rast_id=rast_id, date=date)

    def get_by_study_id(self, study_id, date=None):
        """DEPRICATED: Retrieve all the executions of a particular study id."""

        print(f'DEPRICATION WARNING: .get_by_study_id() will be removed soon.'
             +f' Please use .search(study_id, [date]) instead.', file=sys.stderr)
        self.search(study_id=study_id, date=date)


    def get_by_acronym(self, acronym, date=None):
        """DEPRICATED: Retrieve info using the study acronym. Regular
           expressions are NOT currently supported."""

        print(f'DEPRICATION WARNING: .get_by_acronym() will be removed soon.'
             +f' Please use .search(stud_acr, [date]) instead.', file=sys.stderr)
        self.search(stud_acr=acronym, date=date)

    def get_by_user_sql(self, sql, show_sql=False):
        """Search using arbitary user sql string."""
        exp_string = """filename, date_obs, date_end, xcen, ycen,
            fovx, fovy, tl_id, rast_acr, rast_id """

        sql_string = "select " + exp_string + sql
        if show_sql:
            print(sql_string)
        self.cur.execute(sql_string)
        self.mk_obs_list()

def _make_obs_dict(exp_row, ll_row, rast_row, main_row):
    """Helper function for packing info about an EIS obs into a single dict
    """
    row = {}
    # Info from the experiment database
    row['date_obs'] = tai2utc(exp_row['date_obs'])[0:19]
    row['date_end'] = tai2utc(exp_row['date_end'])[0:19]
    row['filename'] = exp_row['filename']
    row['xcen'] = -9999 if exp_row['xcen'] is None else exp_row['xcen']
    row['ycen'] = -9999 if exp_row['ycen'] is None else exp_row['ycen']
    row['fovx'] = exp_row['fovx']
    row['fovy'] = exp_row['fovy']
    row['tl_id'] = exp_row['tl_id']
    row['rast_acr'] = exp_row['rast_acr']
    row['rast_id'] = exp_row['rast_id']

    # Info from the main database
    row['stud_acr'] = main_row['stud_acr']
    row['study_id'] = main_row['study_id']
    row['jop_id'] = main_row['jop_id']
    row['obstitle'] = main_row['obstitle']
    row['obstitle'] = row['obstitle'][0:-1].rstrip()
    row['obs_dec'] = main_row['obs_dec']
    row['obs_dec'] = row['obs_dec'][0:-1].rstrip()
    row['sci_obj'] = main_row['sci_obj']
    row['sci_obj'] = row['sci_obj'][0:-1].rstrip()
    row['target'] = main_row['target']

    # Info from the raster database
    # row['acronym'] = rast_row['acronym'] # This is a duplicate of rast_acr!!!
    row['title'] = rast_row['title']
    row['rastertype'] = rast_row['rastertype']
    row['scan_fm_nsteps'] = rast_row['scan_fm_nsteps']
    row['scan_fm_stepsize'] = rast_row['scan_fm_stepsize']
    row['sns_nexps'] = rast_row['sns_nexps']
    row['sns_duration'] = rast_row['sns_duration']
    row['ll_id'] = rast_row['ll_id']
    row['n_windows'] = rast_row['n_windows']
    win_pixel_wid = rast_row['data_windows'] # temp variable
    row['wind_height'] = rast_row['wind_height']
    row['slit_index'] = rast_row['slit_index']
    row['nexp'] = rast_row['nexp']
    row['exptime'] = rast_row['exposures']

    # Info from the linelist database
    row['ll_acr'] = ll_row['acronym']
    row['n_lines'] = ll_row['n_lines']
    row['ll_title'] = ll_row['title']
    row['wave'] = ll_row['wavelength']

    # Expand wavelength info
    # Note: One obs can have up to 25 data windows!
    row['wave'] = np.array(row['wave'].split(','), dtype=np.float64) # len = 25
    row['wave'] = row['wave']/100.0
    row['width'] = np.array(win_pixel_wid.split(','), dtype=np.int32) # len = 25
    row['wavemin'] = row['wave'] - row['width']/2.0*0.0223
    row['wavemax'] = row['wave'] + row['width']/2.0*0.0223
    line_list = row['ll_title'].replace(',', ' ')
    line_list = np.array(line_list.split(';'), dtype='U24') # len = n_lines
    # Note: the list of window titles only has N_LINES total elements
    #       Therefore we need to expand it to 25 elements, to be consistent
    row['ll_title'] = np.zeros(25, dtype='U24')
    row['ll_title'][0:row['n_lines']] = line_list # Copy titles back over

    # Clean up exposure time info
    # Note: one obs can have up to 8 different exposures!
    row['exptime'] = np.array(row['exptime'].split(','), dtype=np.float64)
    row['exptime'] = row['exptime']/1000.0 # convert to [s]

    return row

def members(obj):
    """List all the member variables in an object."""
    result = [attr for attr in dir(obj) if not callable(attr)
              and not attr.startswith("__")]
    # print result
    return(result)

def disp_member_info(obj):
    """Print member values in an object"""
    info = vars(obj)
    for item in info.items():
        print(item[0], '=', info[item[0]])


def main():

    t0 = '2015-02-10 15:00'
    t1 = '2015-02-10 20:00'
    dbfile = '/Users/mariska/soft/eis/database/catalog/eis_cat.sqlite'

    d = EISAsRun(dbfile)
    #d.get_by_date(t0, t1)
    #d.get_by_rast_id('51')
    d.get_by_user_sql("""from eis_experiment where
                         rast_id=51 and
                         xcen between -500 and 500 and
                         ycen between -500 and 500 and
                         fovx >= 60""")

    if len(d.eis_str) > 0:
        foo = members(d.eis_str[0])
        print(foo)
        # print vars(d.eis_str[0])
        # disp_member_info(d.eis_str[0])
        print('nrows: ', len(d.eis_str))
        for row in d.eis_str:
            print(row.filename)
            print(row.obstitle)
            print(row.n_lines)
            print(row.wave)
            print(row.width)
            print(row.exptime)
            print('({0:.1f}, {1:.1f})'.format(row.xcen, row.ycen))
    else:
        print('No entries found')

if __name__ == "__main__":
    main()
