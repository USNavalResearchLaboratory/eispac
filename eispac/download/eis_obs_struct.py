#! /usr/bin/env python

"""Routines to access the EIS asrun sqlite catalog (eis_cat.sqlite).

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

import sys
import os
import re
import sqlite3
import numpy as np
from .convert import tai2utc, utc2tai

def regexp(expr, item):
    """Returns 1 if item matches the regular expression, 0 otherwise."""
    r = re.compile(expr, re.IGNORECASE)
    return r.match(item) is not None

class EIS_DB():
    """Connect with the EIS sqlite catalog."""

    def __init__(self, dbfile):
        """Establish connection and load some useful info."""
        self.conn = sqlite3.connect(dbfile)
        if sys.version_info[0] < 3:
            self.conn.text_factory = str
        else:
            self.conn.text_factory = lambda x: str(x, 'latin1')
        self.conn.row_factory = sqlite3.Row
        # Allow us to use REGEX
        self.conn.create_function("regexp", 2, regexp)
        self.cur = self.conn.cursor()

        # Create indexes for tl_id (should speed up mk_obs_list)
        self.cur.execute('CREATE INDEX IF NOT EXISTS main_tl_idx ON eis_main (tl_id);')
        self.cur.execute('CREATE INDEX IF NOT EXISTS exp_tl_idx ON eis_experiment (tl_id);')

        # Create lists of column names (for refernce and easy query generation)
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

        self.load_ll()
        self.load_raster()

        # Create placeholders
        self.eis_str = []
        self.skipped_obs = []

    def load_ll(self):
        """Load useful info from line linelist db."""
        ll_string = "id, acronym, title, n_lines, wavelength"
        self.cur.execute("select " + ll_string + " from eis_linelist_db")
        self.ll = self.cur.fetchall()

        # Check integrity of ID numbers
        self.n_ll_ids = len(self.ll)
        self.ll_sorted = True # assume LL list index == ID num
        for ID in range(self.n_ll_ids):
            if ID != self.ll[ID]['id']:
                self.ll_sorted = False # LL list index != ID num
                break

    def load_raster(self):
        """Load useful info from raster db."""
        rast_string = """id, acronym, title, rastertype, scan_fm_nsteps,
           scan_fm_stepsize, sns_nexps, sns_duration, ll_id, n_windows,
           data_windows, wind_height, slit_index, nexp, exposures"""
        self.cur.execute("select " + rast_string + " from eis_raster_db")
        self.rast = self.cur.fetchall()

        # Check integrity of ID numbers
        self.n_rast_ids = len(self.rast)
        self.rast_sorted = True # assume rast list index == ID num
        for ID in range(self.n_rast_ids):
            if ID != self.rast[ID]['id']:
                self.rast_sorted = False # rast list index != ID num
                break

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

            self.eis_str = []
            self.skipped_obs = []
            for loop_m_row in main_rows:
                if primary_db.lower().startswith('main'):
                    # Search the EXPERIMENT DB for info
                    exp_string = ('filename, date_obs, date_end, xcen, ycen,'
                                 +' fovx, fovy, tl_id, rast_acr, rast_id')
                    self.cur.execute('SELECT '+exp_string+' FROM'
                                    +' eis_experiment WHERE tl_id==?',
                                    (loop_m_row['tl_id'],))
                    exp_rows = self.cur.fetchall()
                
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
                        main_string = ('stud_acr, study_id, jop_id, obstitle,'
                                      +' obs_dec, sci_obj, target')
                        self.cur.execute('SELECT '+main_string+ ' FROM eis_main'
                                        +' WHERE tl_id = ?', (e_row['tl_id'],))
                        m_row, = self.cur.fetchall()

                    # Extract, merge, and append the obs info from all rows
                    self.eis_str.append(EIS_Struct(e_row, ll_row, rast_row, m_row))

    def search(self, quiet=False, print_sqlite=False, **kwargs):
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
            select_cols = ('tl_id, stud_acr, study_id, jop_id, obstitle,'
                          +' obs_dec, sci_obj, target')
            text_cols = ['stud_acr', 'obstitle', 'obs_dec', 'target', 'sci_obj',
                         'join_sb', 'solb_sci', 'eis_sc', 
                         'st_auth', 'observer', 'planner', 'tohbans']
            valid_cols = self.main_cols
        
        query = 'SELECT '+select_cols+' FROM '+primary_db+' WHERE'
        k_vals = tuple()
        next_sep = ' '

        # Parse (and pop out) keys that require special handling (e.g. date)
        # Note: valid dates NOT in a list will be converted as needed
        if 'date' in kwargs:
            input_date = kwargs.pop('date')
            if isinstance(input_date, int):
                input_date = list(input_date) # Note: use with care!
            if isinstance(input_date, str):
                if ',' in input_date:
                    input_date = input_date.split(',')
                else:
                    input_date = list(input_date)
            if isinstance(input_date, (list, tuple)):
                if input_date[0] == '' or input_date[0] is None:
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
            else:
                print(f'WARNING: {key} is not a column in {primary_db}.'
                     +f' Skipping ...')
                
        if print_sqlite:
            print('\nSQLite query:\n  ', query)
            print('\nSearch parameter list:\n  ', k_vals)

        # Run the actual search query and then assemble list of observations
        self.cur.execute(query, k_vals)
        self.mk_obs_list()

        if not quiet:
            print(f'Finished!\n'
                 +f'   {len(self.eis_str)} science files found'
                 +f'   {len(self.skipped_obs)} engineering/other files omitted')

    def query_main(self, **kwargs):
        """DEPRICATED Legacy interface for eis_main. Now an alias for search()
        """
        self.search(**kwargs)

    def get_by_date(self, t0, t1):
        """Retrieve info from eis_experiment for date range. Time is
        expected to be in ISO standard form.
        """
        self.t0 = utc2tai(t0)
        if t1 == '':
            self.t1 = self.t0 + 86400.0 # just add a day
        else:
            self.t1 = utc2tai(t1)

        exp_string = """filename, date_obs, date_end, xcen, ycen,
                     fovx, fovy, tl_id, rast_acr, rast_id"""

        self.cur.execute("select " + exp_string + """ from
                         eis_experiment where date_obs between
                         ? and ?""", (self.t0, self.t1))
        self.mk_obs_list()

    def get_by_rast_id(self, rast_id, date=None):
        """Retrieve info from eis_experiment for raster id."""
        exp_string = """filename, date_obs, date_end, xcen, ycen,
                     fovx, fovy, tl_id, rast_acr, rast_id"""
        if date == None:
            self.cur.execute("select " + exp_string + """ from
                             eis_experiment where rast_id=?""",
                             (rast_id,))
            self.mk_obs_list()
        else:
            t0 = utc2tai(date[0])
            if date[1] == '':
                t1 = t0 + 86400.0 # just add a day
            else:
                t1 = utc2tai(date[1])
            self.cur.execute("select " + exp_string + """ from
                             eis_experiment where rast_id=? and
                             date_obs between ? and ?""",
                             (rast_id, t0, t1))
            self.mk_obs_list()

    def get_by_study_id(self, study_id, date=None):
        """Retrieve all the executions of a particular study id."""

        # First we do the main database
        main_string = """tl_id, stud_acr, study_id, jop_id, obstitle,
                      obs_dec, sci_obj, target"""
        if date == None:
            self.cur.execute("select " + main_string + """ from
                             eis_main where study_id=?""",
                             (study_id,))
            self.mk_obs_list()
        else:
            t0 = utc2tai(date[0])
            if date[1] == '':
                t1 = t0 + 86400.0 # just add a day
            else:
                t1 = utc2tai(date[1])
            self.cur.execute("select " + main_string + """ from
                             eis_main where study_id=? and
                             date_obs between ? and ?""",
                             (study_id, t0, t1))
            self.mk_obs_list()

    def test_get(self):
        """Play method for testing various ideas."""
        self.cur.execute("""select tl_id, stud_acr, study_id from
                         eis_main where tl_id=58439""")
        r = self.cur.fetchall()
        print(r)
        self.cur.execute("""select tl_id, filename from eis_experiment
                         where tl_id=58439""")
        rr = self.cur.fetchall()
        for row in rr:
            print(row)

    def get_by_acronym(self, acronym, date=None):
        """Retrieve info using the study acronym. Regular
           expressions are fine."""
        main_string = """tl_id, stud_acr, study_id, jop_id, obstitle,
                      obs_dec, sci_obj, target"""
        if date == None:
            self.cur.execute("select " + main_string + """ from
                             eis_main where stud_acr regexp ?""",
                             (acronym,))
            self.mk_obs_list()
        else:
            t0 = utc2tai(date[0])
            if date[1] == '':
                t1 = t0 + 86400.0 # just add a day
            else:
                t1 = utc2tai(date[1])
            self.cur.execute("select " + main_string + """ from
                             eis_main where stud_acr regexp ? and
                             date_obs between ? and ?""",
                             (acronym, t0, t1))
            self.mk_obs_list()

    def get_by_user_sql(self, sql, show_sql=False):
        """Retrieve info using arbitary user sql string."""
        exp_string = """filename, date_obs, date_end, xcen, ycen,
            fovx, fovy, tl_id, rast_acr, rast_id """

        sql_string = "select " + exp_string + sql
        if show_sql:
            print(sql_string)
        self.cur.execute(sql_string)
        self.mk_obs_list()


class EIS_Struct(object):
    """A class to collect data from the databases and put into a
     nice form. May need to use __slots__ at some point, but for
     now, keep simple. Note could also use collections.namedtuple
     for this.
     """

    def __init__(self, exp_row, ll_row, rast_row, main_row):
        """For now passing in result of the timerange query. Should
        be provided with the result of any query and then will fill
        all the values beyond those in the experiment database on its
        own.
        """
        # Query results from experiment database
        self.filename = exp_row['filename']
        self.date_obs = tai2utc(exp_row['date_obs'])[0:19]
        self.date_end = tai2utc(exp_row['date_end'])[0:19]
        self.xcen = -9999 if exp_row['xcen'] is None else exp_row['xcen']
        self.ycen = -9999 if exp_row['ycen'] is None else exp_row['ycen']
        self.fovx = exp_row['fovx']
        self.fovy = exp_row['fovy']
        self.tl_id = exp_row['tl_id']
        self.rast_acr = exp_row['rast_acr']
        self.rast_id = exp_row['rast_id']

        # Things from the raster database
        self.acronym = rast_row['acronym']
        self.title = rast_row['title']
        self.rastertype = rast_row['rastertype']
        self.scan_fm_nsteps = rast_row['scan_fm_nsteps']
        self.scan_fm_stepsize = rast_row['scan_fm_stepsize']
        self.sns_nexps = rast_row['sns_nexps']
        self.sns_duration = rast_row['sns_duration']
        self.ll_id = rast_row['ll_id']
        self.n_windows = rast_row['n_windows']
        self.data_windows = rast_row['data_windows']
        self.wind_height = rast_row['wind_height']
        self.slit_index = rast_row['slit_index']
        self.nexp = rast_row['nexp']
        self.exposures = rast_row['exposures']

        # Things from the linelist database
        self.ll_acronym = ll_row['acronym']
        self.ll_title = ll_row['title']
        self.n_lines = ll_row['n_lines']
        self.wavelength = ll_row['wavelength']

        # Things from the main database
        self.stud_acr = main_row['stud_acr']
        self.study_id = main_row['study_id']
        self.jop_id = main_row['jop_id']
        self.obstitle = main_row['obstitle']
        self.obstitle = self.obstitle[0:-1].rstrip()
        self.obs_dec = main_row['obs_dec']
        self.obs_dec = self.obs_dec[0:-1].rstrip()
        self.sci_obj = main_row['sci_obj']
        self.sci_obj = self.sci_obj[0:-1].rstrip()
        self.target = main_row['target']

        self.add_waves()
        self.add_exptime()

    def add_waves(self):
        """Expand wavelength info."""
        self.wave = np.array(self.wavelength.split(','), dtype=np.float64)
        self.wave = self.wave[0:self.n_lines]
        self.width = np.array(self.data_windows.split(','), dtype=np.int32)
        self.width = self.width[0:self.n_lines]
        self.wavemin = self.wave/100.0 - self.width/2.0*0.0223
        self.wavemax = self.wave/100.0 + self.width/2.0*0.0223
        self.wave = self.wave/100.0
        line_list = self.ll_title.replace(',', ' ')
        self.ll_title = np.array(line_list.split(';'), dtype='U32')

    def add_exptime(self):
        """Clean up exposure time info."""
        self.exptime = np.array(self.exposures.split(','), dtype=np.float64)
        self.exptime = self.exptime[0:self.nexp]/1000.0 # convert to seconds

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

    d = EIS_DB(dbfile)
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
