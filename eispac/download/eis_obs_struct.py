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

        # Create indexes for tl_id (should speed up mk_list and mk_list_main)
        self.cur.execute('CREATE INDEX IF NOT EXISTS main_tl_idx ON eis_main (tl_id);')
        self.cur.execute('CREATE INDEX IF NOT EXISTS exp_tl_idx ON eis_experiment (tl_id);')

        # Create lists of column names (for refernce and easy query generation)
        self.cur.execute('PRAGMA table_info(eis_main);')
        self.main_cols = [col['name'] for col in self.cur.fetchall()]
        self.cur.execute('PRAGMA table_info(eis_experiment);')
        self.exp_cols = [col['name'] for col in self.cur.fetchall()]

        self.load_ll()
        self.load_raster()

    def load_ll(self):
        """Load useful info from line linelist db."""
        ll_string = "id, acronym, title, n_lines, wavelength"
        self.cur.execute("select " + ll_string + " from eis_linelist_db")
        self.ll = self.cur.fetchall()

    def load_raster(self):
        """Load useful info from raster db."""
        rast_string = """id, acronym, title, rastertype, scan_fm_nsteps,
           scan_fm_stepsize, sns_nexps, sns_duration, ll_id, n_windows,
           data_windows, wind_height, slit_index, nexp, exposures"""
        self.cur.execute("select " + rast_string + " from eis_raster_db")
        self.rast = self.cur.fetchall()

    def mk_list_exp(self):
        """Build up full info on selection, converting some items to
        nicer form. Assume experiment db was searched first.
        """
        exp_rows = self.cur.fetchall()

        unknown_main_row = {'stud_acr':'unknown_trigger', 'study_id':0, 'jop_id':0, 
                            'obstitle':' ', 
                            'obs_dec':'Triggered raster (see previous obs) ', 
                            'sci_obj':'AR? ', 'target':'Flare? '}

        self.eis_str = []
        for e_row in exp_rows:
            if e_row['rast_id'] <= 0:
                continue
            # Get needed row from raster
            rast_row = None
            rast_id = e_row['rast_id']
            for item in self.rast:
                if item['id'] == rast_id:
                    rast_row = item

            if rast_row is None:
                continue
            # Get needed row from linelist
            ll_id = rast_row['ll_id']
            for item in self.ll:
                if item['id'] == ll_id:
                    ll_row = item

            # Get needed row from from main
            if e_row['tl_id'] == 1:
                # Triggered rasters (NOT IN MAIN!)
                m_row = unknown_main_row
                m_row['obstitle'] = f"Triggered {e_row['rast_acr']} (rast_acr) "
            else:
                main_string = """stud_acr, study_id, jop_id, obstitle,
                            obs_dec, sci_obj, target"""
                self.cur.execute("select " + main_string + """ from eis_main
                                where tl_id = ?""", (e_row['tl_id'],))
                m_row, = self.cur.fetchall()
            self.eis_str.append(EIS_Struct(e_row, ll_row, rast_row, m_row))

    def mk_list_main(self):
        """Build up full info on selection, converting some items to
        nicer form. Assume main db was searched first.
        """
        main_rows = self.cur.fetchall()
        self.eis_str = []
        for m_row in main_rows:
            # Get experiment data
            exp_string = """filename, date_obs, date_end, xcen, ycen,
                         fovx, fovy, tl_id, rast_acr, rast_id"""
            self.cur.execute("select " + exp_string + """ from
                             eis_experiment where tl_id==?""",
                             (m_row['tl_id'],))
            exp_rows = self.cur.fetchall()
            for e_row in exp_rows:
                if e_row['rast_id'] <= 0:
                    continue
                # Get needed row from raster
                rast_row = None
                rast_id = e_row['rast_id']
                for item in self.rast:
                    if item['id'] == rast_id:
                        rast_row = item

                if rast_row is None:
                    continue
                # Get needed row from linelist
                ll_id = rast_row['ll_id']
                for item in self.ll:
                    if item['id'] == ll_id:
                        ll_row = item

                self.eis_str.append(EIS_Struct(e_row, ll_row, rast_row, m_row))

    def query_main(self, **kwargs):
        """Retrieve info from eis_main using various search criteria."""

        select_main = """tl_id, stud_acr, study_id, jop_id, obstitle,
                      obs_dec, sci_obj, target"""

        start_text_cols = ['stud_acr', 'target']
        mid_text_cols = ['sci_obj', 'obstitle']
        k_vals = tuple()
        query = 'SELECT '+select_main+' FROM eis_main WHERE'
        next_sep = ' '

        # First, parse (and pop out) keys that require special handling
        if 'date' in kwargs:
            input_date = kwargs.pop('date')
            if isinstance(input_date, (list, tuple)):
                if input_date[0] == '' or input_date[0] is None:
                    # don't bother searching by date
                    pass
                else:
                    t0 = utc2tai(input_date[0])
                    if input_date[1] == '' or input_date[1] is None:
                        t1 = t0 + 86400.0 # just add a day
                    else:
                        t1 = utc2tai(input_date[1])
                    k_vals = (t0, t1,)
                    query = query+next_sep+'date_obs BETWEEN ? and ?'
                    next_sep = ' AND '

        # Loop over all remaining eis_main keys/cols
        for key in kwargs:
            if key in self.main_cols:
                if key in start_text_cols:
                    k_vals = k_vals + (kwargs[key]+'%',)
                    query = query+next_sep+key+' LIKE ?'
                elif key in mid_text_cols:
                    k_vals = k_vals + ('%'+kwargs[key]+'%',)
                    query = query+next_sep+key+' LIKE ?'
                else:
                    k_vals = k_vals + (kwargs[key],)
                    query = query+next_sep+key+'=?'
                next_sep = ' AND '

        # Run the actual search query and then assemble list of rasters
        self.cur.execute(query, k_vals)
        self.mk_list_main()

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
        self.mk_list_exp()

    def get_by_rast_id(self, rast_id, date=None):
        """Retrieve info from eis_experiment for raster id."""
        exp_string = """filename, date_obs, date_end, xcen, ycen,
                     fovx, fovy, tl_id, rast_acr, rast_id"""
        if date == None:
            self.cur.execute("select " + exp_string + """ from
                             eis_experiment where rast_id=?""",
                             (rast_id,))
            self.mk_list_exp()
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
            self.mk_list_exp()

    def get_by_study_id(self, study_id, date=None):
        """Retrieve all the executions of a particular study id."""

        # First we do the main database
        main_string = """tl_id, stud_acr, study_id, jop_id, obstitle,
                      obs_dec, sci_obj, target"""
        if date == None:
            self.cur.execute("select " + main_string + """ from
                             eis_main where study_id=?""",
                             (study_id,))
            self.mk_list_main()
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
            self.mk_list_main()

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
            self.mk_list_main()
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
            self.mk_list_main()

    def get_by_user_sql(self, sql, show_sql=False):
        """Retrieve info using arbitary user sql string."""
        exp_string = """filename, date_obs, date_end, xcen, ycen,
            fovx, fovy, tl_id, rast_acr, rast_id """

        sql_string = "select " + exp_string + sql
        if show_sql:
            print(sql_string)
        self.cur.execute(sql_string)
        self.mk_list_exp()


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
