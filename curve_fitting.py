"""
module for curve fitting assays from a pubchem

this stack overflow has a demo on how to curve fit a function: https://stackoverflow.com/questions/55078451/how-to-use-curvefit-in-python
this guy also: https://gist.github.com/yannabraham/5f210fed773785d8b638
 """
import os
import config
import glob
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import seaborn as sns
import matplotlib.pyplot as plt
import extractor as ext
import ntpath
import sqlite3 as sql
import shutil
import gzip
from curvep import curveP


# CONCISE_DATA_DIR = os.path.join(config.Config.BOX_PATH, 'DATA', 'pubchem', 'bioassay', 'concise_csv', 'all')
# FAILED_CSV_DIR = os.path.join(config.Config.BOX_PATH, 'DATA', 'Nada', 'failed_aids')
# ANOTHER_FAILED_CSV_DIR = os.path.join(config.Config.BOX_PATH, 'DATA', 'Nada', 'another_failed_aids')
# PASSED_CSV_DIR = os.path.join(config.Config.BOX_PATH, 'DATA', 'Nada', 'AID')
# ASSAY_FOLDERS = glob.glob(os.path.join(CONCISE_DATA_DIR, '*'))
# TARGET_ASSAYS = sorted(pd.read_table(os.path.join('data', 'dr_aids.txt'), header=None)[0].values.tolist())
# FAILED_CSV = glob.glob(os.path.join(FAILED_CSV_DIR, '*.csv'))
# ANOTHER_FAILED_CSV = glob.glob(os.path.join(ANOTHER_FAILED_CSV_DIR, '*.csv'))
# PASSED_CSV = glob.glob(os.path.join(PASSED_CSV_DIR, '*.csv'))
# ALL_ASSAY_CSV = ANOTHER_FAILED_CSV + PASSED_CSV
# CONVERTED_CSV = glob.glob(os.path.join(config.Config.BOX_PATH, 'DATA', 'Nada', 'csv_from_json', '*.csv'))
# CONCISE_CSVS = glob.glob(os.path.join("/Volumes", "TOSHIBA EXT", "data", "nada", "concise_dr", "*.csv"))

#DR_AIDS = pd.read_table('data/dr_aids.txt', header=None, names=['AIDS'])['AIDS'].values.tolist()

def hill_curve(conc: float, ac50: float, top: float, slope: float) -> float:
    """ returns """
    denom = np.power(conc, slope) + np.power(ac50, slope)
    return top * (np.power(conc, slope) / denom)


def find_folder(pubchem_aid: int) -> str:
    """ PubChem assays are stored in folders containing 1000 assays.  E.g.,
     0000001_0001000, 0001001_0002000, etc.  This function takes a PubChem AID as input
     and finds the folder that aid would fall in.  """

    # first handle the unique case where
    # the aid falls on the first 1000,
    # e.g., 687000 should be in folder
    # 0686001_0687000
    if (pubchem_aid % 1000) == 0:
        upper_bound = pubchem_aid
        lower_bound = upper_bound - 999
    else:
        # because the folders are stored in 1000 increments,
        # taking the modulus of 1000 gives the remainder needed
        # to round down.  Then, add one.
        lower_bound = pubchem_aid - (pubchem_aid % 1000) + 1
        upper_bound = lower_bound + 1000 - 1

        # throw a test to check the aid is in the
        # folder range
        assert (lower_bound <= pubchem_aid) and (pubchem_aid <= upper_bound)

    # both strings need to be 7 digits long
    string = str(lower_bound).zfill(7) + '_' + str(upper_bound).zfill(7)
    return string


def read_aid(pubchem_aid: int, just_columns=False) -> pd.DataFrame:
    """ given an aid, will find the appropriate folder and read the csv """

    aid_folder = find_folder(pubchem_aid)
    full_path = os.path.join(CONCISE_DATA_DIR, aid_folder, f'{pubchem_aid}.concise.csv')


    # DELETE ME ONLY USED FOR TESTING
    if not os.path.exists(full_path):
        return pd.DataFrame()

    if just_columns:
        return pd.read_csv(full_path, error_bad_lines=False, nrows=0).columns


    # DELETE ME ONLY USED FOR TESTING
    if os.path.exists(full_path) and pd.read_csv(full_path, error_bad_lines=False).empty:
        return pd.DataFrame()

    # the first n rows are a header
    # that describes the data
    # dose response example is 434931
    assay_results = pd.read_csv(full_path, error_bad_lines=False)

    # find index where header stops
    # or find where column is an integer
    #idx = assay_results[assay_results.PUBCHEM_RESULT_TAG == '1'].index[0]
    idx = assay_results[assay_results.PUBCHEM_RESULT_TAG.astype(str).str.isnumeric()].index[0]
    assay_results = assay_results.loc[idx:]

    return assay_results


def read_concise_aid(pubchem_aid_file: str, just_columns=False) -> pd.DataFrame:
    """ given an aid, will find the appropriate folder and read the csv """


    if just_columns:
        return pd.read_csv(pubchem_aid_file, error_bad_lines=False, nrows=0).columns

    # the first n rows are a header
    # that describes the data
    # dose response example is 434931
    assay_results = pd.read_csv(pubchem_aid_file, error_bad_lines=False)

    # find index where header stops
    # or find where column is an integer
    #idx = assay_results[assay_results.PUBCHEM_RESULT_TAG == '1'].index[0]
    idx = assay_results[assay_results.PUBCHEM_RESULT_TAG.astype(str).str.isnumeric()].index[0]
    assay_results = assay_results.loc[idx:]

    return assay_results


def find_unique_columns() -> list:
    """ finds and prints the unique columns across the set of concise csv """
    all_columns = []

    for aid in TARGET_ASSAYS:
        assay_data = read_aid(aid)
        if not assay_data.empty:
            all_columns.extend(list(assay_data.columns))
            print(aid, list(assay_data.columns))

    return list(set(all_columns))

def find_empty_csv() -> int:
    """ in the csv_from_json folder there is a number of csv files that are empty
     this function just counts and returns that number """

    empty_counter = 0

    for assay in CONVERTED_CSV:

        aid = os.path.basename(assay).split('.')[0]
        try:
            df = pd.read_csv(assay, nrows=5)
        except pd.errors.ParserError as e:
            print(f"Parser error on {aid}")
            continue

        if df.empty:
            empty_counter = empty_counter + 1
    print(f"There are {empty_counter} assays with no data!")
    return empty_counter

def find_columns_nada(to_csv=False, unique=False) -> list:
    """ finds and prints the unique columns across the set of failed csvs

    to_csv: str, path to write the unique columns to, if False return the list
    """
    all_columns = []

    for aid_file in ALL_ASSAY_CSV:
        # just read the header so we dont
        # load all the data into memory
        # to speed things up
        try:
            assay_data_columns = pd.read_csv(aid_file, index_col=0, nrows=0)
            all_columns.extend(list(assay_data_columns.columns))
        except pd.errors.ParserError as e:
            print(f"Error on {aid_file}")

    all_columns = list(set(all_columns)) if unique else list(all_columns)

    if not to_csv:
        return all_columns

    f = open(to_csv, 'w')

    for column in all_columns:
        f.write(str(column) + '\n')
    f.close()


def move_concise_files():
    """ move concise files from local to the box """

    for pubchem_aid in DR_AIDS:
        aid_folder = find_folder(pubchem_aid)
        local_file_concise = os.path.join(config.Config.LOCAL_CONCISE, aid_folder, f'{pubchem_aid}.concise.csv')

        new_file = os.path.join(config.Config.BOX_PATH, 'DATA', 'Nada', 'concise_dr', f'{pubchem_aid}.concise.csv')
        shutil.copyfile(local_file_concise, new_file)


def move_json_files():
    """ move json files from local to the box """

    for pubchem_aid in DR_AIDS:
        aid_folder = find_folder(pubchem_aid)
        local_file_json = os.path.join(config.Config.LOCAL_JSON, aid_folder, f'{pubchem_aid}.json.gz')

        new_file = os.path.join(config.Config.BOX_PATH, 'DATA', 'Nada', 'json_dr', f'{pubchem_aid}.json')

        with gzip.open(local_file_json, 'rb') as f_in:
            with open(new_file, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

def build_sql_db(DB_FILE):
    """ build a sqlite database from all the converted CSV files

     """

    con = sql.connect(DB_FILE)

    TOTAL_ASSAYS = len(CONVERTED_CSV)
    COUNTER = 1

    # create a counter for id
    LAST_ID = 1

    FAILED_CSV_SQL = []

    # for the converted ones just
    # do as normal.  Not every dataframe
    # has data.  Need to check why
    for assay in CONVERTED_CSV:

        try:
            df = pd.read_csv(assay)
            if not df.empty:

                ids = list(range(LAST_ID, LAST_ID+df.shape[0]))
                df['ID'] = ids
                df.to_sql('dose_response', con=con, if_exists='append', index=False)
                LAST_ID = ids[-1]
            else:
                print(f"{assay} if empty!")
                FAILED_CSV_SQL.append(assay)
        except:
            FAILED_CSV_SQL.append(assay)


        print(f"{COUNTER / TOTAL_ASSAYS * 100}% Completed")
        COUNTER += 1


    error_file = os.path.join(os.path.dirname(DB_FILE), 'failed_aids.csv')

    with open(error_file, 'w') as f:
        for aid in FAILED_CSV_SQL:
            f.write(str(aid) + "\n")

def add_concise_sql(DB_FILE):
    """ add all the corresponding consise data files to the sqlite database """


    LAST_ID = 1
    ALL_COLUMNS = []

    # need to first read all the columns
    # inorder to create the table in sql
    # and make sure its frame has the same columns
    for aid_file in CONCISE_CSVS:
        print(aid_file)
        assay_data_columns = read_concise_aid(aid_file, just_columns=True)
        ALL_COLUMNS.extend(list(assay_data_columns))

    ALL_COLUMNS = list(set([col.upper() for col in ALL_COLUMNS]))

    with sql.connect(DB_FILE) as con:

        for aid_file in CONCISE_CSVS:

            aid = int(os.path.basename(aid_file).replace('.concise.csv', ''))
            print(aid)
            # add aid
            concise_data = read_concise_aid(aid_file)
            concise_data['PUBCHEM_AID'] = aid
            concise_data.columns = [col.upper() for col in concise_data.columns]
            if not concise_data.empty:
                ids = list(range(LAST_ID, LAST_ID + concise_data.shape[0]))
                concise_data['ID'] = ids

                for col in ALL_COLUMNS:
                    if col not in concise_data.columns:
                        concise_data[col] = np.nan

                concise_data.to_sql('concise', con=con, if_exists='append', index=False)

                LAST_ID = ids[-1]


def add_targets(DB_FILE):
    """ add all the target datat o the sql database """

    target_file = 'data/target_info.csv'
    target_data = pd.read_csv(target_file, index_col=0)

    with sql.connect(DB_FILE) as con:

        target_data.to_sql('targets', con=con, if_exists='replace', index=True)





def fit_assay_to_hill(df: pd.DataFrame) -> pd.DataFrame:
    """ takes assay information as a dataframe as will fit each SID to a hill curve and extract
    the apprioriate parameters """

    # infinity covariance matrix
    # https://stackoverflow.com/questions/27230285/numpy-polyfit-gives-useful-fit-but-infinite-covariance-matrix

    data = []

    for sid, sid_data in df.groupby('SID'):
        x = sid_data['Concentration'].values
        y = sid_data['Response'].values
        direction = sid_data['Direction'].values[0]
        if direction == 'Descending':
            y = y*-1


        init_params = np.array([1, y.max(), 1.0])

        try:
            # changing maxfev could possible solve
            # the number of NaNs?
            # seems like it was the slope that was
            # the parameter that was not being able to
            # be found....fixed with bounds
            fit_params, pcov = curve_fit(hill_curve,
                                         x,
                                         y,
                                         init_params,
                                         bounds=([-np.inf, min(0, y.max()), 0.3], [np.inf, 120, 8]), # constrains found in TCPL vignette
                                         maxfev=5000
                                         )
            model_preds = hill_curve(x, *fit_params)
            error = model_preds - y
            SE = np.square(error)  # squared errors
            MSE = np.mean(SE)  # mean squared errors
            RMSE = np.sqrt(MSE)  # Root Mean Squared Error, RMSE
            Rsquared = 1.0 - (np.var(error) / np.var(y))

            vals = [sid] + list(fit_params) + [SE, MSE, RMSE, Rsquared]

            # if the direction fit is negative, convert the top
            # to negative
            if direction == "Descending":
                vals[2] = -1*vals[2]

        except RuntimeError as e:
            # runtime error is usually
            # that the parameters for the fit
            # were not found
            vals = [sid, None, None, None, None, None, None, None]

        data.append(vals)

    return pd.DataFrame(data, columns=['SID', 'AC50', 'TOP', 'SLOPE', 'SE', 'MSE', 'RMSE', 'R^2'])


def check_direction(df: pd.DataFrame) -> bool:
    """  check the analysis direction to find if the chemical results in a positive or negastive response
     NOTE: does not account for complete loss of signal at high concentrations """
    data = df.sort_values('log(Concentration)')
    curve_direction = data.iloc[0].Response - data.iloc[-1].Response
    # curve direction is > 0
    # means descendning
    data['Direction'] = 'Ascending' if curve_direction <= 0 else 'Descending'


    return data

def plot_dosresponse(df, sid):
    df = df[df.SID == sid]
    df['Concentration_Log'] = np.log10(df['Concentration'].astype(float))
    df['Response'] = df['Response'].astype(float)

    ax = sns.scatterplot(data=df, x='Concentration_Log', y='Response')

    ac50, top, slope = df.AC50.iloc[0], df.TOP.iloc[0], df.SLOPE.iloc[0]

    print(ac50, top, slope)
    xs = np.linspace(df.Concentration_Log.min(), df.Concentration_Log.max(), 100)
    #xs = np.linspace(0.001, df.Concentration.max(), 100)
    ys = hill_curve(10**xs, ac50, top, slope)
    #ax.plot(xs, ys)
    ax.plot(xs, ys)

    plt.show()


def fit_curves(DB_FILE):
    """ fits tries to fit every unique assay available in the dose response table
    to a hill curve....it will only fit "ACTIVE" chemicals i.e., CIDs defined as active
    in the concise assays

    the curve fitting is processed in a similar manner to the TCPL package in R.
    the vignette can be found here: https://cran.r-project.org/web/packages/tcpl/vignettes/Data_processing.html#uploading-and-processing-data
     """
    # number of concentration
    # a chemical has to have
    # to be fit for a curve

    MINIMUM_CONC_POINTS = 4

    # check if table already exists,
    # and remove if it does, then
    # get all unique AIDs in the
    # dose_response database
    with sql.connect(DB_FILE) as con:

        result = con.execute('DROP TABLE IF EXISTS hill_models;').fetchone()
        print(result)
        TOTAL_ASSAYS = pd.read_sql_query('SELECT DISTINCT PUBCHEM_AID FROM concise', con=con).PUBCHEM_AID.values.tolist()


    COUNTER = 1

    # create a counter for id
    LAST_ID = 1

    FAILED_CSV_SQL = []




    for assay in TOTAL_ASSAYS:

        # get all active CIDS for the assay
        actives_query = 'SELECT c.PUBCHEM_CID as CID, c.PUBCHEM_AID as AID, c.PUBCHEM_SID as SID ' \
                        'FROM concise c ' \
                        'WHERE c.PUBCHEM_AID == {} AND c.PUBCHEM_ACTIVITY_OUTCOME == "Active" AND ' \
                        'c.PUBCHEM_CID is not null AND c.PUBCHEM_SID is not null'.format(assay)

        active_cmps = pd.read_sql_query(actives_query, con=config.Config.DB_URI)
        active_cmps['SID'] = active_cmps['SID'].astype(int)

        # now get the dose responses
        sid_list = [str(sid) for sid in active_cmps.SID]

        sid_string = ", ".join(map(str, sid_list))
        sid_query = f'({sid_string})'

        dr_query = 'SELECT  SID, AID, concentration as Concentration, response as Response ' \
                   'FROM dose_response ' \
                   'WHERE AID == {} AND SID in {} '.format(assay, sid_query)

        dose_responses = pd.read_sql_query(dr_query, con=config.Config.DB_URI)
        if dose_responses.empty:
            print(f"skipping {assay}")
            continue

        from pandas.api.types import is_numeric_dtype
        if not is_numeric_dtype(dose_responses['Response']):
            print(dose_responses)
            print(f"skipping {assay}")
            continue
        # remove compounds where there
        # are not enough data points
        dose_responses = dose_responses.groupby(['AID', 'SID']).filter(lambda sid: sid.Concentration.nunique() >= MINIMUM_CONC_POINTS)
        if dose_responses.empty:
            print(f"skipping {assay}")
            continue

        # avg responses per SID
        # TODO: could possibly find a way to do this w/o averaging?
        dose_responses['log(Concentration)'] = np.log10(dose_responses.Concentration)
        dose_responses['log(Concentration)'] = dose_responses['log(Concentration)'].round(2)

        # average multiple responses per concentration
        # then apply curveP to correct points
        response_avg = dose_responses.groupby(['AID', 'SID', 'log(Concentration)'])['Response'].mean().reset_index()
        response_avg_corrected = response_avg.groupby(['AID', 'SID']).apply(curveP).reset_index(drop=True)
        if response_avg_corrected.empty:
            break
        print(response_avg_corrected.head())
        data = response_avg_corrected.groupby(['AID', 'SID']).apply(check_direction).reset_index(drop=True)
        data['Concentration'] = 10**data['log(Concentration)']
        ac50s = fit_assay_to_hill(data)
        ac50s['AID'] = assay

        ids = list(range(LAST_ID, LAST_ID + ac50s.shape[0]))

        ac50s['ID'] = ids

        ac50s.to_sql('hill_models', con=con, if_exists='append', index=False)

        LAST_ID = ids[-1]








if __name__ == '__main__':

    #TARGET_AID = 1346983
    #TARGET_SID = 144205764
    # TARGET_AID = 411
    # TARGET_SID = 7977150
    # df = pd.read_csv(os.path.join(ANOTHER_FAILED_CSV_DIR, f'{TARGET_AID}.csv'), index_col=0)
    # df = df.query("SID == @TARGET_SID")
    # df['Response'] = df['Response']*-1
    # ac50s = fit_assay_to_hill(df)
    # df = df.merge(ac50s).dropna()
    # #df = df[df.TOP > 80]
    # #print(df.sort_values(['RMSE'], ascending=[True]).iloc[:10])
    # plot_dosresponse(df, TARGET_SID)
    #sql_file = os.path.join(config.Config.BOX_PATH, 'DATA', 'Nada', 'SQLite', )
    sql_file = os.path.join("/Volumes", "TOSHIBA EXT", "data", "nada", "SQLIte", "pubchem_dr.db")

    # ac50_dir = os.path.join(config.Config.BOX_PATH, 'DATA', 'Nada', 'ac50s')
    # fit_curves(ac50_dir)
    #build_sql_db(sql_file)
    #add_concise_sql(sql_file)
    fit_curves(sql_file)
    #targets(sql_file)