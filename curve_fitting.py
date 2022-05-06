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


CONCISE_DATA_DIR = os.path.join(config.Config.BOX_PATH, 'DATA', 'pubchem', 'bioassay', 'concise_csv', 'all')
FAILED_CSV_DIR = os.path.join(config.Config.BOX_PATH, 'DATA', 'Nada', 'failed_aids')
ANOTHER_FAILED_CSV_DIR = os.path.join(config.Config.BOX_PATH, 'DATA', 'Nada', 'another_failed_aids')
PASSED_CSV_DIR = os.path.join(config.Config.BOX_PATH, 'DATA', 'Nada', 'AID')
ASSAY_FOLDERS = glob.glob(os.path.join(CONCISE_DATA_DIR, '*'))
TARGET_ASSAYS = sorted(pd.read_table(os.path.join('data', 'dr_aids.txt'), header=None)[0].values.tolist())
FAILED_CSV = glob.glob(os.path.join(FAILED_CSV_DIR, '*.csv'))
ANOTHER_FAILED_CSV = glob.glob(os.path.join(ANOTHER_FAILED_CSV_DIR, '*.csv'))
PASSED_CSV = glob.glob(os.path.join(PASSED_CSV_DIR, '*.csv'))
ALL_ASSAY_CSV = ANOTHER_FAILED_CSV + PASSED_CSV

DR_AIDS = pd.read_table('data/dr_aids.txt', header=None, names=['AIDS'])['AIDS'].values.tolist()

def hill_curve(conc: float, ac50: float, top: float, slope: float) -> float:
    """ returns """
    denom = np.power(conc, slope) + np.power(ac50, slope)
    return top * (np.power(conc, slope) / denom)


def find_folder(pubchem_aid: int) -> str:
    """ PubChem assays are stored in folders containing 1000 assays.  E.g.,
     0000001_0001000, 0001001_0002000, etc.  This function takes a PubChem AID as input
     and finds the folder that aid would fall in.  """

    # because the folders are stored in 1000 increments,
    # taking the modulus of 1000 gives the remainder needed
    # to round down.  Then, add one.
    lower_bound = pubchem_aid - (pubchem_aid % 1000) + 1
    upper_bound = lower_bound + 1000 - 1

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


def find_unique_columns() -> list:
    """ finds and prints the unique columns across the set of concise csv """
    all_columns = []

    for aid in TARGET_ASSAYS:
        assay_data = read_aid(aid)
        if not assay_data.empty:
            all_columns.extend(list(assay_data.columns))
            print(aid, list(assay_data.columns))

    return list(set(all_columns))


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

        new_dir = os.path.join(config.Config.BOX_PATH, 'DATA', 'Nada', 'concise_dr')



def build_sql_db(DB_FILE):
    """ build a sqlite database from all the converted CSV files

     """

    con = sql.connect(DB_FILE)

    TOTAL_ASSAYS = len(PASSED_CSV) + len(ANOTHER_FAILED_CSV)
    COUNTER = 1

    # create a counter for id
    LAST_ID = 1

    FAILED_CSV_SQL = []

    # these files come from PUG-REST
    # the columns are not correct, I
    # reported the bug to the Pug-REST people
    # so we need to switch the concentration
    # and response units
    for assay in PASSED_CSV:

        try:
            df = pd.read_csv(assay, index_col=0).rename(columns={'Concentration Unit': 'Response Unit',
                                                                 'Response Unit': 'Concentration Unit'})
            ids = list(range(LAST_ID, LAST_ID+df.shape[0]))
            df['ID'] = ids
            df.to_sql('dose_response', con=con, if_exists='append', index=False)
        except:
            FAILED_CSV_SQL.append(assay)

        print(f"{COUNTER / TOTAL_ASSAYS * 100}% Completed")
        COUNTER += 1
        LAST_ID = ids[-1]
    # for the converted ones just
    # do as normal.  Not every dataframe
    # has data.  Need to check why
    for assay in ANOTHER_FAILED_CSV:

        try:
            df = pd.read_csv(assay)
            if not df.empty:

                ids = list(range(LAST_ID, LAST_ID+df.shape[0]))
                df['ID'] = ids
                df.to_sql('dose_response', con=con, if_exists='append', index=False)
        except:
            FAILED_CSV_SQL.append(assay)

        print(f"{COUNTER / TOTAL_ASSAYS * 100}% Completed")
        COUNTER += 1
        LAST_ID = ids[-1]


def add_concise_sql():
    """ add all the corresponding consise data files to the sqlite database """

    sqlite_file = os.path.join(config.Config.BOX_PATH, 'DATA', 'Nada', 'SQLite', 'pubchem_dr.db')

    LAST_ID = 1
    ALL_COLUMNS = []

    # need to first read all the columns
    # inorder to create the table in sql
    # and make sure its frame has the same columns
    for aid in DR_AIDS:
        assay_data_columns = read_aid(aid, just_columns=True)
        ALL_COLUMNS.extend(list(assay_data_columns))

    ALL_COLUMNS = list(set(ALL_COLUMNS))

    with sql.connect(sqlite_file) as con:

        for aid in DR_AIDS:

            concise_data = read_aid(aid)

            if not concise_data.empty:
                ids = list(range(LAST_ID, LAST_ID + concise_data.shape[0]))
                concise_data['ID'] = ids

                for col in ALL_COLUMNS:
                    if col not in concise_data.columns:
                        concise_data[col] = np.nan

                concise_data.to_sql('concise', con=con, if_exists='append', index=False)

                LAST_ID = ids[-1]


def fit_assay_to_hill(df: pd.DataFrame) -> pd.DataFrame:
    """ takes assay information as a dataframe as will fit each SID to a hill curve and extract
    the apprioriate parameters """

    # infinity covariance matrix
    # https://stackoverflow.com/questions/27230285/numpy-polyfit-gives-useful-fit-but-infinite-covariance-matrix

    data = []

    for sid, sid_data in df.groupby('SID'):
        x = sid_data['Concentration'].values
        y = sid_data['Response'].values
        init_params = np.array([1, 1, 1.0])

        try:
            fit_params, pcov = curve_fit(hill_curve, x, y, init_params)
            model_preds = hill_curve(x, *fit_params)
            error = model_preds - y
            SE = np.square(error)  # squared errors
            MSE = np.mean(SE)  # mean squared errors
            RMSE = np.sqrt(MSE)  # Root Mean Squared Error, RMSE
            Rsquared = 1.0 - (np.var(error) / np.var(y))

            vals = [sid] + list(fit_params) + [SE, MSE, RMSE, Rsquared]

        except RuntimeError as e:
            vals = [sid, None, None, None, None, None, None, None]

        data.append(vals)

    return pd.DataFrame(data, columns=['SID', 'AC50', 'TOP', 'SLOPE', 'SE', 'MSE', 'RMSE', 'R^2'])





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

    add_concise_sql()