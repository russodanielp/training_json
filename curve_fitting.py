"""
module for curve fitting assays from a pubchem

this stack overflow has a demo on how to curve fit a function: https://stackoverflow.com/questions/55078451/how-to-use-curvefit-in-python

 """
import os
import config
import glob
import pandas as pd


CONCISE_DATA_DIR = os.path.join(config.Config.BOX_PATH, 'DATA', 'pubchem', 'bioassay', 'concise_csv', 'all')
ASSAY_FOLDERS = glob.glob(os.path.join(CONCISE_DATA_DIR, '*'))
TARGET_ASSAYS = sorted(pd.read_table(os.path.join('data', 'dr_aids.txt'), header=None)[0].values.tolist())


def hill_curve(conc: float, ac50: float, top: float, slope: float) -> float:
    """ returns """
    denom = conc**slope + ac50**slope
    return top * (conc**slope / denom)


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


def read_aid(pubchem_aid: int) -> pd.DataFrame:
    """ given an aid, will find the appropriate folder and read the csv """

    aid_folder = find_folder(pubchem_aid)
    full_path = os.path.join(CONCISE_DATA_DIR, aid_folder, f'{pubchem_aid}.concise.csv')

    # DELETE ME ONLY USED FOR TESTING
    if not os.path.exists(full_path):
        return pd.DataFrame()

    # the first n rows are a header
    # that describes the data
    # dose response example is 434931
    assay_results = pd.read_csv(full_path, error_bad_lines=False)

    # find index where header stops
    idx = assay_results[assay_results.PUBCHEM_RESULT_TAG == '1'].index[0]
    assay_results = assay_results.loc[idx:]

    return assay_results


def find_unique_columns():

    all_columns = []

    for aid in TARGET_ASSAYS:
        assay_data = read_aid(aid)
        if not assay_data.empty:
            all_columns.extend(list(assay_data.columns))

    return list(set(all_columns))

if __name__ == '__main__':
    print(find_unique_columns())


