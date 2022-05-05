import pandas as pd
import requests
from functools import reduce
import os
import sqlite3 as sql
import config

def get_targets(aid_list):
    """ function to get target information for a list of aids """
    # convert list of identifers to str
    aid_list = list(map(str, aid_list))

    # make the base URL for the PubChem POST Request
    url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/targets/ProteinGI,ProteinName,GeneID,GeneSymbol/json"

    headers = {'Content-Type': 'multipart/form-data'}
    data = {'aid': ','.join(aid_list)}

    response = requests.post(url, data=data)

    return response.json()

def load_targets_sql():
    """ get target information from pubchem and load it in the SQLite database """
    dr_aids = pd.read_table('data/dr_aids.txt', header=None, names=['AIDS'])
    target_json = get_targets(dr_aids['AIDS'].values)

    targets = pd.DataFrame(target_json['InformationList']['Information'])
    dataframes = []

    for col in ['GI', 'GeneID', 'ProteinName', 'GeneSymbol']:
        new_df = targets.set_index('AID')[col].explode().reset_index()
        dataframes.append(new_df)

    targets_clean = reduce(lambda x, y: pd.merge(x, y, on='AID'), dataframes)
    targets_clean = targets_clean[targets_clean[["GI", "GeneID", "ProteinName", "GeneSymbol"]].notnull().any(1)]

    sqlite_file = os.path.join(config.Config.BOX_PATH, 'DATA', 'Nada', 'SQLite', 'pubchem_dr.db')
    con = sql.connect(sqlite_file)

    targets_clean['ID'] = list(range(1, targets_clean.shape[0]+1))

    targets_clean.to_sql('targets', con=con, if_exists='replace', index=False)

if __name__ == '__main__':
    load_targets_sql()