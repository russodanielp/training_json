import dask.dataframe as dd
import pandas as pd
from curve_fitting import CONCISE_DATA_DIR
import config, os, glob

a50_files = glob.glob(os.path.join(config.Config.BOX_PATH, 'DATA', 'Nada', 'ac50s', '*.csv'))


# dtypes = {
#             'PUBCHEM_ASSAYDATA_COMMENT': 'object',
#             'PUBCHEM_ACTIVITY_URL': 'object'
#
#             #'PUBCHEM_SID': 'int'
#
# }

df = dd.read_csv(a50_files[:],
                 #dtype=dtypes,
                 include_path_column='PUBCHEM_AID',
                 error_bad_lines=False)


# store file name as AID
df['PUBCHEM_AID'] = df.PUBCHEM_AID.apply(lambda path: int(os.path.basename(path).split('.')[0].split('_')[0]), meta=('PUBCHEM_AID', 'int'))


lowest_rmse = df.RMSE.min().compute()
highest_rsme = df.RMSE.max().compute()
print(lowest_rmse, highest_rsme)