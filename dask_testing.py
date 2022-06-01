import dask.dataframe as dd
import pandas as pd
from curve_fitting import CONCISE_DATA_DIR
import config, os, glob

CONCISE_CSVS = glob.glob(os.path.join(config.Config.BOX_PATH, 'DATA', 'Nada', 'concise_dr', '*.csv'))

IN_VIVO_DATA = pd.read_csv('../pc-bayes/new_data/er_data.csv').rename(columns={'CID': 'PUBCHEM_CID', 'AGONIST_CLASS': 'Tox'})
IN_VIVO_DATA['PUBCHEM_CID'] = IN_VIVO_DATA['PUBCHEM_CID'].astype('int')

TARGET_INFO = pd.read_csv('data/target_info.csv', index_col=0).rename(columns={'AID': 'PUBCHEM_AID'})



dtypes = {
            'PUBCHEM_ASSAYDATA_COMMENT': 'object',
            'PUBCHEM_ACTIVITY_URL': 'object'

            #'PUBCHEM_SID': 'int'

}

df = dd.read_csv(CONCISE_CSVS[:],
                 dtype=dtypes,
                 include_path_column='PUBCHEM_AID',
                 error_bad_lines=False)
df = df[df.PUBCHEM_RESULT_TAG.astype(str).str.isnumeric()]

df = df[df.PUBCHEM_CID.notnull()]

df['PUBCHEM_SID'] = df.PUBCHEM_SID.astype(int)
df['PUBCHEM_CID'] = df.PUBCHEM_CID.astype(int)

# store file name as AID
df['PUBCHEM_AID'] = df.PUBCHEM_AID.apply(lambda path: int(os.path.basename(path).split('.')[0]), meta=('PUBCHEM_AID', 'int'))

# merge tox and target
df = df.merge(IN_VIVO_DATA[['PUBCHEM_CID', 'Tox']], how='left', on='PUBCHEM_CID')
df = df.merge(TARGET_INFO[['PUBCHEM_AID', 'GeneSymbol']], how='left', on='PUBCHEM_AID')

df = df[df['Tox'].notnull()]
df = df[df.GeneSymbol.notnull()]
df = df[df['PUBCHEM_ACTIVITY_OUTCOME'].notnull()]
df['PUBCHEM_ACTIVITY_OUTCOME'] = df['PUBCHEM_ACTIVITY_OUTCOME'].astype(str)

cid_counts_active = df[['PUBCHEM_CID', 'GeneSymbol', 'PUBCHEM_ACTIVITY_OUTCOME']].groupby(['PUBCHEM_CID', 'GeneSymbol']).apply(lambda data: (data['PUBCHEM_ACTIVITY_OUTCOME'] == 'Active').sum(), ).compute()
#cid_counts_inactive = df[['PUBCHEM_CID', 'GeneSymbol', 'PUBCHEM_ACTIVITY_OUTCOME']].groupby(['PUBCHEM_CID', 'GeneSymbol']).apply(lambda data: (data['PUBCHEM_ACTIVITY_OUTCOME'] == 'Inactive').sum(), ).compute()
tested_targets = df[['PUBCHEM_CID', 'GeneSymbol', 'PUBCHEM_AID']].groupby(['PUBCHEM_CID', 'GeneSymbol'])['PUBCHEM_AID'].nunique().compute()

# df = pd.concat([cid_counts_active, cid_counts_inactive, tested_targets], axis=1)
# df.columns = ['Active', 'Inactive', 'N_AIDS']

df = pd.concat([cid_counts_active, tested_targets], axis=1)
df.columns = ['Active', 'N_AIDS']

df = df.reset_index().merge(IN_VIVO_DATA[['PUBCHEM_CID', 'Tox']])

df.to_csv('data/for_bayes.csv')

# aids = []
# act_toxs = []
# aid_groups = df[['PUBCHEM_AID', 'PUBCHEM_ACTIVITY_OUTCOME', 'Tox']].groupby('PUBCHEM_AID')
#
# act_tox_func = lambda data: ((data.PUBCHEM_ACTIVITY_OUTCOME == 'Active') & (data.Tox == 1)).sum()
# inact_tox_func = lambda data: ((data.PUBCHEM_ACTIVITY_OUTCOME == 'Inactive') & (data.Tox == 1)).sum()
# act_nontox_func = lambda data: ((data.PUBCHEM_ACTIVITY_OUTCOME == 'Active') & (data.Tox == 0)).sum()
# inact_nontox_func = lambda data: ((data.PUBCHEM_ACTIVITY_OUTCOME == 'Inctive') & (data.Tox == 0)).sum()
#
# act_toxs = aid_groups.apply(act_tox_func).compute()
# inact_toxs = aid_groups.apply(inact_tox_func).compute()
# act_nontoxs = aid_groups.apply(act_nontox_func).compute()
# inact_nontoxs= aid_groups.apply(inact_tox_func).compute()
#
# confusion_matrix = pd.concat([act_toxs, inact_toxs, act_nontoxs, inact_nontoxs], axis=1)
# confusion_matrix.columns = ['ActTox', 'InactTox', 'ActNonTox', 'InActNonTox']
# print(confusion_matrix.sort_values('ActTox'))