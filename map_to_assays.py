"""
Find top targets from posterior
Gather assays from top targets
Gather dose responses for active compounds
"""

import pandas as pd
import config
import seaborn as sns
import matplotlib.pyplot as plt
import sqlite3 as sql
import dask.dataframe as dd

# take top N targets for modeling
N_TARGETS = 1

p_ppv = pd.read_csv('data/p_ppv.csv', index_col=0)
p_sens = pd.read_csv('data/p_sens.csv', index_col=0)

tox_rel = (p_ppv + p_sens) / 2

print(tox_rel.columns)

top_targets_ = tox_rel.quantile(0.05).sort_values(ascending=False).iloc[:N_TARGETS].index.tolist()

# add strings
top_targets = ["\"" + t + "\"" for t in top_targets_]

targets = ", ".join(map(str, top_targets))
targets_array = f'({targets})'

query = 'SELECT tG.PUBCHEM_AID, tG.GeneSymbol ' \
        'FROM targets tG ' \
        'WHERE GeneSymbol IN {}' \
        ''.format(targets_array)


query = 'SELECT DISTINCT tG.PUBCHEM_AID, tG.GeneSymbol, c.PUBCHEM_CID as CID, dR.SID, dR.Concentration, dR.Response ' \
        'FROM targets tG ' \
        'INNER JOIN dose_response dR on dR.AID == tG.PUBCHEM_AID ' \
        'INNER JOIN concise c on c.PUBCHEM_SID == dR.SID ' \
        'WHERE GeneSymbol IN {} limit 1000'.format(targets_array)

genes_aids = dd.read_sql_table(query, config.Config.DB_URI)

# cids_per_target = (
#                   genes_aids
#                   .filter(['GeneSymbol', 'CID', 'SID'])
#                   .groupby(['GeneSymbol', 'PUBCHEM_CID'])
#                   .count()
#                   .sort_values('PUBCHEM_AID', ascending=False)
#                   )
#
print(genes_aids)