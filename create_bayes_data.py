import pandas as pd
import numpy as np
import config, os




# this is just dummy data for testing
# for now

IN_VIVO_DATA = pd.read_csv('../pc-bayes/new_data/er_data.csv').rename(columns={'CID': 'PUBCHEM_CID', 'AGONIST_CLASS': 'TOXICITY'})
IN_VIVO_DATA['PUBCHEM_CID'] = IN_VIVO_DATA['PUBCHEM_CID'].astype('int')


# get data from the SQL database

targets = pd.read_sql_table('targets', config.Config.DB_URI, index_col='index', columns=['PUBCHEM_AID', 'GeneSymbol'])

cids = IN_VIVO_DATA['PUBCHEM_CID'].sort_values().values.tolist()
cids = ", ".join(map(str, cids))
cid_array = f'({cids})'

query = 'SELECT PUBCHEM_CID, PUBCHEM_AID, PUBCHEM_ACTIVITY_OUTCOME ' \
        'FROM concise ' \
        'WHERE PUBCHEM_CID IN {}'.format(cid_array)

assay_activities = pd.read_sql_query(query, config.Config.DB_URI)

# one AID->CID assay pair
# can result in several bioactivity outcomes
# this can be due to several reasons (e.g.,
# multiple tests for an AID, or the CID is
# from several SIDs, like in mixtures).
# for now, just map the CIDs to one response
# giving preference to actives

# convert assays responses to be 1, NaN, or -1
assay_activities['PUBCHEM_ACTIVITY_OUTCOME_INT'] = (
                                                        assay_activities['PUBCHEM_ACTIVITY_OUTCOME']
                                                        .replace('Inactive', 0)
                                                        .replace('Active', 1)
                                                        .replace('Probe', 1)
                                                        .replace('Inconclusive', np.nan)
                                                        .replace('Unspecified', np.nan)
                                                    )

# now give preference to active
# responses.
assay_responses_reduced = (
                           assay_activities
                           .groupby(['PUBCHEM_AID', 'PUBCHEM_CID'])
                           ['PUBCHEM_ACTIVITY_OUTCOME_INT'].max()
                           .reset_index()
                           .rename(columns={'PUBCHEM_ACTIVITY_OUTCOME_INT': 'OUTCOME'})
                           )


# add target info and IN VIVO data
data = (
        assay_responses_reduced
        .merge(targets, how='left', on='PUBCHEM_AID')
        .merge(IN_VIVO_DATA[['PUBCHEM_CID', 'TOXICITY']], on='PUBCHEM_CID')
        )


# now convert the data for bayes.
# basically we want to see if a
# chemical is active for a target

assays_per_gene = (
                    data
                    .groupby(['GeneSymbol'])
                    ['PUBCHEM_AID'].nunique()
                    .sort_values(ascending=False)
                   )

gene_responses_per_cmp = (
                            data
                            .groupby(['PUBCHEM_CID', 'GeneSymbol'])
                            ['OUTCOME'].mean()
                            .reset_index()
                            .rename(columns={'OUTCOME': 'OUTCOME_MEAN'})
                            .sort_values('OUTCOME_MEAN', ascending=False)
                            .dropna(subset=['OUTCOME_MEAN'])
                            .assign(GENE_OUTCOME = lambda x: (x['OUTCOME_MEAN'] > 0.1).astype(int)) # cmp is active if > 10%
                            .merge(IN_VIVO_DATA[['PUBCHEM_CID', 'TOXICITY']], how='left', on='PUBCHEM_CID')
                         )


# confusion matrix time

def confusion_matrix(df):

    tps = ((df.TOXICITY == 1) & (df.GENE_OUTCOME == 1)).sum()
    fps = ((df.TOXICITY == 0) & (df.GENE_OUTCOME == 1)).sum()
    tns = ((df.TOXICITY == 0) & (df.GENE_OUTCOME == 0)).sum()
    fns = ((df.TOXICITY == 1) & (df.GENE_OUTCOME == 0)).sum()
    return tps, fps, tns, fns


for_bayes = (
             gene_responses_per_cmp
             .groupby('GeneSymbol')
             .apply(confusion_matrix)
             .apply(pd.Series)
             .set_axis(['TP', 'FP', 'TN', 'FN'], axis=1)
             .reset_index()
             .sort_values('TP', ascending=False)
            )

for_bayes.to_csv('data/for_bayes.csv', index=False)

