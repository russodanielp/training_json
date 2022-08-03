import pandas as pd
import config
import seaborn as sns
import matplotlib.pyplot as plt

# take top N targets for modeling
N_TARGETS = 20

p_ppv = pd.read_csv('data/p_ppv.csv', index_col=0)
p_sens = pd.read_csv('data/p_sens.csv', index_col=0)

tox_rel = (p_ppv + p_sens) / 2

print(tox_rel.columns)

top_targets_ = tox_rel.quantile(0.05).sort_values(ascending=False).iloc[:N_TARGETS].index.tolist()

# add strings
top_targets = ["\"" + t + "\"" for t in top_targets_]

targets = ", ".join(map(str, top_targets))
targets_array = f'({targets})'

query = 'SELECT PUBCHEM_AID, GeneSymbol ' \
        'FROM targets ' \
        'WHERE GeneSymbol IN {}'.format(targets_array)

genes_aids = pd.read_sql_query(query, config.Config.DB_URI)

counts = (
          genes_aids
          .groupby('GeneSymbol')
          .count()
          .sort_values('PUBCHEM_AID', ascending=False)
          )
sns.set_context("talk")
sns.set_style("ticks")
sns.boxplot(data=tox_rel[top_targets_], color="grey")

plt.ylabel("Toxicity Relevance\n[(PPV + Sens) /2]")
locs, labels = plt.xticks()
plt.setp(labels, rotation=45, ha="right")
plt.show(transparent=True)



sns.barplot(data=counts.loc[top_targets_].reset_index(),
            x='GeneSymbol', y='PUBCHEM_AID',
            facecolor="grey", edgecolor="black")
locs, labels = plt.xticks()
plt.setp(labels, rotation=45, ha="right")

plt.ylabel("# PubChem DR Assays")
plt.xlabel("")
plt.show(transparent=True)