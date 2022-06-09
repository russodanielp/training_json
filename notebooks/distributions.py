
#%%
# imports
import pandas as pd
import dask.dataframe as dd
from dask import delayed
import config
import seaborn as sns
import dask.array as da
import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LinearRegression
sns.set_style('dark')

#%%

save="string"
print(save)
print("Testing..Working")

#%%

print(save)
#%%

add_limit = False

query = "select AID as aid, Concentration as conc from dose_response " \
        "WHERE conc != 0 "
limit = "limit 1000000 "

if add_limit:
    query = query + limit

data = dd.read_sql_table('dose_response', config.Config.DB_URI, index_col='ID', columns=['Concentration', 'ID'])
#data = dd.from_pandas(pd.DataFrame(np.random.rand(10000,), columns=['Concentration']), npartitions=2)


data = data[data.Concentration > 0]
data['Concentration'] = da.log10(data.Concentration)

h, bins = da.histogram(data.Concentration, bins=50, range=[data.Concentration.min(), data.Concentration.max()])

bins = bins.compute()
h = h.compute()

#data = pd.read_sql_query(query, con=config.Config.DB_URI, chunksize=10000)
sns.set_theme(style="whitegrid")
fig, ax = plt.subplots()

ax.bar(range(len(bins)-1), h, width=1)

ax.set_xlabel("Log(Concentration)")

# for scaling the axis

lm = LinearRegression()

#lm.fit(list(range(len(bins))), bins)

m, b = np.polyfit(list(range(len(bins))), bins, 1)

ax.set_xticklabels(['{:.2f}'.format(x) for x in ax.get_xticks()*m + b])
print(bins)

plt.show()
#
# dataframes = [delayed(d)(filename) for filename in data]
# df = dd.from_delayed(dataframes)

#%%
import numpy as np

# #%%
# sns.histplot(data=data, x='conc', log_scale=True, bins=100)
# plt.show()


