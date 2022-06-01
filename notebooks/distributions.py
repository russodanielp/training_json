
#%%
# imports
import pandas as pd
import dask.dataframe as dd
from dask import delayed
import config
import seaborn as sns
import dask.array as da
import matplotlib.pyplot as plt
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

data = data[data.Concentration > 0]
data['Concentration'] = da.log10(data.Concentration)

h, bins = da.histogram(data.Concentration, bins=100, range=[data.Concentration.min(), data.Concentration.max()])

bins = bins.compute()
h = h.compute()

#data = pd.read_sql_query(query, con=config.Config.DB_URI, chunksize=10000)


plt.bar(bins[:-1], h, width=1)

plt.show()
#
# dataframes = [delayed(d)(filename) for filename in data]
# df = dd.from_delayed(dataframes)

#%%
import numpy as np

# #%%
# sns.histplot(data=data, x='conc', log_scale=True, bins=100)
# plt.show()


