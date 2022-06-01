import pandas as pd
import config

mmp_aid = 720635

query = 'SELECT * from CONCISE where PUBCHEM_AID == {}'.format(mmp_aid)

q2 = "SELECT PUBCHEM_AID, PUBCHEM_SID, COUNT(DISTINCT PUBCHEM_ACTIVITY_OUTCOME) " \
    "FROM concise " \
"GROUP BY PUBCHEM_AID, PUBCHEM_SID "

data = pd.read_sql_query(q2, config.Config.DB_URI)



print(data)

# max_responses = data.groupby(['PUBCHEM_CID']).count()
# print(max_responses)
#
#
