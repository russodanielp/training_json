import json
import gzip
import shutil
import pandas as pd
import os


f = open('failed_aids.txt')
failed_aids = []

for line in f:
    line.strip()
    failed_aids.append(int(line))

f.close()

for aid in failed_aids:

    json_file_csv = r'E:\JSON\Bioassays\failed_aids\{}.csv'.format(aid)

    if os.path.exists(json_file_csv):
        print('Skipped AID: {}'.format(aid))
        continue

    json_file_in = r'E:\JSON\Bioassays\failed_aids\{}.json.gz'.format(aid)
    json_file_out = r'E:\JSON\Bioassays\failed_aids\{}.json'.format(aid)

    with gzip.open(json_file_in, 'rb') as f_in:
        with open(json_file_out, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

    f = open(json_file_out)
    json_file = json.load(f)
    f.close()

    results = json_file['PC_AssaySubmit']['assay']['descr']['results']
    data = json_file['PC_AssaySubmit']['data']

    headers = {}
    headers_list = []
    tid_list = []
    dictionary = {}

    for tid in results:
        for key in tid:
            if (key == 'name') and (tid['type'] == 1):
                headers[tid['tid']] = tid.get('name')
                tid_list.append(tid['tid'])
                headers_list.append(tid['name'])

    for i in range(len(data)):
        for key, value in data[i].items():
            sid = data[i]['sid']
            dictionary[sid] = {}
            for j in range(len(data[i]['data'])):
                tid_num = data[i]['data'][j]['tid']
                for num in tid_list:
                    if tid_num == num:
                        dictionary[sid][num] = data[i]['data'][j]['value']['fval']

    df = pd.DataFrame(dictionary).T

    for key, value in headers.items():
        df.rename(columns={key: value}, inplace=True)

    df.to_csv(r'E:\JSON\Bioassays\failed_aids\{}.csv'.format(aid))
    print(aid)


# json_file = r'E:\JSON\Bioassays\failed_aids'
#
# with gzip.open('2823.json.gz', 'rb') as f_in:
#     with open('2823.json', 'wb') as f_out:
#         shutil.copyfileobj(f_in, f_out)

