import json
import gzip
import shutil
import pandas as pd
import os
import json
import config
import glob
import ntpath

def parse_json(json_file: str) -> pd.DataFrame:
    """ parse a json file and convert to a dataframe in a similar fashion as the ones returned by PugREST
    the columns should be AID SID Concentration	Concentration Unit Response	Response Unit"""

    if json_file.split('.')[-1] == 'gz':
        json_file_out = json_file.replace('.gz', '')
        with gzip.open(json_file, 'rb') as f_in:
            with open(json_file_out, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        json_file = json_file_out


    f = open(json_file)
    json_data = json.load(f)
    f.close()


    # json files have two parts, results and
    # descriptions
    results = json_data['PC_AssaySubmit']['assay']['descr']['results']
    data = json_data['PC_AssaySubmit']['data']

    # it looks like respones that have a "tc" for test concentration
    # property are dose-response activities

    # dict that contains tids that are
    # concentration-activity points
    # stores unit information as well
    conc_tids = {}

    for tid in results:
        if 'tc' in tid.keys():

            tid_id = tid["tid"]
            conc_tids[tid_id] = [tid["unit"], tid["tc"]["concentration"], tid["tc"]["unit"]]

    result_dic = {}

    for compound in data:
        sid = compound["sid"]

        # remove points not in the conc_tids dictionary
        cmp_responses = filter(lambda datum: datum["tid"] in conc_tids.keys(), compound["data"])
        cmp_data = {}
        for response in cmp_responses:
            cmp_data[response["tid"]] = response["value"]["fval"]
        result_dic[sid] = cmp_data

    responses = pd.DataFrame(result_dic)
    responses.index.name = "tid"
    responses.columns.name = "SID"

    responses_melt = responses.T.reset_index().melt(id_vars=["SID"]).rename(columns={"value": "Response"}).dropna()

    # create a dataframe that has the response meta information
    # like concentration units, etc
    responses_meta = pd.DataFrame(conc_tids, index=["Response Unit", "Concentration", "Concentration Unit"]).T
    responses_meta.index.name = "tid"
    responses_meta = responses_meta.reset_index()

    final_frame = responses_melt.merge(responses_meta, on='tid', how='left')

    aid = json_data['PC_AssaySubmit']['assay']['descr']['aid']["id"]
    final_frame["AID"] = aid
    return final_frame[["AID", "SID", "Concentration", "Concentration Unit", "Response", "Response Unit"]].sort_values("SID")

def old():
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


def convert_json(target_dir):
    """ convery the assays from failed sids into csvs like return by pubchem """
    JSON_FILES = glob.glob(os.path.join(config.Config.BOX_PATH, 'DATA', 'Nada', 'failed_aids', '*.json'))

    for json_file in JSON_FILES:
        base_filename = ntpath.basename(json_file).split('.')[0]
        csv_file = os.path.join(target_dir, '{}.csv'.format(base_filename))
        if not os.path.exists(csv_file):
            df = parse_json(json_file)
            df.to_csv(csv_file)


if __name__ == '__main__':
    convert_json(os.path.join(config.Config.BOX_PATH, 'DATA', 'Nada', 'another_failed_aids'))