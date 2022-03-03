import json
import os, urllib.request

json_file = r"E:\pubchem\bioassay\Concise\JSON\all\0000001_0001000\1.json"


f = open(json_file)

data = json.load(f)

f.close()

assay_descr_info = data["PC_AssaySubmit"]["assay"]["descr"]

aid_info = {
    "name": assay_descr_info["name"],
    "aid": assay_descr_info["aid"]["id"],
    "description": assay_descr_info["description"]

}

f = open('{}_parsed.csv'.format(aid_info['aid']), 'w')

for key, val in aid_info.items():
    f.write(key + "," + str(val) + "\n")

f.close()

###### PugRest #####

target_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/1/description/json"

result = urllib.request.urlopen(target_url)
json_data = json.loads(result.read().decode())

print(json_data)