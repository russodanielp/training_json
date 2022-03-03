import os
import glob
from zipfile import ZipFile
import gzip
import shutil

# set as your directory
BIO_ASSAY_DATA_DIR = "E:\\pubchem\\bioassay\\Concise\\CSV\\Data"

if not os.path.join(BIO_ASSAY_DATA_DIR, 'all'):
    os.mkdir(os.path.join(BIO_ASSAY_DATA_DIR, 'all'))

bioassay_data_files = glob.glob(os.path.join(BIO_ASSAY_DATA_DIR, '*.zip'))

for f in bioassay_data_files:

    with ZipFile(f, 'r') as zipObj:
       # Extract all the contents of zip file in current directory
       zipObj.extractall(os.path.join(BIO_ASSAY_DATA_DIR, 'all'))


BIOASSAY_JSON_FILES = glob.glob(os.path.join(BIO_ASSAY_DATA_DIR, "all", "*", "*.json.gz"))

for f in BIOASSAY_JSON_FILES:
    with gzip.open(f, 'rb') as f_in:
        with open(f.replace('.gz', ''), 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)