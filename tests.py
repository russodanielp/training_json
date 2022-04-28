import unittest
import os
import glob
import pandas as pd
import config
from curve_fitting import find_folder, CONCISE_DATA_DIR, TARGET_ASSAYS



class TestSum(unittest.TestCase):

    def test_find_folder(self):
        for aid in TARGET_ASSAYS:
            aid_folder = find_folder(aid)
            full_path = os.path.join(CONCISE_DATA_DIR, aid_folder, f'{aid}.concise.csv')
            self.assertTrue(os.path.exists(full_path))

if __name__ == '__main__':
    unittest.main()