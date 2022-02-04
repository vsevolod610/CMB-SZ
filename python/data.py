""""
Extract data from data.csv

require:    data.py
make:       SZ_data.txt

coments:    Pay attention to spaces in rows
"""


import csv
import numpy as np


with open('data.csv', newline='') as File:
    reader = csv.reader(File, delimiter=';')
    data = [row for row in reader]
    data = np.array(data)
    data = data.astype(float)

nu_data = data[:, 0]
SZ_data = data[:, 1]
SZ_data_errors = data[:, [2, 3]].T