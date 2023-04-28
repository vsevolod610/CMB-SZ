# -*- coding: utf-8 -*-
"""
Show: generate data
"""

# setting path
import sys
sys.path.append('../core')


import numpy as np
import matplotlib.pyplot as plt

from data import SZ_data_read
from show_data import data_show
from start_analyze import N


path_pic_fit = '../data/set/pic_fit/pic_fit-{}'
path_to_SZ_data = '../data/set/szdata/szdata{}.txt'


if __name__ == "__main__":
    for i in range(N):
        print('step = {:<4}'.format(i+1))
        data_show(path_to_SZ_data.format(i + 1), save=path_pic_fit.format(i + 1))
