#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Show: generate data
"""

# Imports
import os

import paths # manage path imports by paths.py
from show_data import data_show


path_SZdata = paths.src + '/main/input/szdata/szdata{}.txt'
path_pic = paths.src + '/tests/pics/gen/pic_data-{}.png'


dir_SZdata = os.path.dirname(path_SZdata)
N = len(os.listdir(dir_SZdata)) - 1 # without .gitkeep
for i in range(1, N + 1):
    if os.path.exists(path_SZdata.format(i)):
        print(f"step = {i:<4}")
        data_show(path_SZdata.format(i), save=path_pic.format(i))
