# -*- coding: utf-8 -*-
"""
Analyze example cluster
"""

# Imports
import numpy as np
import matplotlib.pyplot as plt

import paths # manage path imports by paths.py
from core_sz.sz_mcmc import sz_analyze
from core_sz.model import model_100 as model


# Paths
data_dir = paths.src + '/../data/example_cluster/'
path_SZdata = data_dir + 'SZ_data.txt'
path_prior = data_dir + 'prior.dat'
path_startSZ = data_dir + 'startSZ.sss'

path_pic_chain = 'pics/mcmc/pic_chain'
path_pic_params = 'pics/mcmc/pic_params'
path_pic_fit = 'pics/mcmc/pic_fit'

paths_input = [path_SZdata, path_prior, path_startSZ]
path_pics = [path_pic_chain, path_pic_params, path_pic_fit]

# Mcmc
z, *T0_params = sz_analyze(model, paths_input, pics=path_pics)
#print(*T0_params)
#plt.show()
