# -*- coding: utf-8 -*-
"""
Analyze example cluster
"""

# setting path
import sys
sys.path.append('../core')


import numpy as np
import matplotlib.pyplot as plt

from sz_analyze import sz_analyze as SZ_mcmc

#method = 'T0'
#method = 'lazy'
method = 'simple'
#method = 'Tz'

# data
path_SZ_data = '../data/one/SZ_data.txt'
path_startSZ = '../data/startSZ.sss'
path_prior = '..//data/one/prior.dat'

path_pic_chain = '../data/one/pic_chain'
path_pic_consum = '../data/one/pic_consum'
path_pic_fit = '../data/one/pic_fit'

paths = [path_SZ_data, path_startSZ, path_prior]
path_pics = [path_pic_chain, path_pic_consum, path_pic_fit]
z, *T0_params = SZ_mcmc(method, paths, pics=path_pics)

print(*T0_params)

plt.show()
