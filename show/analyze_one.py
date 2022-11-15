# -*- coding: utf-8 -*-
"""
Analyze example cluster
"""
import sys

# setting path
sys.path.append('../core')


import numpy as np
import matplotlib.pyplot as plt

from mcmc_analyze import mcmc_analyze

# data
path_SZ_data = '../data/SZ_data.txt'
path_startSZ = '../data/startSZ.sss'
path_prior = '../data/prior.dat'

path_pic_chain = '../pic_chain'
path_pic_consum = '../pic_consum'
path_pic_fit = '../pic_fit'

paths = [path_SZ_data, path_startSZ, path_prior]
nwalkers, nsteps = 200, 200
path_pics = [path_pic_chain, path_pic_consum, path_pic_fit]
z, T0_fit, T0_sp, T0_sm = mcmc_analyze(paths, nwalkers, nsteps, pic=path_pics)

print(T0_fit, T0_sp, T0_sm)

plt.show()
