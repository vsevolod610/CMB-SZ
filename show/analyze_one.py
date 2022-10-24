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
path_to_SZ_data = '../data/SZ_data.txt'
path_to_startSZ = '../data/startSZ.sss'
path_to_prior = '../data/prior.dat'

paths = [path_to_SZ_data, path_to_startSZ, path_to_prior]
nwalkers, nsteps = 200, 200
z, T0_fit, T0_sp, T0_sm = mcmc_analyze(paths, nwalkers, nsteps, pic=True)

print(T0_fit, T0_sp, T0_sm)

plt.show()
