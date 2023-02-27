# -*- coding: utf-8 -*-
"""
Create set input-files for N clusters
comment: Давайте в SZ data записывать и частоты
"""

import numpy as np

from data import SZ_data_write, prior_write
from model import gauss_model
from trans_function import filtrs


np.random.seed(123)

# Paths
path_to_NSZ_data = '../data/N/szdata/szdata{}.txt'
path_to_Npriors = '../data/N/priors/prior{}.dat'

# constants
N = 1000
T0 = 2.7255 # K
Te_min = 1 # keV
Te_max = 10 # keV
beta_min = - 0.5 * 1/300
beta_max = 0.5 * 1/300
Tau_min = 0.5
Tau_max = 2
sz_rerr = np.array([0.22, 0.08, 0.1, 0.98, 0.30])
z_mu = 0.0
z_sigma = 0.3

# Generate params
Te = np.random.uniform(Te_min, Te_max, N)
beta = np.random.uniform(beta_min, beta_max, N)
Tau = np.random.uniform(Tau_min, Tau_max, N)
z = np.abs(np.random.normal(z_mu, z_sigma, N))

# Write files: SZdata.txt
for i in range(N):
    path = path_to_NSZ_data.format(i + 1)
    comment = 'for freq. {} GHz | '.format(str(filtrs))
    comment += 'T0 = {}, kTe = {:.2}, beta = {:.1e}, Tau = {:.2}'
    comment = comment.format(T0, Te[i], beta[i], Tau[i])
    sz = gauss_model(T0, Te[i], beta[i], Tau[i])

    # asymetric errors
    sz_sigma = np.abs(
            sz * np.random.normal(sz_rerr, sz_rerr / 5, (2, len(filtrs))))

    # symetric errors
    sz_sigma1 = np.abs(sz * np.random.normal(sz_rerr, sz_rerr / 5, len(filtrs)))
    sz_sigma = np.array([sz_sigma1, sz_sigma1])

    sz_data = np.array([sz, *sz_sigma]).T

    SZ_data_write(path, z[i], sz_data, comment)

# Write files: prior.dat
for i in range(N):
    path = path_to_Npriors.format(i + 1)
    Te_mean = np.random.normal(Te[i], 0.01 * Te[i])
    Te_sigma = np.random.normal(0.05 * Te[i], 0.05 * 0.05 * Te[i], 2)
    prior_gauss_data = {'Te': [Te_mean, *Te_sigma]}

    prior_write(path, prior_gauss_data)

