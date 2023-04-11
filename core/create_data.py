# -*- coding: utf-8 -*-
"""
Create set input-files for N clusters
comment: Давайте в SZ data записывать и частоты
"""

import numpy as np

from data import SZ_data_write, prior_write
from model import *
from trans_function import filtrs


np.random.seed(123)

def norm(N, mu, sigma, left, right):
    data = np.random.normal(mu, sigma, N)

    for k, r in enumerate(data):
        while not left <= r <= right:
            r = np.random.normal(mu, sigma)
        data[k] = r

    return data

# Paths
path_to_SZ_data = '../data/set/szdata/szdata{}.txt'
path_to_priors = '../data/set/priors/prior{}.dat'

# constants
N = 1000
T0 = 2.7255 # K
Te_min = 1 # keV
Te_max = 10 # keV
beta_min = - 0.5 * 1/300
beta_max = 0.5 * 1/300
Tau_min = 0.5
Tau_max = 2
yerr_rel = np.array([0.4637668, 0.1400576, 0.0686723, 0.0610072, 0.2664959])
yerr_rel_s = np.array([0.0906936, 0.0298204, 0.0195390, 0.0198462, 0.0810123])
rscale, rscale_s = 0.9655, 0.4632
z_mu = 0.0
z_sigma = 0.3

model = simple_model

# Generate params
Te = np.random.uniform(Te_min, Te_max, N)
beta = np.random.uniform(beta_min, beta_max, N)
Tau = np.random.uniform(Tau_min, Tau_max, N)
z = np.abs(np.random.normal(z_mu, z_sigma, N))

if __name__ == "__main__":
    # Write files: SZdata.txt
    for i in range(N):
        path = path_to_SZ_data.format(i + 1)
        comment = 'for freq. {} GHz | '.format(str(filtrs))
        comment += 'T0 = {}, kTe = {:.2}, beta = {:.1e}, Tau = {:.2}'
        comment = comment.format(T0, Te[i], beta[i], Tau[i])
        #sz = model(T0, Te[i], beta[i], Tau[i])
        sz = model(T0, Te[i])

        # symetric errors
        #rscale_g = np.random.normal(rscale, rscale_s)
        rscale_g = norm(1, rscale, rscale_s, 0.3, 1.3)
        scale = max(sz) - min(sz)
        yerr_rel_g = np.random.normal(yerr_rel, yerr_rel_s)
        sz_sigma1 = np.abs(rscale_g * scale * yerr_rel_g)
        sz_sigma = np.array([sz_sigma1, sz_sigma1])
        sz_g = np.random.normal(sz, 0.8 * sz_sigma1)

        sz_data = np.array([sz_g, *sz_sigma]).T

        SZ_data_write(path, z[i], sz_data, comment)

    # Write files: prior.dat
    for i in range(N):
        path = path_to_priors.format(i + 1)
        Te_mean = np.random.normal(Te[i], 0.01 * Te[i])
        Te_sigma = np.random.normal(0.05 * Te[i], 0.05 * 0.05 * Te[i], 2)
        prior_gauss_data = {'Te': [Te_mean, *Te_sigma]}

        prior_write(path, prior_gauss_data)

