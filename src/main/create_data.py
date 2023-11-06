#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Create set input-files for N clusters
"""

# Imports
import numpy as np

import paths # manage path imports by paths.py
from core_sz.data import write_SZdata, write_prior
from core_sz.trans.trans_data import filtrs
from config import config


# random seeds
seed_params = 123
seed_data   = 5769
seed_prior  = 9675


# FIXME: what if left = right
def norm(mu, sigma, left, right, N=1):
    data = np.random.normal(mu, sigma, N)

    for k, r in enumerate(data):
        while not left <= r <= right:
            r = np.random.normal(mu, sigma)
        data[k] = r

    return data


# Paths
path_SZdata = paths.src + '/main/input/szdata/szdata{}.txt'
path_priors = paths.src + '/main/input/priors/prior{}.dat'


# Constants
N = config['N']
T0_param = 2.7255 # K
Te_min, Te_max = 1, 10 # keV
beta_min, beta_max = -0.001, 0.001 # 1/300 +-= 0.003
Tau_min, Tau_max = 0.5, 2
z_mu, z_sigma = 0.0, 0.3

# For Wawes: 70, 100, 143, 217, 353 GHz
rel_error   = np.array([0.4637668, 0.1400576, 0.0686723, 0.0610072, 0.2664959])
rel_error_s = np.array([0.0906936, 0.0298204, 0.0195390, 0.0198462, 0.0810123])
scale_coeff, sigma_coeff = 0.9655, 0.4632


# Generate params
np.random.seed(seed_params)
T0 = np.full(N, T0_param)
Te = np.random.uniform(Te_min, Te_max, N)
beta = np.random.uniform(beta_min, beta_max, N)
Tau = np.random.uniform(Tau_min, Tau_max, N)
z = np.abs(np.random.normal(z_mu, z_sigma, N))
#Tz = T0 * (1 + z)

model = config['model']
params_include = config['params_include']
params = np.array(params_include(T0, Te, beta, Tau, z=z)).T
# beta, Tau, = np.zeros(N), np.ones(N) # if you use model_000
#params = np.array([T0, Te, beta, Tau]).T



def gen_data(cluster):
    # data -> scale of errors (+ rand) -> errors (+ rand) -> signal (+ rand)
    sz = model(*cluster)
    scale_sz = max(sz) - min(sz)
    scale_error = scale_sz * norm(scale_coeff, sigma_coeff, left=0.3, right=1.3)
    error = scale_error * np.random.normal(rel_error, rel_error_s)
    signal = np.random.normal(sz, 0.8 * np.abs(error))
    data = np.array([signal, error, error]).T #symmetric errors
    #data = np.array([signal, errors + (sz - sigma), errors - (sz - sigma]).T
    return data

def gen_prior(Te: float):
    Te_mean = np.random.normal(Te, 0.001* Te)
    Te_sigma = 2 * [np.random.normal(0.05 * Te, 0.0025 * Te)] # symmetric errors
    #Te_sigma = np.random.normal(0.05 * Te, 0.0025 * Te, 2)
    prior_gauss_data = {'Te': [Te_mean, *Te_sigma]}
    return prior_gauss_data


# Generate data
np.random.seed(seed_data)
sz_data = tuple(map(gen_data, params))

# Generate prior
np.random.seed(seed_prior)
prior_gauss_data = tuple(map(gen_prior, Te))


if __name__ == "__main__":
    # Write files: SZdata.txt
    for i in range(N):
        path = path_SZdata.format(i + 1)
        comment = f"for freq. {str(filtrs)} GHz | "
        comment += "T0 = {}, kTe = {:.2}, beta = {:.1e}, Tau = {:.2}"
        comment = comment.format(T0[i], Te[i], beta[i], Tau[i])
        write_SZdata(path, z[i], sz_data[i], comment)

    # Write files: prior.dat
    for i in range(N):
        path = path_priors.format(i + 1)
        write_prior(path, prior_gauss_data[i])

