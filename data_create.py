# -*- coding: utf-8 -*-
""""
Create set input-files for N clusters
"""

import numpy as np

from model import gauss_model
from SZfunction import SZfunction
from trans_function import spec_trans, filtrs


np.random.seed(123)

# Paths
path_to_NSZ_data = './data/N/SZdata/SZdata{}.txt'
path_to_Npriors = './data/N/priors/prior{}.dat'

# constants
N = 2000
T0 = 2.7255 # K
Te_min = 1 # keV
Te_max = 10 # keV
beta_min = - 0.5 * 1/300
beta_max = 0.5 * 1/300
Tau_min = 0.1
Tau_max = 2
sz_relerr = dict(zip(filtrs,[0.22, 0.08, 0.1, 0.98, 0.30]))
z_sigma = 0.3

# Generate params
Te = np.random.uniform(Te_min, Te_max, N)
beta = np.random.uniform(beta_min, beta_max, N)
Tau = np.random.uniform(Tau_min, Tau_max, N)
z = np.abs(np.random.normal(0.0, z_sigma, N))


for i in range(N):
    # SZdata.txt
    with open(path_to_NSZ_data.format(i + 1), 'w') as File:
        File.write('SZ data for freq. {} GHz | '.format(str(filtrs[1:-1])))
        line_params = 'T0 = {}, kTe = {:.2}, beta = {:.1e}, Tau = {:.2}'
        File.write(line_params.format(T0, Te[i], beta[i], Tau[i]))
        File.write('\n{:.4}'.format(z[i]))
        File.write('\n{}\n'.format(len(filtrs)))
        sz = dict(zip(filtrs, gauss_model(T0, Te[i], beta[i], Tau[i])))
        for wave in filtrs:
            sz_sigma = np.abs(sz[wave] * np.random.normal(sz_relerr[wave], sz_relerr[wave] / 5, 2))
            s = '{: .2e}  {:.2e}  {:.2e}\n'.format(sz[wave], sz_sigma[0], sz_sigma[1])
            File.write(s)
    
    # prior.dat
    with open(path_to_Npriors.format(i + 1), 'w') as File:
        File.write('1 ! Num of priors\n')
        s = '1 2 Te {:.2} {:.2} {:.2} ! 3 index of arg Te name of arg'
        Te_mean = np.random.normal(Te[i], 0.01 * Te[i])
        Te_sigma = np.random.normal(0.05 * Te[i], 0.05 * 0.05 * Te[i], 2)
        s = s.format(Te_mean, Te_sigma[0], Te_sigma[1])
        File.write(s)


