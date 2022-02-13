# -*- coding: utf-8 -*-
""""
Create set Input files for N clusters

require:    SZfunction.py, trans_function.py
make:       ../INput/SZdatas/SZ_data{}.txt, ../INput/priors/prior{}.dat

coments:    Remember about '1 ! Num of priors'
"""

import numpy as np
from SZfunction import SZfunction
from trans_function import spec_trans, filtrs

np.random.seed(123)

N = 77
T0 = 2.7255 # K
Te_min = 1 # keV
Te_max = 10 # keV
beta_min = - 0.5 * 1/300
beta_max = 0.5 * 1/300
Tau_min = 0.1
Tau_max = 2
sz_relerr = dict(zip(filtrs,[0.22, 0.08, 0.1, 0.98, 0.30]))

Te = np.random.uniform(Te_min, Te_max, N)
beta = np.random.uniform(beta_min, beta_max, N)
Tau = np.random.uniform(Tau_min, Tau_max, N)

def SZtr(T0, kTe, beta, Tau):
    sz = dict.fromkeys(filtrs)
    for wave in filtrs:
        nu = spec_trans[wave]['nu']
        tr = spec_trans[wave]['tr']
        sz[wave] = np.trapz(SZfunction(T0, kTe, beta, Tau, nu) * tr, nu)
    return sz


for i in range(N):
    with open("../INput/SZdatas/SZ_data{}.txt".format(i + 1), 'w') as File:
        File.write('SZ data for freq. {} GHz | '.format(str(filtrs[1:-1])))
        line_params = 'T0 = {}, kTe = {:.2}, beta = {:.1e}, Tau = {:.2}'
        File.write(line_params.format(T0, Te[i], beta[i], Tau[i]))
        File.write('\n{}\n'.format(len(filtrs)))
        sz = SZtr(T0, Te[i], beta[i], Tau[i])
        for wave in filtrs:
            sz_sigma = np.abs(sz[wave] * np.random.normal(sz_relerr[wave], sz_relerr[wave] / 5, 2))
            s = '{: .2e}  {:.2e}  {:.2e}\n'.format(sz[wave], sz_sigma[0], sz_sigma[1])
            File.write(s)
    
    with open("../INput/priors/prior{}.dat".format(i + 1), 'w') as File:
        File.write('1 ! Num of priors\n')
        s = '1 2 Te {:.2} {:.2} {:.2} ! 2 index of arg Te name of arg'
        Te_mean = np.random.normal(Te[i], 0.01 * Te[i])
        Te_sigma = np.random.normal(0.05 * Te[i], 0.05 * 0.05 * Te[i], 2)
        s = s.format(Te_mean, Te_sigma[0], Te_sigma[1])
        File.write(s)


