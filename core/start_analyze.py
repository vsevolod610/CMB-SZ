# -*- coding: utf-8 -*-
"""
Analyze example cluster
"""

import numpy as np
import matplotlib.pyplot as plt

from mcmc_analyze import SZ_mcmc


def result_read(path):
    res = []
    with open(path) as file:
        for k, line in enumerate(file, start=1):
            if 'None' in line.split():
                continue
            res.append(line.split())
        res = np.array(res, dtype=float)
    z = res[:, 0]
    T0 = res[:, 1]
    T0_err = res[:, [2, 3]]

    return z, T0, T0_err


def result_write(z, T0_params, path):
    if not None in T0_params:
        s = "{:<8}   {:<10.6} {:<10.6} {:<10.6}".format(z, *T0_params)
    else:
        s = "{:<8}   {:>30}".format(z, 'None')

    with open(path, 'a') as file:
        file.write(s + "\n")


method = 'T0'
#method = 'Tz'

# paths
path_SZ_data = '../data/N/szdata/szdata{}.txt'
path_startSZ = '../data/startSZ.sss'
path_prior = '../data/N/priors/prior{}.dat'

path_result = '../data/N/result-{}.txt'.format(method)
path_pics = ['../data/N/pic_chain/pic_chain-{}',
             '../data/N/pic_consumer/pic_consumer-{}',
             '../data/N/pic_fit/pic_fit-{}']

N = 200

if __name__ == "__main__":
    for i in range(N):
        paths = [path_SZ_data.format(i + 1),
                 path_startSZ,
                 path_prior.format(i + 1)]

        #nwalkers, nsteps = 100, 100
        nwalkers, nsteps = 'Read', 'Read'

        pics = [s.format(i + 1) for s in path_pics]
        z, *T0_params = SZ_mcmc(method, paths, nwalkers, nsteps, pics=pics)
        result_write(z, T0_params, path_result)
        print(i + 1, z, T0_params)

plt.show()