# -*- coding: utf-8 -*-
"""
Analyze example cluster
comment: добавить сохранение картинок
"""
#import sys

## setting path
#sys.path.append('./core')


import numpy as np
import matplotlib.pyplot as plt

from mcmc_analyzeTz import mcmc_analyze


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


# paths
path_SZ_data = '../data/N/szdata/szdata{}.txt'
path_startSZ = '../data/startSZ.sss'
path_prior = '../data/N/priors/prior{}.dat'

path_result = '../data/N/resultTz.txt'

N = 20

if __name__ == "__main__":
    for i in range(N):
        paths = [path_SZ_data.format(i + 1),
                 path_startSZ,
                 path_prior.format(i + 1)]

        nwalkers, nsteps = 100, 100
        z, *T0_params= mcmc_analyze(paths, nwalkers, nsteps, pic=False)
        result_write(z, T0_params, path_result)
        print(i + 1, z, T0_params)

plt.show()
