# -*- coding: utf-8 -*-
"""
Analyze set of cluster
"""

# Imports
import numpy as np

import paths # manage path imports by paths.py
from core_sz.sz_mcmc import sz_analyze
from config import config


def read_result(path):
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


def write_result(z, T0_params, path):
    if not None in T0_params:
        s = "{:<10}   {:<10.6} {:<10.6} {:<10.6}".format(z, *T0_params)
    else:
        s = "{:<10}   {:>30}".format(z, 'None')

    with open(path, 'a') as file:
        file.write(s + "\n")


# config
N = config['N']
case = config['case']
model = config['model']


# paths
paths_input = [
        paths.src + '/main/input/szdata/szdata{}.txt',
        paths.src + '/main/input/priors/prior{}.dat',
        config['path_startSZ']]
        # feauture not a bag: path_startSZ is not format(i)

paths_pics = [
        paths.src + '/main/output/pic_chain/pic_chain-{}', 
        paths.src + '/main/output/pic_params/pic_params-{}', 
        paths.src + '/main/output/pic_fit/pic_fit-{}' ]

path_result = paths.src + f'/main/output/result-{case}.txt'
#pathCH = '../data/set/OUT/chainconsum{}.pkl'


if __name__ == "__main__":
    for i in range(1, N + 1):
        iterize = lambda s: s.format(i)
        path_input = list(map(iterize, paths_input))
        path_pics = list(map(iterize, paths_pics))

        z, *T0_params = sz_analyze(model, path_input, pics=path_pics)
        write_result(z, T0_params, path_result)
        print(f"\n {i:>{len(str(N))}}: {z=}, {T0_params}\n")

