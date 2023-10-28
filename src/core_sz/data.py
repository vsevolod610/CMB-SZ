#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Read Data files: SZ_data.txt, prior.dat, startSZ.sss
"""

import numpy as np


#-----------DATA: SZ_data.txt---------------

def read_SZdata(path):
    with open(path) as file:
        lines = file.readlines()
        z = float(lines[1].split()[0])

        data = []
        for line in lines[3:]:
            data.append(line.split()[:3])
        data = np.array(data, dtype=float)

    x = np.array([70.0, 100.0, 143.0, 217.0, 353.0])
    y = data[:, 0]
    yerr = data[:, [1, 2]]

    return x, y, yerr, z

def write_SZdata(path, z, sz_data, comment=""):
    with open(path, 'w') as file:
        file.write("SZ data" + comment)
        file.write(f"\n{z:<8.4} | z")
        file.write(f"\n{len(sz_data)} | num\n")
        for dot in sz_data:
            s = '{: .2e}  {:.2e}  {:.2e}\n'.format(*dot)
            file.write(s)

#-----------DATA: prior.dat-----------------

def read_prior(path):
    with open(path) as file:
        lines = file.readlines()

        prior_gauss = dict()
        for line in lines[1:]:
            param = int(line.split()[1]) - 1 # if 2 then Te
            prior_gauss[param] = np.array(line.split()[3:6], dtype=float)

    return prior_gauss

def write_prior(path, prior_gauss_data):
    with open(path, 'w') as file:
        file.write(f"{len(prior_gauss_data)} ! Num of priors\n")
        for param in list(prior_gauss_data):
            s = '1 2 {} {:.6} {:.6} {:.6} ! 3 index of arg Te name of arg'
            s = s.format(param, *prior_gauss_data[param])
            file.write(s)

#-----------DATA: startSZ.sss-----------------

def read_startSZ(path):
    with open(path) as file:
        lines = file.readlines()
        nsteps = int(lines[3].split()[0])
        nwalkers = int(lines[4].split()[0])

        data = []
        for line in lines[6:]:
            data.append(line.split()[1:5])
        data = np.array(data, dtype=float)

        init = data[:, [0, 3]]

        prior_box = data[:, [1,2]]
        ndim = len(data)
        prior_box = dict(zip(np.arange(ndim),prior_box))

    return nwalkers, nsteps, init, prior_box

