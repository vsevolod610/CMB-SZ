# -*- coding: utf-8 -*-
"""
Read Data
comment: формат файлов SZ_data и prior устарел
         не стоит здесь определять params_names
         стоит определиться что сначала nsteps или nwalkers
"""

import numpy as np


def SZ_data_read(path):
    # DATA: import data from SZ_data.txt
    data = []
    with open(path) as file:
        for k, line in enumerate(file):
            if k == 1: z = float(line.split()[0])
            if k < 3: continue
            data.append(line.split()[:3])
        data = np.array(data, dtype=float)

    x = np.array([70.0, 100.0, 143.0, 217.0, 353.0])
    y = data[:, 0]
    yerr = data[:, [1, 2]]

    return x, y, yerr, z

def SZ_data_write(path, z, sz_data, comment=""):
    with open(path, 'w') as file:
        file.write('SZ data' + comment)
        file.write('\n{:<8.4} | z'.format(z))
        file.write('\n{} | num\n'.format(len(sz_data)))
        for dot in sz_data:
            s = '{: .2e}  {:.2e}  {:.2e}\n'.format(*dot)
            file.write(s)


def prior_read(path):
    # DATA: impot prior from prior.dat
    prior_gauss = dict()
    with open(path) as file:
        for k, line in enumerate(file):
            if k < 1: continue 
            param = int(line.split()[1]) - 1 # if 2 then Te
            prior_gauss[param] = np.array(line.split()[3:6], dtype=float)

    return prior_gauss

def prior_write(path, prior_gauss_data):
    with open(path, 'w') as file:
        file.write('{} ! Num of priors\n'.format(len(prior_gauss_data)))
        for param in list(prior_gauss_data):
            s = '1 2 {} {:.2} {:.2} {:.2} ! 3 index of arg Te name of arg'
            s = s.format(param, *prior_gauss_data[param])
            file.write(s)


def startSZ_read(path):
    # DATA: import start position from startSZ.sss
    with open(path) as file:
        text = file.readlines()
        nsteps = int(text[3].split()[0])
        nwalkers = int(text[4].split()[0])
        data = []
        for k, line in enumerate(text):
            if k < 6: continue # skip first 3 line
            data.append(line.split()[1:5])
        data = np.array(data, dtype=float)

        init = data[:, [0, 3]]

        prior_box = data[:, [1,2]]
        ndim = len(data)
        prior_box = dict(zip(np.arange(ndim),prior_box))

    return nwalkers, nsteps, init, prior_box

#   # Example
#   # Paths to files
#   path_to_SZ_data = "./data/SZ_data.txt"
#   path_to_prior = "./data/prior.dat"
#   path_to_startSZ = "./data/startSZ.sss"

#   # Names
#   params_names = [r'T0', r'Te', r'beta', r'Tau']
#   params_names_latex = [r'$T_0$', r'$T_e$', r'$\beta$', r'$\mathfrac{T}']
#   ndim = len(params_names)

#   # Set data
#   x, y, yerr, z = SZ_data_read(path_to_SZ_data)
#   nwalkers, nsteps, init, prior_box = startSZ_read(path_to_startSZ)
#   prior_data = dict()
#   prior_data['box'] = prior_box
#   prior_data['gauss'] = prior_read(path_to_prior)
#   #ndim = len(init)

