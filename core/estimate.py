# -*- coding: utf-8 -*-
"""
set of T0_i --> Estimate T0
comment: нодо бы перенести этот файл в корень
"""

import numpy as np
import matplotlib.pyplot as plt
from chainconsumer import ChainConsumer

from mcmc import mcmc
from start_analyze import result_read, path_result, method


x, y, yerr = result_read(path_result)


def modelT0(T0, z):
    return T0 + (0 * z)


def modelTz(T0, z):
    return T0 * (1 + z)


if method == 'T0':
    model = modelT0
if method == 'lazy':
    model = modelT0
if method == 'Tz':
    model = modelTz

init = np.array([[2.7255], [0.2]]).T
nwalkers = 200
nsteps = 200
amputete = int(0.2 * nsteps)
ndim = len(init)

path_pics = ['../res/pic1.pdf', '../res/pic2.pdf', '../res/pic3.pdf']

# mcmc for estimation
summary = mcmc(data=(x, y, yerr),
               model_params=(model, init, None, None),
               settings=(nwalkers, nsteps, amputete),
               prnt=False, show=False, save=path_pics)

T0 = summary['0'][1]
T0p, T0m = summary['0'][2] - T0, T0 - summary['0'][0]
print('T0 = ', T0, ' +', T0p, ' -', T0m)
