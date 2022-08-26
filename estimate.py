# -*- coding: utf-8 -*-
"""
set of T0_i --> Estimate T0
"""

import numpy as np
import matplotlib.pyplot as plt
import emcee
from multiprocessing import Pool
from chainconsumer import ChainConsumer

from poly_mcmc import path_result, read_result


x, y, yerr = read_result(path_result)


def model(T0, z):
    return T0 + 0 * z


def log_probability(params, x, y, yerr):
    m = model(*params, x)
    N = len(y)
    sigma2 = np.zeros(N)
    for i in range(N):
        sigma2[i] = (yerr[i,1] if m[i] > y[i] else yerr[i, 0]) ** 2
    lp_value = -0.5 * np.sum([(y[i] - m[i]) ** 2 / sigma2[i] for i in range(N)])
    return lp_value


# setting of mcmc
nwalkers = 100
nsteps = 200
amputete = 100

ndim = 1
names_of_params = [r'T0']
pos = [2.72] +  0.5 * np.random.randn(nwalkers, ndim)

# mcmc mechanism
with Pool() as pool:
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, args=(x, y, yerr), pool=pool)
    sampler.run_mcmc(pos, nsteps, progress=True)

flat_sample = sampler.chain[:, amputete : , :].reshape((-1, ndim))
c = ChainConsumer()
c.add_chain(flat_sample, parameters=names_of_params)
summary = c.analysis.get_summary(parameters=names_of_params)[names_of_params[0]]
print(summary)


# pic's
fig = c.plotter.plot(display=False, legend=False, figsize=(6, 6))

plt.show()
