# -*- coding: utf-8 -*-
""""
set of T0_i temperature --> Estimate T0

require:    ../Result/Result.txt
make:       none
            

coments:    
"""

import numpy as np
import matplotlib.pyplot as plt
import emcee
from chainconsumer import ChainConsumer


# data
res = []
with open('../Result/Result.txt') as f:
    for k, line in enumerate(f, start=1):
        if 'None' in line.split():
            continue
        res.append(line.split())
    res = np.array(res, dtype=float)
n = res[:, 0]
y = res[:, 1]
yerr = res[:, [2, 3]]

#____________________________________mcmc_________________________________

def log_prior_box(v, vleft, vright):
    if vleft < v < vright:
        return 0.0
    return -np.inf

def log_prior_gauss(v, mean, sigma_p, sigma_m):
    sigma = sigma_p if (v - mean) else sigma_m
    return - 0.5 * (v - mean) ** 2 / sigma ** 2


def log_probability(v, y, yerr, prior=lambda v: 0):
    if not np.isfinite(prior(v)):
        return -np.inf
    sigma2 = np.array([(yerr[i, 0] if v > yi else yerr[i, 1]) for i, yi in enumerate(y)]) ** 2
    s2 = np.array([(yerr[i, 0] + yerr[i, 1]) / 2 for i, yi in enumerate(y)]) ** 2
    return - 0.5 * np.sum([(yi - v) ** 2 / sigma2[i] + 2 * np.pi * s2[i]for i, yi in enumerate(y)]) + prior(v)

# setting of mcmc
nwalkers = 100
nsteps = 200
amputete = 50 # int(0.5 * nsteps)
ndim = 1
names_of_params = [r'T0']
pos = [2.72] +  0.5 * np.random.randn(nwalkers, ndim)
prior = lambda v: log_prior_box(v, 0, 20)
log_prob = lambda v: log_probability(v, y, yerr, prior)
# mcmc mechanism
sampler = emcee.EnsembleSampler(nwalkers, ndim, log_prob)
sampler.run_mcmc(pos, nsteps, progress=True)
flat_samples = sampler.chain[:, amputete : , :].reshape((-1, ndim))
c = ChainConsumer()
c.add_chain(flat_samples, parameters=names_of_params)
summary = c.analysis.get_summary(parameters=names_of_params)[names_of_params[0]]
print(summary)

#____________________________________pic__________________________________

# chainconsum pic
fig = c.plotter.plot(display=False, legend=False, figsize=(6, 6))
# MCchain pic
fig, ax = plt.subplots(figsize=(10, 7), sharex=True)
samples = sampler.get_chain()
ax.plot(samples[:, :, 0], "k", alpha=0.3)
ax.set_xlabel(r'steps', fontsize=12)
ax.set_ylabel(r'$T_0$', fontsize=12)
# data pic
fig, ax = plt.subplots(figsize=(8, 8))
ax.errorbar(n, y, np.transpose(yerr), capsize=3.5, mew=1.5, fmt='.k')                 
x = np.linspace(0, 77, 100)
ax.plot(x, summary[1] + 0 * x)
ax.fill_between(x, summary[0] + 0 * x, summary[2] + 0 * x,
            facecolor='b', alpha=0.3, color='b', linewidth=2, linestyle='-')
ax.set_xlabel(r'number of cluster', fontsize=12)
ax.set_ylabel(r'$T_0$', fontsize=12)
# chi2 pic
fig, ax = plt.subplots(figsize=(8, 8))
T = np.linspace(2.5, 3., 100)
ax.plot(T, [-2 * log_prob(t) for t in T])
ax.set_xlabel(r'$T_0$', fontsize=12)
ax.set_ylabel(r'$\chi^2$', fontsize=12)
