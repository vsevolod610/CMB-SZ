# -*- coding: utf-8 -*-
""""
set of T0_i temperature --> Estimate T0
"""

import numpy as np
import matplotlib.pyplot as plt
import emcee
from multiprocessing import Pool
from chainconsumer import ChainConsumer



# data
res = []
with open('./Result/Result.txt') as f:
    for k, line in enumerate(f, start=1):
        if 'None' in line.split():
            continue
        res.append(line.split())
    res = np.array(res, dtype=float)
x = res[:, 0]
y = res[:, 2]
yerr = res[:, [3, 4]]


#____________________________________mcmc_________________________________

def log_prior_box(v, vleft, vright):
    if vleft < v < vright:
        return 0.0
    return -np.inf

def log_prior_gauss(v, mean, sigma_p, sigma_m):
    sigma = sigma_p if (v - mean) else sigma_m
    return - 0.5 * (v - mean) ** 2 / sigma ** 2


def log_probability(params, x, y, yerr, prior=lambda v: 0):
    if not np.isfinite(prior(params)):
        return -np.inf
    m = model(*params, x)
    sigma2 =  np.array([(yerr[i, 1] if m[i] > yi else yerr[i, 0]) for i, yi in enumerate(y)]) ** 2
    return - 0.5 * np.sum([(yi - m[i]) ** 2 / sigma2[i] for i, yi in enumerate(y)]) + prior(params)


def model(T0, z):
    return T0 + 0 * z

# setting of mcmc
nwalkers = 200
nsteps = 500
amputete = 100
ndim = 1
names_of_params = [r'T0']
pos = [2.72] +  0.5 * np.random.randn(nwalkers, ndim)

# mcmc mechanism
with Pool() as pool:
    sampler = emcee.EnsembleSampler( nwalkers, ndim, log_probability, args=(x, y, yerr), pool=pool)
    sampler.run_mcmc(pos, nsteps, progress=True)
flat_sample = sampler.chain[:, amputete : , :].reshape((-1, ndim))
c = ChainConsumer()
c.add_chain(flat_sample, parameters=names_of_params)
summary = c.analysis.get_summary(parameters=names_of_params)[names_of_params[0]]
print(summary)

#____________________________________pic__________________________________

# chainconsum pic
fig = c.plotter.plot(display=False, legend=False, figsize=(6, 6))
fig.savefig('T_0 dist.pdf', bbox_inches='tight')

# MCchain pic
fig, ax = plt.subplots(figsize=(10, 7), sharex=True)
samples = sampler.get_chain()
ax.plot(samples[:, :, 0], "k", alpha=0.3)
ax.set_xlabel(r'steps', fontsize=12)
ax.set_ylabel(r'$T_0$', fontsize=12)
fig.savefig('MCchain.pdf', bbox_inches='tight')

# data pic
fig, ax = plt.subplots(figsize=(8, 8))
ax.errorbar(x, y, np.transpose(yerr), capsize=3.5, mew=1.5, fmt='.k', alpha=0.5)             
xline = np.linspace(min(x), max(x), 100)
ax.plot(xline, 0 * xline, 'k')
ax.plot(xline, model(summary[1], xline))
ax.fill_between(xline, model(summary[0], xline), model(summary[2], xline),
            facecolor='b', alpha=0.3, color='b', linewidth=2, linestyle='-')
ax.set_xlabel(r'z', fontsize=12)
ax.set_ylabel(r'$T_z$', fontsize=12)
fig.savefig('Data.pdf', bbox_inches='tight')
