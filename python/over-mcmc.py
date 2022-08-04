# -*- coding: utf-8 -*-
"""
MCMC realize
"""

import numpy as np
import matplotlib.pyplot as plt
import emcee
from multiprocessing import Pool
from chainconsumer import ChainConsumer

from model import model, spec_model, gauss_model
from trans_function import filtrs


#____________________________DATA_______________________________________

# Paths to files
path_to_SZ_data = "./SZ_data.txt"
path_to_prior = "./prior.dat"
path_to_startSZ = "./startSZ.sss"


# DATA: import data from SZ_data.txt
SZ_data = []
with open(path_to_SZ_data) as file:
    for k, line in enumerate(file):
        if k == 1: z = float(line) # redshift
        if k < 3: continue # skip first 3 lines
        SZ_data.append(line.split())
    SZ_data = np.array(SZ_data, dtype=float)

x = [70.0, 100.0, 143.0, 217.0, 353.0]
if x != filtrs:
    print("Error! fequenses form data and spec_trans don't match")

y = SZ_data[:, 0]
yerr = SZ_data[:, [1, 2]]


# DATA: impot priors from prior.dat
priors_data = dict()
with open(path_to_prior) as file:
    for k, line in enumerate(file):
        if k < 1: continue # skip first line
        prior_param = line.split()[2]
        priors_data[prior_param] = np.array(line.split()[3:6], dtype=float) 


# DATA: import start position from startSZ.sss
with open(path_to_startSZ) as file:
    text = file.readlines()
    nsteps = int(text[3].split()[0])
    nwalkers = int(text[4].split()[0])
    init = []
    for k, line in enumerate(text):
        if k < 6: continue # skip first 3 line
        init.append(line.split()[1:5])
    init = np.array(init, dtype=float)


#____________________________MCMC_______________________________________

# MCMC: model
params_names = [r'T0', r'Te', r'beta', r'Tau']
params_names_latex = [r'$T_0$', r'$T_e$', r'$\beta$', r'$\mathfrac{T}']
ndim = len(params_names)


def log_prior_box(v, vleft, vright):
    if vleft < v < vright:
        return 0.0
    return -np.inf


def log_prior_gauss(v, mean, sigma_p, sigma_m):
    sigma = sigma_p if (v - mean) else sigma_m
    return - 0.5 * (v - mean) ** 2 / sigma ** 2


def prior(params):
    T0, Te, beta, Tau = params
    prior_value = 0
    # box priors
    for k, p in enumerate(params):
        prior_value += log_prior_box(p, init[k, 1], init[k, 2])
    # gauss priors
    for key in priors_data.keys():
        if key == "Te":
            v = priors_data[key]
            prior_value += log_prior_gauss(Te, v[0], v[1], v[2])
    
    return prior_value 


def log_probability(params, x, y, yerr, prior=lambda v: 0):
    if not np.isfinite(prior(params)):
        return -np.inf
    #m = model(*params, x)
    #m = spec_model(*params)
    m = gauss_model(*params)
    N = len(y)
    if np.shape(yerr) == (N,):
        sigma2 = yerr ** 2
    if np.shape(yerr) == (N, 2):
        sigma2 = np.zeros(N)
        for i in range(N):
            sigma2[i] = (yerr[i,1] if m[i] > y[i] else yerr[i, 0]) ** 2
    lp_value = -0.5 * np.sum([(y[i] - m[i]) ** 2 / sigma2[i] for i in range(N)])
    lp_value += prior(params)
    return lp_value


# MCMC: settings
#nwalkers = 100
#nsteps = 200
amputete = int(0.5 * nsteps)
pos = init[:, 0] +  init[:, 3] * np.random.randn(nwalkers, ndim)

# check init posintion in prior-box
for k, param in enumerate(pos):
    p = param
    while not np.isfinite(prior(p)):
        p = [2.72, 6.0, 0.0, 1.5] + np.array([0.5, 2.0, 0.011, 0.5]) * np.random.rand(ndim)
    pos[k] = p

# MCMC: create chain
with Pool() as pool:
    # args(x, y, yerr, prior) - on priors
    # args(x, y, yerr) - off priors
    sampler = emcee.EnsembleSampler(
            nwalkers, ndim, log_probability, args=(x, y, yerr, prior), pool=pool
            )
    sampler.run_mcmc(pos, nsteps, progress=True)

# MCMC: analysis
# amputetion first steps of chain
flat_sample = sampler.chain[:, amputete:, :].reshape((-1, ndim))
c = ChainConsumer()
c.add_chain(flat_sample, parameters=params_names)
summary = c.analysis.get_summary(parameters=params_names)
# summary print
print(*[" {}: {}".format(k, summary[k]) for k in summary.keys()], sep='\n')


#____________________________PIC________________________________________

# PIC: chain-walkers by spteps 
fig, ax = plt.subplots(nrows=ndim, figsize=(10, 7), sharex=True)
samples = sampler.get_chain()
for i, row in enumerate(ax, start=0):
    row.plot(samples[:, :, i], "k", alpha=0.3)
    row.set_ylabel(params_names[i], fontsize=12)
row.set_xlabel(r'steps', fontsize=12)

# PIC: ChainConsumer parameters plot
fig = c.plotter.plot(display=False, legend=False, figsize=(6, 6))

# PIC: data and fit
fig, ax = plt.subplots(figsize=(8, 8))
sample_last = samples[-1, :, :]
for w in sample_last:
    ax.plot(x, model(*w, x), 'b', alpha=0.09)
ax.errorbar(x, y, yerr.T, label='data', capsize=3.5, mew=1.5, fmt='.k', alpha=0.5)
params_fit = [summary[p][1] for p in params_names]
ax.plot(x, model(*params_fit, x), 'g')
params_chi = max(sample_last, key=lambda s: log_probability(s, x, y, yerr, prior))
ax.plot(x, model(*params_chi, x), 'r', label='best fit')

ax.set_xlabel(r"$\nu$, GHz")
ax.set_ylabel(r"SZ signal $\mu$K")
ax.legend(frameon=False)

plt.show()
