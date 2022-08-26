
"""
MCMC realize
"""

import numpy as np
import matplotlib.pyplot as plt
import emcee
from multiprocessing import Pool
from chainconsumer import ChainConsumer

from model import SZmodel, spec_model, gauss_model
from data import x, y, yerr, init, prior_data, ndim, nsteps, nwalkers, params_names


def log_prior_box(v, vleft, vright):
    if vleft < v < vright:
        return 0.0
    return -np.inf


def log_prior_gauss(v, mean, sigma_p, sigma_m):
    sigma = sigma_p if (v - mean) else sigma_m
    return - 0.5 * (v - mean) ** 2 / sigma ** 2


def prior(params, prior_data):
    T0, Te, beta, Tau = params
    prior_value = 0

    # box priors
    for k, p in enumerate(params):
        left, right = prior_data['box'][k]
        prior_value += log_prior_box(p, left, right)

    # gauss priors
    for key in prior_data['gauss'].keys():
        if key == "Te":
            mu, sigmap, sigmam = prior_data['gauss'][key]
            prior_value += log_prior_gauss(Te, mu, sigmam, sigmap)

    return prior_value


def log_probability(params, x, y, yerr, prior_data):
    prior_value = prior(params, prior_data)
    if not np.isfinite(prior_value):
        return -np.inf
    m = gauss_model(*params)
    N = len(y)
    sigma2 = np.zeros(N)
    for i in range(N):
        sigma2[i] = (yerr[i,1] if m[i] > y[i] else yerr[i, 0]) ** 2
    lp_value = -0.5 * np.sum([(y[i] - m[i]) ** 2 / sigma2[i] for i in range(N)])
    lp_value += prior_value
    return lp_value


def SZmcmc(x, y, yerr, nwalkers, nsteps, ndim, init, prior_data):
    # MCMC: settings
    pos = init[:, 0] + init[:, 1] * np.random.randn(nwalkers, ndim)

    # check init posintion in prior-box
    for k, param in enumerate(pos):
        while not np.isfinite(prior(param, prior_data)):
            param = init[:, 0] + init[:, 1] * np.random.rand(ndim)
        pos[k] = param

    # MCMC: create chain
    with Pool() as pool:
        sampler = emcee.EnsembleSampler(nwalkers, ndim, 
                log_probability, args=(x, y, yerr, prior_data), pool=pool)
        sampler.run_mcmc(pos, nsteps, progress=True)

    return sampler


def pic_chain(sampler):
    fig, ax = plt.subplots(nrows=ndim, figsize=(10, 7), sharex=True)
    samples = sampler.get_chain()

    for i, row in enumerate(ax, start=0):
        row.plot(samples[:, :, i], "k", alpha=0.3)
        row.set_ylabel(params_names[i], fontsize=12)

    row.set_xlabel(r'steps', fontsize=12)
    return fig


def pic_fit(sampler, x, y, yerr, prior_data):
    fig, ax = plt.subplots(figsize=(8, 8))
    samples = sampler.get_chain()
    sample_last = samples[-1, :, :]
    params_chi = max(sample_last, 
            key=lambda s: log_probability(s, x, y, yerr, prior_data))

    for w in sample_last:
        ax.plot(x, SZmodel(*w, x), 'b', alpha=0.09)
    ax.errorbar(x, y, yerr.T, label='data', 
            capsize=3.5, mew=1.5, fmt='.k', alpha=0.5)
    ax.plot(x, 0 * x, 'k')
    ax.plot(x, SZmodel(*params_chi, x), 'r', label='best fit')

    ax.set_xlabel(r"$\nu$, GHz")
    ax.set_ylabel(r"SZ signal $\mu$K")
    ax.legend(frameon=False)
    return fig


if __name__ == "__main__":
    # MCMC: analysis
    sampler = SZmcmc(x, y, yerr, nwalkers, nsteps, ndim, init, prior_data)
    amputete = int(0.5 * nsteps)
    flat_sample = sampler.chain[:, amputete:, :].reshape((-1, ndim))
    c = ChainConsumer()
    c.add_chain(flat_sample, parameters=params_names)
    summary = c.analysis.get_summary(parameters=params_names)

    print("\nMCMC results:")
    print(*[" {:>4}: {}".format(k, summary[k]) for k in summary.keys()], sep='\n')

    # Pics
    fig = pic_chain(sampler)
    fig = c.plotter.plot(display=False, legend=False, figsize=(6, 6))
    fig = pic_fit(sampler, x, y, yerr, prior_data)

    plt.show()
