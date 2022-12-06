# -*- coding: utf-8 -*-
"""
Analyze clusters
"""

import gc
import numpy as np
import matplotlib.pyplot as plt
from chainconsumer import ChainConsumer

from model import gauss_model, alt_model
from mcmc_kern import mcmc_kern, pic_chain, pic_fit
from data import SZ_data_read, startSZ_read, prior_read


def mcmc_analyze(model, mcmc_data, nwalkers=100, nsteps=100, amputate=0.5,
                 prnt=False, pics=False):
    #data
    init, x, y, yerr, prior_data = mcmc_data
    amputete_step = int(amputate * nsteps)
    ndim = len(init)

    # mcmc
    sampler = mcmc_kern(model, nwalkers, nsteps, 
                        init, x, y, yerr, prior_data)
    flat_sample = sampler.chain[:, amputete_step : , :].reshape((-1, ndim))
    c = ChainConsumer()
    c.add_chain(flat_sample)
    summary = c.analysis.get_summary()
    if prnt:
        print(summary)

    # pics
    if pics:
        fig0, ax = pic_chain(sampler, params_names=None)
        fig1 = c.plotter.plot(display=False, legend=False, figsize=(6, 6))
        fig2, ax = pic_fit(sampler, model, x, y, yerr, prior_data)
        if type(pics) is list:
            fig0.savefig(pics[0])
            fig1.savefig(pics[1])
            fig2.savefig(pics[2])

            # garved collector
            plt.clf()
            plt.close()
            plt.close(fig0)
            plt.close(fig1)
            plt.close(fig2)
            plt.close('all')
            del sampler
            del flat_sample
            del c
            del fig0
            del fig1
            del fig2
            gc.collect()

    return summary


def SZ_mcmc(method, paths, nwalkers='Read', nsteps='Read',
                 amputate=0.5, pics=False):
    #data
    path_to_SZ_data, path_to_startSZ, path_to_prior = paths
    x, y, yerr, z = SZ_data_read(path_to_SZ_data)

    nwalkers_r, nsteps_r, init, prior_box = startSZ_read(path_to_startSZ)
    if nwalkers == 'Read': nwalkers = nwalkers_r
    if nsteps == 'Read': nsteps = nsteps_r

    prior_gauss = prior_read(path_to_prior)
    prior_data = dict()
    prior_data['box'] = prior_box
    prior_data['gauss'] = prior_gauss


    # method
    if method == 'T0':
        model = gauss_model
    if method == 'Tz':
        model = alt_model
        prior_data['const'] = {0: z}

    # mcmc analyze
    mcmc_data = init, x, y, yerr, prior_data
    summary = mcmc_analyze(model, mcmc_data, nwalkers, nsteps, pics=pics)
    summary = summary['0']
    
    # summary
    T0_fit  = summary[1]
    if summary[0] is None or summary[2] is None:
        T0_sp, T0_sm = None, None
    else:
        T0_sp, T0_sm = (summary[2] - T0_fit), (T0_fit - summary[0])

    return z, T0_fit, T0_sp, T0_sm

