# -*- coding: utf-8 -*-
"""
Analyze clusters
"""

import numpy as np
import matplotlib.pyplot as plt
from chainconsumer import ChainConsumer

from model import gauss_model as model
from mcmc_kern import mcmc_kern, pic_chain, pic_fit
from data import SZ_data_read, startSZ_read, prior_read


def mcmc_analyze(paths, nwalkers='Read', nsteps='Read', amputate=0.5,
                                                        pic=False):
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

    amputete_step = int(amputate * nsteps)
    ndim = len(init)

    # mcmc analyze
    sampler = mcmc_kern(model, nwalkers, nsteps, 
                        init, x, y, yerr, prior_data)
    flat_sample = sampler.chain[:, amputete_step : , :].reshape((-1, ndim))
    c = ChainConsumer()
    c.add_chain(flat_sample)
    summary = c.analysis.get_summary()['0']
    
    T0_fit  = summary[1]
    if summary[0] is None or summary[2] is None:
        T0_sp, T0_sm = None, None
    else:
        T0_sp, T0_sm = (summary[2] - T0_fit), (T0_fit - summary[0])

    if pic:
        # pics
        fig0, ax = pic_chain(sampler, params_names=None)
        fig1 = c.plotter.plot(display=False, legend=False, figsize=(6, 6))
        fig2, ax = pic_fit(sampler, model, x, y, yerr, prior_data)
        if type(pic) is list:
            fig0.savefig(pic[0])
            fig1.savefig(pic[1])
            fig2.savefig(pic[2])

    return z, T0_fit, T0_sp, T0_sm

