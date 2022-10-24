# -*- coding: utf-8 -*-
"""
Analyze clusters
comment: нодо бы перенести этот файл в корень
"""

import numpy as np
import matplotlib.pyplot as plt
from chainconsumer import ChainConsumer

from model import gauss_model as model
from mcmc_kern import mcmc_kern, pic_chain, pic_fit
from data import SZ_data_read, startSZ_read, prior_read


def mcmc_analyze(paths, nwalkers_force=None, nsteps_force=None, pic=False):
    #data
    path_to_SZ_data, path_to_startSZ, path_to_prior = paths
    x, y, yerr, z = SZ_data_read(path_to_SZ_data)

    nsteps, nwalkers, init, prior_box = startSZ_read(path_to_startSZ)
    if nwalkers_force: nwalkers = nwalkers_force
    if nsteps_force: nsteps = nsteps_force

    prior_gauss = prior_read(path_to_prior)
    prior_data = dict()
    prior_data['box'] = prior_box
    prior_data['gauss'] = prior_gauss

    amputete = int(0.2 * nsteps)
    ndim = len(init)

    # mcmc analyze
    sampler = mcmc_kern(model, nwalkers, nsteps, 
                        init, x, y, yerr, prior_data)
    flat_sample = sampler.chain[:, amputete : , :].reshape((-1, ndim))
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
        fig, ax = pic_chain(sampler, params_names=None)
        fig = c.plotter.plot(display=False, legend=False, figsize=(6, 6))
        fig, ax = pic_fit(sampler, model, x, y, yerr, prior_data)

    return z, T0_fit, T0_sp, T0_sm

