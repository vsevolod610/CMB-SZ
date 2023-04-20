# -*- coding: utf-8 -*-
"""
Analyze cluster
"""

import gc
import numpy as np
from chainconsumer import ChainConsumer

from model import *
from data import SZ_data_read, startSZ_read, prior_read
from mcmc import mcmc


def sz_analyze(method, paths, nwalkers='Read', nsteps='Read',
                 amputate=0.3, pics=False):

    #data
    path_to_SZ_data, path_to_startSZ, path_to_prior = paths

    x, y, yerr, z = SZ_data_read(path_to_SZ_data)
    nwalkers_r, nsteps_r, init, prior_box = startSZ_read(path_to_startSZ)
    if nwalkers == 'Read': nwalkers = nwalkers_r
    if nsteps == 'Read': nsteps = nsteps_r
    amputate_step = int(amputate * nsteps)

    prior_data = dict()
    prior_gauss = prior_read(path_to_prior)
    prior_data['box'] = prior_box
    prior_data['gauss'] = prior_gauss

    # method
    if method == 'T0':
        model = gauss_model
    if method == 'lazy':
        model = lazy_model
    if method == 'simple':
        model = simple_model
    if method == 'Tz':
        model = alt_model
        prior_data['const'] = {0: z}

    # mcmc analyze
    summary = mcmc(data=(x, y, yerr), 
                   model_params=(model, init, prior_data, None), 
                   settings=(nwalkers, nsteps, amputate_step), 
                   prnt=True, show=False, save=pics)
    summary = summary['0']
    
    # summary
    T0_fit  = summary[1]
    if summary[0] is None or summary[2] is None:
        T0_sp, T0_sm = None, None
    else:
        T0_sp, T0_sm = (summary[2] - T0_fit), (T0_fit - summary[0])

    return z, T0_fit, T0_sp, T0_sm

