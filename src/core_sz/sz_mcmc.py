# -*- coding: utf-8 -*-
"""
Mcmc: analyze a cluster
"""

# Imports
import numpy as np

import paths # manage path imports by paths.py
from core_sz.data import read_SZdata, read_prior, read_startSZ
from core_mcmc.mcmc_quick import mcmc_quick


# FIXME: params_names = None
def sz_analyze(model, input_paths, amputate_percent=0.3, show=False, pics=False):
    # Paths
    path_SZdata, path_prior, path_startSZ = input_paths

    # Read
    x, y, yerr, z = read_SZdata(path_SZdata)
    prior_gauss = read_prior(path_prior)
    nwalkers, nsteps, init, prior_box = read_startSZ(path_startSZ)

    amputate = int(amputate_percent * nsteps)

    # Priors
    prior_data = dict()
    prior_data['box'] = prior_box
    prior_data['gauss'] = prior_gauss
    prior_data['const'] = {0: z}

    # mcmc analyze
    data_params = {
            'x' : x, 
            'y' : y, 
            'yerr' : yerr
            }
    model_params = {
            'model' : model, 
            'init' : init, 
            'prior_data' : prior_data, 
            'params_names' : None
            }
    settings_params = {
            'nwalkers' : nwalkers, 
            'nsteps' : nsteps, 
            'amputate' : amputate}

    summary = mcmc_quick(
            data_params, model_params, settings_params, 
            prnt=True, show=show, save=pics
            )
    
    # summary
    T_summary = summary['0']
    T_left, T_fit, T_right  = T_summary
    T_sp = T_right - T_fit if T_right is not None else None
    T_sm = T_fit - T_left  if T_left  is not None else None

    return z, T_fit, T_sp, T_sm

