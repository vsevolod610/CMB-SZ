# -*- coding: utf-8 -*-
"""
set of T_i --> Estimate T
"""

# Imports
import numpy as np

import paths # manage path imports by paths.py
from core_mcmc.mcmc_quick import mcmc_quick
from start_analyze import read_result, path_result
from config import config

# Paths
path_pics = [paths.src + '/main/output/res/pic_chian.png',
             paths.src + '/main/output/res/pic_params.png',
             paths.src + '/main/output/res/pic_fit.png']

# Data
x, y, yerr = read_result(path_result)

# Model
modelT = config['modelT']

# MCMC 
data_params = {
        'x' : x, 
        'y' : y, 
        'yerr' : yerr
        }
model_params = {
        'model' : modelT, 
        'init' : np.array([[2.7255], [0.2]]).T, 
        'prior_data' : None,
        'params_names' : None
        }
settings_params = {
        'nwalkers' : 200, 
        'nsteps' : 200, 
        'amputate' : 50}

summary = mcmc_quick(
        data_params, model_params, settings_params, 
        prnt=False, show=False, save=path_pics
        )

Tl, T, Tr = summary['0']
print(f"T = {T} +{Tr - T} -{T - Tl}")
