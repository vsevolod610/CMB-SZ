# -*- coding: utf-8 -*-
"""
это все нужно сделать через логирование
"""

from create_data import N as N_create_data
from start_analyze import N as N_start_analyze

print("\nN: \n", 
      "N from create_data.py = %d\n"%(N_create_data),
      "N from start_analyze.py = %d\n"%(N_start_analyze))


from create_data import model as model_create_data
from start_analyze import method as method_start_analyze
from model import *

if model_create_data is gauss_model:
    modelname_create_data = 'gauss_model'
elif model_create_data is lazy_model:
    modelname_create_data = 'lazy_model'
elif model_create_data is simple_model:
    modelname_create_data = 'simple_model'
else:
    modelname_create_data = 'error'

print("\nModel: \n",
      "model from create_data.py = %s\n"%(modelname_create_data), 
      "method from start_analyze.py = %s\n"%(method_start_analyze))


from start_analyze import nwalkers as nwalkers_start_analyze
from start_analyze import nsteps as nsteps_start_analyze
from data import startSZ_read
from start_analyze import path_startSZ
nwalkers_r, nsteps_r, init, prior_box = startSZ_read(path_startSZ)

print("\n Nwalkers & Nsteps: \n", 
      "nwalkers & nsteps from startSZ = %s, %s\n"%(nwalkers_r, nsteps_r),
      "nwalkers & nsteps from start_analyze.py = %s, %s\n"
      %(nwalkers_start_analyze, nsteps_start_analyze))


print("done checking!")
