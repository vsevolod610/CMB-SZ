#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Intergrate spectral-transmission function with fuction of SZ-effect
Define models to analyze
"""

# Imports: built-in
import numpy as np
from scipy import integrate

# Imports: from project
from .trans_function import filtrs, trans_data, trans_gauss, gaussian
from .sz_effect import sz_effect


# Account tranmsission functions: mode: 'off', 'row', 'gauss'
# NB! always return array, not number
def sz_signal(T0, Te, beta, Tau, Tz, z=0, corrs=True, mode='guass'):

    if mode == 'off':
        nu = np.array(filtrs, dtype=float)
        signal = sz_effect(T0, Te, beta, Tau, Tz, nu, z=z, corrs=corrs) 
    elif mode == 'row':
        signal = []
        for wave in filtrs:
            nu = trans_data[wave]['nu']
            tr = trans_data[wave]['tr']
            sz = sz_effect(T0, Te, beta, Tau, Tz, nu, z=z, corrs=corrs) 
            signal.append(np.trapz(sz * tr, nu))
        signal = np.array(signal)
    elif mode == 'gauss':
        signal = []
        for wave in filtrs:
            mu = trans_gauss[wave]['mu']
            sigma2 = trans_gauss[wave]['sigma2']

            def f(nu):
                sz = sz_effect(T0, Te, beta, Tau, Tz, nu, z=z, corrs=corrs)
                tr = gaussian(nu, mu, sigma2)
                return sz * tr

            nu_start = mu - 4 * np.sqrt(sigma2)
            nu_stop = mu + 4 * np.sqrt(sigma2)
            sz, _ = integrate.quad(f, nu_start, nu_stop, 
                                   epsabs=1e+3, epsrel=1e-2, limit=2)
            signal.append(sz)
        signal = np.array(signal)

    return signal


### Models
# model(*params, *const, x)

# ------------- test only -----------------

# very light test model
def model_000(T0, Te, z=0):
    beta = 0.00
    Tau = 1.0
    Tz = T0 * (1 + z)
    return sz_signal(T0, Te, beta, Tau, Tz, corrs=False, mode='off')

# ------------- default -------------------

# light test model
def model_100(T0, Te, beta, Tau, z=0, x=None):
    Tz = T0 * (1 + z)
    return sz_signal(T0, Te, beta, Tau, T0, corrs=False, mode='off')

# test model
def model_110(T0, Te, beta, Tau, z=0, x=None):
    Tz = T0 * (1 + z)
    return sz_signal(T0, Te, beta, Tau, Tz, corrs=True, mode='off')

# heavy model
def model_111(T0, Te, beta, Tau, z=0, x=None):
    Tz = T0 * (1 + z)
    return sz_signal(T0, Te, beta, Tau, Tz, corrs=True, mode='row')

# default model
def model_112(T0, Te, beta, Tau, z=0, x=None): 
    Tz = T0 * (1 + z)
    return sz_signal(T0, Te, beta, Tau, Tz, corrs=True, mode='gauss')

# ----------- alternative -----------------

# default alternative model
def model_212(Tz, Te, beta, Tau, z, x=None):
    T0 = 2.7255
    return sz_signal(T0, Te, beta, Tau, Tz, z=z, corrs=True, mode='gauss')

# test alternative model
def model_210(Tz, Te, beta, Tau, z, x=None):
    T0 = 2.7255
    return sz_signal(T0, Te, beta, Tau, Tz, z=z, corrs=True, mode='off')
