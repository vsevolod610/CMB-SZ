#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Read .fits files and set spectral-transmission functions
"""

# Imports: built-in
import numpy as np
from scipy import integrate

# Imports: from project
from .trans_data import filtrs, trans_data, trans_gauss, gaussian


def sptr_integrate_off(func):
    nu = np.array(filtrs, dtype=float)
    signal = func(nu)
    return signal

def sptr_integrate_row(func):
    signal = []
    for wave in filtrs:
        nu = trans_data[wave]['nu']
        tr = trans_data[wave]['tr']
        sz = np.trapz(func(nu) * tr, nu)
        signal.append(sz)
    return np.array(signal)

def sptr_integrate_gauss(func):
    signal = []
    for wave in filtrs:
        mu = trans_gauss[wave]['mu']
        sigma2 = trans_gauss[wave]['sigma2']

        f = lambda nu: func(nu) * gaussian(nu, mu, sigma2)

        nu_start = mu - 4 * np.sqrt(sigma2)
        nu_stop = mu + 4 * np.sqrt(sigma2)
        sz, _ = integrate.quad(f, nu_start, nu_stop, 
                               epsabs=1e+3, epsrel=1e-2, limit=2)
        signal.append(sz)
    return np.array(signal)
 
    
# Account tranmsission functions: mode: 'off', 'row', 'gauss'
# NB! always return array, not number
def trans_integrate(func, mode='off'):
    if mode == 'off': return sptr_integrate_off(func)
    elif mode == 'row': return sptr_integrate_row(func)
    elif mode == 'gauss': return sptr_integrate_gauss(func)
    else: raise NameError(f"{mode = } not in ('off', 'row', 'gauss')")


