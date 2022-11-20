# -*- coding: utf-8 -*-
"""
SZfunction with trans-function
comment: x=None сомнительная практика
"""


import numpy as np
import matplotlib.pyplot as plt 
from scipy import integrate

from trans_function import filtrs, trans_data, trans_gauss, gauss
from sz_function import SZfunction

h = 2 * 3.141 * 1.05e-27
c = 3e10
k = 1.38e-16
m = 0.91e-27
coef1 = (h * 1e9) / k
coef2 = (1e3 * 1.6e-12) / (m * c ** 2)


#NB! all model return array, not number

def sz_model(T0, Te, beta, Tau, nu, rel_corrs=True):
    theta = coef2 * Te
    x = coef1 * nu / T0
    sz = [SZfunction(T0, theta, beta, Tau, x_i, rel_corrs) for x_i in x]
    return np.array(sz)


def trans_model(T0, Te, beta, Tau, x=None):
    theta = coef2 * Te
    m = []
    for wave in filtrs:
        nu = trans_data[wave]['nu']
        tr = trans_data[wave]['tr']
        x = coef1 * nu / T0
        sz = np.array([SZfunction(T0, theta, beta, Tau, x_i) for x_i in x])
        m.append(np.trapz(sz * tr, nu))
    return np.array(m)


def gauss_model(T0, Te, beta, Tau, x=None):
    theta = coef2 * Te
    m = []
    for wave in filtrs:
        mu = trans_gauss[wave]['mu']
        sigma2 = trans_gauss[wave]['sigma2']
        f = lambda nu: (SZfunction(T0, theta, beta, Tau, coef1 * nu / T0) 
                        * gauss(nu, mu, sigma2))
        nu_start = mu - 4 * np.sqrt(sigma2)
        nu_stop = mu + 4 * np.sqrt(sigma2)
        m.append(integrate.quad(f, nu_start, nu_stop, 
                                epsabs=1e+3,epsrel=1e-2, limit=2)[0])
    return np.array(m)

def alt_model(Tz, Te, beta, Tau, z, x=None):
    T0 = 2.7255
    theta = coef2 * Te
    m = [] 
    for wave in filtrs:
        mu = trans_gauss[wave]['mu']
        sigma2 = trans_gauss[wave]['sigma2']
        f = lambda nu: (SZfunction(T0, theta, beta, Tau, coef1 * nu / Tz * (1 + z)) * gauss(nu, mu, sigma2))
        nu_start = mu - 4 * np.sqrt(sigma2)
        nu_stop = mu + 4 * np.sqrt(sigma2)
        m.append(integrate.quad(f, nu_start, nu_stop, 
                                epsabs=1e+3,epsrel=1e-2, limit=2)[0])
    return np.array(m)