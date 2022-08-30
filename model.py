# -*- coding: utf-8 -*-
"""
SZfunction with trans-function
"""


import numpy as np
import matplotlib.pyplot as plt 
from scipy import integrate

from data import x, y, yerr
from trans_function import filtrs, spec_trans, gauss_trans, gauss
from SZfunction import SZfunction


def SZmodel(T0, Te, beta, Tau, nu):
    return [SZfunction(T0, Te, beta, Tau, nu_i) for nu_i in nu]


def spec_model(T0, Te, beta, Tau, nu=0):
    m = []
    for wave in filtrs:
        nu = spec_trans[wave]['nu']
        tr = spec_trans[wave]['tr']
        sz = np.array([SZfunction(T0, Te, beta, Tau, nu_i) for nu_i in nu])
        m.append(np.trapz(sz * tr, nu))
    return np.array(m)


def gauss_model(T0, Te, beta, Tau, nu=0):
    m = []
    for wave in filtrs:
        mu = gauss_trans[wave]['mu']
        sigma2 = gauss_trans[wave]['sigma2']
        f = lambda nu: SZfunction(T0, Te, beta, Tau, nu) * gauss(nu, mu, sigma2)
        nu_start = mu - 4 * np.sqrt(sigma2)
        nu_stop = mu + 4 * np.sqrt(sigma2)
        m.append(integrate.quad(f, nu_start, nu_stop, epsabs=1e+3,epsrel=1e-2, limit=2)[0])
    return np.array(m)


if __name__ == "__main__":
    # example
    exampe_params = [2.7255, 6.9, 0.0, 1.4]
    #T0, Te, beta, Tau = exampe_params

    # models analisys
    sz = SZmodel(*exampe_params, x)
    spec_sz = spec_model(*exampe_params)
    gauss_sz = gauss_model(*exampe_params)
    print("\nWaves = {} GHz".format(filtrs))
    print("only clear SZfunction:", *[" {:.0f}".format(i) for i in sz])
    print("with clear trans-function:", *[" {:.0f}".format(i) for i in spec_sz])
    print("with gauss trans-function:", *[" {:.0f}".format(i) for i in gauss_sz])

    # Picture:
    fig, ax = plt.subplots(figsize=(8, 8))
    xline = np.linspace(min(x), max(x), 100)
    yline = SZmodel(*exampe_params, xline)

    ax.plot(x, 0 * x, 'k')
    ax.errorbar(x, y, yerr.T, capsize=3.5, mew=1.5, fmt='.k', alpha=0.5, label='data')
    ax.plot(xline, yline, label='model')
    ax.plot(x, spec_sz, 'o', label='with clear trans-function')
    ax.plot(x, gauss_sz, 'o', label='with gauss trans-function')

    ax.set_xlabel(r"$\nu$, GHz")
    ax.set_ylabel(r"SZ signal $\mu$K")
    ax.legend(frameon=False)

    plt.show()
