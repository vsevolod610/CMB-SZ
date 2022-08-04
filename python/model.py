# -*- coding: utf-8 -*-
"""
SZfunction with trans-function
"""


import numpy as np
import matplotlib.pyplot as plt 
from scipy import integrate

from trans_function import filtrs, spec_trans, gauss_trans, gauss
from SZfunction import SZfunction


def model(T0, Te, beta, Tau, nu):
    return [SZfunction(T0, Te, beta, Tau, nu_i) for nu_i in nu]


def spec_model(T0, Te, beta, Tau):
    m = []
    for wave in filtrs:
        nu = spec_trans[wave]['nu']
        tr = spec_trans[wave]['tr']
        sz = np.array([SZfunction(T0, Te, beta, Tau, nu_i) for nu_i in nu])
        m.append(np.trapz(sz * tr, nu))
    return np.array(m)


def gauss_model(T0, Te, beta, Tau):
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
    example_data = np.array([
            [-8.44e+04,  1.90e+04,  1.90e+04],
            [-7.05e+04,  5.70e+03,  6.33e+03],
            [-4.45e+04,  4.43e+03,  4.43e+03],
            [-6.47e+03,  5.70e+03,  6.97e+03],
            [ 1.01e+05,  3.04e+04,  3.04e+04],])
    x = np.array([70, 100, 143, 217, 353], dtype=float)
    y = example_data[:, 0]
    yerr = example_data[:, [1, 2]]
    exampe_params = [2.7255, 6.9, 0.0, 1.4]
    #T0, Te, beta, Tau = exampe_params

    # models analisys
    sz = model(*exampe_params, x)
    spec_sz = spec_model(*exampe_params)
    gauss_sz = gauss_model(*exampe_params)
    print("Waves = {} GHz".format(filtrs))
    print("only clear SZfunction:", *[" {:.0f}".format(i) for i in sz])
    print("with clear trans-function:", *[" {:.0f}".format(i) for i in spec_sz])
    print("with gauss trans-function:", *[" {:.0f}".format(i) for i in gauss_sz])

    # Picture:
    fig, ax = plt.subplots(figsize=(8, 8))
    xline = np.linspace(min(x), max(x), 100)
    yline = model(*exampe_params, xline)

    ax.plot(x, 0 * x, 'k')
    ax.errorbar(x, y, yerr.T, capsize=3.5, mew=1.5, fmt='.k', alpha=0.5, label='data')
    ax.plot(xline, yline, label='model')
    ax.plot(x, spec_sz, 'o', label='with clear trans-function')
    ax.plot(x, gauss_sz, 'o', label='with gauss trans-function')

    ax.set_xlabel(r"$\nu$, GHz")
    ax.set_ylabel(r"SZ signal $\mu$K")
    ax.legend(frameon=False)

    plt.show()
