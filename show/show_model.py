# -*- coding: utf-8 -*-
"""
Show: models
"""
import sys

# setting path
sys.path.append('../core')


import numpy as np
import matplotlib.pyplot as plt 
from scipy import integrate

from show_data import x, y, yerr
from trans_function import filtrs
from model import sz_model, trans_model, gauss_model


if __name__ == "__main__":
    # example
    exampe_params = [2.7255, 6.9, 0.0, 1.4]
    #T0, Te, beta, Tau = exampe_params

    # models analisys
    sz = sz_model(*exampe_params, x)
    trans_sz = trans_model(*exampe_params)
    gauss_sz = gauss_model(*exampe_params)
    print("\nWaves = {} GHz".format(filtrs))
    print("only clear SZfunction:    ", *[" {:5.0f}".format(i) for i in sz])
    print("with clear trans-function:", *[" {:5.0f}".format(i) for i in trans_sz])
    print("with gauss trans-function:", *[" {:5.0f}".format(i) for i in gauss_sz])

    # Picture:
    fig, ax = plt.subplots(figsize=(8, 8))
    xline = np.linspace(min(x), max(x), 100)
    yline = sz_model(*exampe_params, xline)

    ax.plot(x, 0 * x, 'k')
    ax.errorbar(x, y, yerr.T, label='data', 
            capsize=3.5, mew=1.5, fmt='.k', alpha=0.5)
    ax.plot(xline, yline, label='model')
    ax.plot(x, trans_sz, 'o', label='with clear trans-function')
    ax.plot(x, gauss_sz, 'o', label='with gauss trans-function')

    ax.set_xlabel(r"$\nu$, GHz")
    ax.set_ylabel(r"SZ signal $\mu$K")
    ax.legend(frameon=False)

    plt.show()
