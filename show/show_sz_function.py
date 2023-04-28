# -*- coding: utf-8 -*-
"""
Show: SZ function
comment: добавить исследование неучитывамых порядков
"""

# setting path
import sys
sys.path.append('../core')


import numpy as np
import matplotlib.pyplot as plt

from model import sz_model
from model import nu_default

if __name__ == "__main__":
    # example
    exampe_params = [2.7255, 6.9, 0.0, 1.4]
    #T0, Te, beta, Tau = exampe_params

    # Picture
    fig, ax = plt.subplots(figsize=(8, 8))
    xline = np.linspace(min(nu_default), max(nu_default), 100)
    yline = sz_model(*exampe_params, xline)
    yline_norel = sz_model(*exampe_params, xline, rel_corrs=False)

    ax.plot(xline, 0 * xline, '-k')
    ax.plot(xline, yline, label='model')
    ax.plot(xline, yline_norel, label='model without relativ corr.')

    ax.set_xlabel(r"$\nu$, GHz")
    ax.set_ylabel(r"SZ signal $\mu$K")
    ax.legend(frameon=False)

    plt.show()
