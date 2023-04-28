# -*- coding: utf-8 -*-
"""
Show: example data
l = [s.format(params_names[i], *prior_data['box'][i]) for i in range(ndim)] - не дело
"""

# setting path
import sys
sys.path.append('../core')


import gc
import numpy as np
import matplotlib.pyplot as plt

from data import SZ_data_read, startSZ_read, prior_read


def data_show(path, prnt=False, show=False, save=False):
    path_to_SZ_data = path

    #data
    x, y, yerr, z = SZ_data_read(path_to_SZ_data)

    if prnt:
        # SZ data
        print("\nSZ data:")
        print(" z = {}".format(z))
        N = len(x)
        s = " {:>5}: {:>8}  {:>+8} {:>8}"
        f = lambda x: s.format(*x)
        l = map(f, zip(x, y, *([1, -1] * yerr).T))
        print(*l, sep = '\n')

    # Pic
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.errorbar(x, y, [yerr[:, 1], yerr[:, 0]], label='data',
            capsize=3.5, mew=1.5, fmt='.k', alpha=0.5)
    ax.plot(x, 0 * y, 'k')

    ax.set_xlabel(r"$\nu$, GHz")
    ax.set_ylabel(r"SZ signal $\mu$K")
    ax.legend(frameon=False)

    if save:
        fig.savefig(save, bbox_inches='tight')

        # garved collector
        plt.close('all')
        gc.collect()
    if show:
        plt.show()


if __name__ == "__main__":
    # paths
    path_to_SZ_data = '../data/one/SZ_data.txt'
    path_save = './save.png'

    data_show(path_to_SZ_data, prnt=True, show=True, save=False)
