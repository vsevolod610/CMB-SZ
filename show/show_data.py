# -*- coding: utf-8 -*-
"""
Show: example data
l = [s.format(params_names[i], *prior_data['box'][i]) for i in range(ndim)] - не дело
"""

# setting path
import sys
sys.path.append('../core')


import numpy as np
import matplotlib.pyplot as plt

from data import SZ_data_read, startSZ_read, prior_read


def data_show(paths, prnt=False, show=False, save=False):
    path_to_SZ_data, path_to_startSZ, path_to_prior = paths

    #data
    x, y, yerr, z = SZ_data_read(path_to_SZ_data)
    nwalkers, nsteps, init, prior_box = startSZ_read(path_to_startSZ)
    prior_gauss = prior_read(path_to_prior)
    prior_data = dict()
    prior_data['box'] = prior_box
    prior_data['gauss'] = prior_gauss

    ndim = len(init)
    params_names = list(map(str, np.arange(ndim)))

    if prnt:
        # SZ data
        print("\nSZ data:")
        print(" z = {}".format(z))
        N = len(x)
        s = " {:>5}: {:>8}  {:>+8} {:>8}"
        f = lambda x: s.format(*x)
        l = map(f, zip(x, y, *([1, -1] * yerr).T))
        print(*l, sep = '\n')

        # Prior
        print("\nPriors")
        print("  box:")
        s = "    {:<4}: [{:^6}, {:>5}]"
        l = [s.format(params_names[i], *prior_data['box'][i]) for i in range(ndim)]
        print(*l, sep='\n')

        print("  gauss:")
        names = list(prior_data['gauss'].keys())
        s = "    {:<4}: {} {:>+3} -{:>3}"
        l = [s.format(name, *prior_data['gauss'][name]) for name in names]
        print(*l, sep='\n')

        # Start
        print("\nStart:")
        s = "  {:<4}: {:>4} +- {:<5}"
        l = [s.format(params_names[i], *init[i]) for i in range(ndim)]
        print(*l, sep='\n')

    # Pic
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.errorbar(x, y, yerr.T, label='data',
            capsize=3.5, mew=1.5, fmt='.k', alpha=0.5)
    ax.plot(x, 0 * y, 'k')

    ax.set_xlabel(r"$\nu$, GHz")
    ax.set_ylabel(r"SZ signal $\mu$K")
    ax.legend(frameon=False)

    if save:
        fig.savefig(save, bbox_inches='tight')
    if show:
        plt.show()


if __name__ == "__main__":
    # paths
    path_to_SZ_data = '../data/SZ_data.txt'
    path_to_startSZ = '../data/startSZ.sss'
    path_to_prior = '../data/prior.dat'

    paths = [path_to_SZ_data, path_to_startSZ, path_to_prior]
    path_save = './save.png'

    data_show(paths, prnt=True, show=False, save=path_save)
