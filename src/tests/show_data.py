#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Show: data from SZ_data.txt
"""

# Imports
import gc
import numpy as np
import matplotlib.pyplot as plt

import paths # manage path imports by paths.py
from core_sz.data import read_SZdata


def data_show(path, prnt=False, show=False, save=False):
    # Data
    x, y, yerr, z = read_SZdata(path)

    # Print
    if prnt:
        # SZ data
        print(f"SZ data: {z = }")
        print(f" {'nu':>5}  {'y':>8}  {'yerr+':>8} {'yerr-':>8}")
        for xi, yi, yerrpi, yerrmi in zip(x, y, yerr[:, 0], -1 * yerr[:, 1]):
            # FIXME:  {yerrpi:>+8} if you wnat
            print(f" {xi:>5}: {yi:>8}  {yerrpi:>8} {yerrmi:>8}")

    # Pictures
    if show or save:
        fig, ax = plt.subplots(figsize=(8, 8))
        ax.axhline(y=0, color='k')
        ax.errorbar(x, y, [yerr[:, 1], yerr[:, 0]], label='data',
                    capsize=3.5, mew=1.5, fmt='.k', alpha=0.5)

        ax.set_xlabel(r"$\nu$, GHz")
        ax.set_ylabel(r"SZ signal $\mu$K")
        ax.legend(frameon=False)
        if save: fig.savefig(save, bbox_inches='tight')
    if show: plt.show()

    # garved collector
    plt.close('all')
    gc.collect()


if __name__ == "__main__":
    path_SZdata = paths.src + '/../data/example_cluster/SZ_data.txt'
    path_save = 'pics/data.png'
    data_show(path_SZdata, prnt=True, show=True, save=path_save)
