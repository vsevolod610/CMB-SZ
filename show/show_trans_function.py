# -*- coding: utf-8 -*-
"""
Show: Spectral transmission functions
"""

# setting path
import sys
sys.path.append('../core')

import numpy as np
import matplotlib.pyplot as plt

from trans_function import filtrs, trans_data, trans_gauss, gauss


if __name__ == "__main__":
    # trans function analysis
    print("\n trans funtions: ")
    for wave in filtrs:
        nu = trans_data[wave]['nu']
        print(
            "Wave {:>3} GHz: nu = [{:>6.2f} ... {:>6.2f}] of {} dots"
            .format(wave, min(nu), max(nu), len(nu)))

    # gauss approximation analysis
    print("\n gauss approximation")
    for wave in filtrs:
        mu = trans_gauss[wave]['mu']
        sigma = np.sqrt(trans_gauss[wave]['sigma2'])
        print(
            "Wave {:>3} GHz: mu = {:>6.2f}, sigma = {:>5.2f}"
            .format(wave, mu, sigma))

    # pictures with gauss approximation
    fig, ax = plt.subplots(figsize=(8, 8))
    for wave in filtrs:
        # trans function data
        x = trans_data[wave]['nu']
        y = trans_data[wave]['tr']
        ax.plot(x, y, label='{} GHz'.format(wave))

        # trans function gauss approximation
        mu = trans_gauss[wave]['mu']
        sigma2 = trans_gauss[wave]['sigma2']
        xline = np.linspace(min(x) - 10, max(x) + 10, 100)
        ax.plot(xline, gauss(xline, mu, sigma2), 'k')

        ax.set_xlabel(r'$\nu$ GHz')
        ax.legend(frameon=False)

    plt.show()
