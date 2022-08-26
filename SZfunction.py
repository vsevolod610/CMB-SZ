# -*- coding: utf-8 -*-
"""
Set SZfunction
"""


import numpy as np
import matplotlib.pyplot as plt

from data import x, y, yerr

h = 2 * 3.141 * 1.05e-27
c = 3e10
k = 1.38e-16
m = 0.91e-27
coef1 = (h * 1e9) / k
coef2 = (1e3 * 1.6e-12) / (m * c ** 2)


def SZfunction(T0, kTe, beta, Tau, nu, rel_corrs=True):
    x = coef1 * nu / T0
    theta = coef2 * kTe

    X = x * (np.exp(x) + 1) / (np.exp(x) - 1)
    S = x / np.sinh(x / 2)

    Y0 = X - 4
    Y1 = (- 10 + 47 / 2 * X - 42 / 5 * X ** 2 + 7 / 10 * X ** 3 
            + S ** 2 * (- 21 / 5 + 7 / 5 * X))
    Y2 = (- 15 / 2 + 1023 / 8 * X - 868 / 5 * X ** 2 + 329 / 5 * X ** 3 
            - 44 / 5 * X ** 4 + 11 / 30 * X ** 5 
            + S ** 2 * (- 434 / 5 + 658 / 5 * X - 242 / 5 * X ** 2 + 143 / 30 * X ** 3) 
            + S ** 4 * (- 44 / 5 + 187 / 60 * X))
    Y3 = (15 / 2 + 2505 / 8 * X - 7098 / 5 * X ** 2 + 14253 / 10 * X ** 3 
            - 18594 / 35 * X ** 4 + 12059 / 140 * X ** 5 - 128 / 21 * X ** 6 
            + 16 / 105 * X ** 7 
            + S ** 2 * (- 7098 / 10 + 14253 / 5 * X - 102267 / 35 * X ** 2 
                    + 156767 / 140 * X ** 3 - 1216 / 7 * X ** 4 + 64 / 7 * X ** 5) 
            + S ** 4 * (- 18594 / 35 + 205003 / 280 * X - 1920 / 7 * X ** 2 
                    + 1024 / 35 * X ** 3) 
            + S ** 6 * (- 544 / 21 + 992 / 105 * X))
    Y4 = (- 135 / 32 + 30375 / 128 * X - 62391 / 10 * X ** 2 + 614727 / 40 * X ** 3 
            - 124389 / 10 * X ** 4 + 355703 / 80 * X ** 5 - 16568 / 21 * X ** 6 
            + 7516 / 105 * X ** 7 - 22 / 7 * X ** 8 + 11 / 210 * X ** 9 
            + S ** 2 * (- 62391 / 20 + 614727 / 20 * X - 1368279 / 20 * X ** 2 
                    + 4624139 / 80 * X ** 3 - 157396 / 7 * X ** 4 + 30064 / 7 * X ** 5
                    - 2717 / 7 * X ** 6 + 2761 / 210 * X ** 7) 
            + S ** 4 * (- 124389 / 10 + 6046951 / 160 * X - 248520 / 7 * X ** 2 
                    + 481024 / 35 * X ** 3 - 15972 / 7 * X ** 4 + 18689 / 140 * X ** 5) 
            + S ** 6 * (- 70414 / 21 + 465992 / 105 * X - 11792 / 7 * X ** 2 
                    + 19778 / 105 * X ** 3) 
            + S ** 8 * (- 682 / 7 + 7601 / 210 * X))
    C1 = 10 - 47 / 5 * X + 7 / 5 * X ** 2 + 7 / 10 * S ** 2
    C2 = (25 - 1117 / 10 * X + 847 / 10 * X ** 2 - 183 / 10 * X ** 3 + 11 / 10 * X ** 4 
            + S ** 2 * (847 / 20 - 183 / 5 * X + 121 / 20 * X ** 2) 
            + 11 / 10 * S ** 4)
    D0 = - 2 / 3 + 11 / 30 * X
    D1 = - 4 + 12 * X - 6 * X ** 2 + 19 / 30 * X ** 3 + S ** 2 * (- 3 + 19 / 15 * X)

    if rel_corrs:
        R = theta ** 2 * Y1 + theta ** 3 * Y2 + theta ** 4 * Y3 + theta ** 5 * Y4 + \
            + beta ** 2 * (1 / 3 * Y0 + theta * (5 / 6 * Y0 + 2 / 3 * Y1)) + \
            + beta * (1 + theta * C1 + theta ** 2 * C2) + \
            + beta ** 2 * (D0 + theta * D1)
    else:
        R = 0

    return 1e+6 * T0 * Tau * (theta * Y0 + R)


if __name__ == "__main__":
    # example
    exampe_params = [2.7255, 6.9, 0.0, 1.4]
    #T0, Te, beta, Tau = exampe_params

    # Picture
    fig, ax = plt.subplots(figsize=(8, 8))
    xline = np.linspace(min(x), max(x), 100)
    yline = [SZfunction(*exampe_params, xi) for xi in xline]
    yline_norel = [SZfunction(*exampe_params, xi, rel_corrs=False) for xi in xline]

    ax.plot(x, 0 * x, 'k')
    ax.errorbar(x, y, yerr.T, capsize=3.5, mew=1.5, fmt='.k', alpha=0.5, label='data')
    ax.plot(xline, yline, label='model')
    ax.plot(xline, yline_norel, label='model without relativ corr.')

    ax.set_xlabel(r"$\nu$, GHz")
    ax.set_ylabel(r"SZ signal $\mu$K")
    ax.legend(frameon=False)

    plt.show()


