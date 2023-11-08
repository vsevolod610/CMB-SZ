#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Show: Spectral transmission functions
"""

# Imports
import numpy as np
import matplotlib.pyplot as plt

import paths # manage path imports by paths.py
from core_sz.trans.trans_data import filtrs, trans_data, trans_gauss, gaussian


# trans function analysis
print("Trans funtions:")
for wave in filtrs:
    nu = trans_data[wave]['nu']
    # FIXME: python 3.12 can make multiline f-string
    print(f"Wave {wave:>3} GHz: nu = [{min(nu):>6.2f} ... {max(nu):>6.2f}] of {len(nu)} dots")

# gauss approximation analysis
print("\n Gauss approximation:")
for wave in filtrs:
    mu = trans_gauss[wave]['mu']
    sigma = np.sqrt(trans_gauss[wave]['sigma2'])
    print(f"Wave {wave:>3} GHz: mu = {mu:>6.2f}, sigma = {sigma:>5.2f}")


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

    ax.plot(xline, gaussian(xline, mu, sigma2), 'k')

ax.set_xlabel(r'$\nu$ GHz')
ax.legend(frameon=False)
fig.savefig('pics/trans_function.png', bbox_inches='tight')
plt.show()

