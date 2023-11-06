#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Show: SZ function
"""

# Imports
import numpy as np
import matplotlib.pyplot as plt

import paths # manage path imports by paths.py
from core_sz.sz_signal import sz_effect


# example params: T0, Te, beta, Tau, Tz
exampe_params = [2.7255, 6.9, 0.0, 1.4, 2.7255]
nu_min, nu_max = 70, 353


# Picture
fig, ax = plt.subplots(figsize=(8, 8))
xline = np.linspace(nu_min, nu_max, 100)
yline0 = sz_effect(*exampe_params, xline, corrs=False)
yline1 = sz_effect(*exampe_params, xline, corrs=True)

ax.axhline(y=0, color='k')
ax.plot(xline, yline0, label='SZ effect (no corr)')
ax.plot(xline, yline1, label='SZ effect')

ax.set_xlabel(r"$\nu$, GHz")
ax.set_ylabel(r"SZ signal $\mu$K")
ax.legend(frameon=False)

fig.savefig('pics/sz_effect.png', bbox_inches='tight')
plt.show()
