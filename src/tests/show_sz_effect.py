#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Show: SZ function
"""

# Imports
import numpy as np
import matplotlib.pyplot as plt

import paths # manage path imports by paths.py
from core_sz.sz_signal import sz_effect, sz_full_effect


# example params: T0, Te, beta, Tau, Tz
exampe_params = [2.7255, 6.9, 0.0001, 1.4, 2.7255]
nu_min, nu_max = 70, 353

def show(yline, origin, comment):
    y_prnt = f"[{yline[0]:5.0f} ... {yline[-1]:5.0f}]"
    delta_y = np.abs((yline[-1] - origin[-1]) / origin[-1])
    if not (yline is origin):
        print(f"{comment:<30}: {y_prnt} <- delta = {delta_y:.2%}")
    else:
        print(f"{comment:<30}: {y_prnt}")


# Picture
fig, ax = plt.subplots(figsize=(8, 8))
xline = np.linspace(nu_min, nu_max, 100)
yline0 = sz_effect(xline, *exampe_params, corrs=False)
yline1 = sz_effect(xline, *exampe_params, corrs=True)

# print
show(yline0, yline1, comment="SZ effect (no corr)")
show(yline1, yline1, comment="SZ effect")

ax.axhline(y=0, color='k')
ax.plot(xline, yline0, label='SZ effect (no corr)')
ax.plot(xline, yline1, label='SZ effect')

# SZ full-effect
if input("Do You Want To evaluate SZ full-effect? [y/n] ") == "y": 
    mu2_0, mu2_1 = 1, -1

    # data
    xline_full = np.linspace(nu_min, nu_max, 5)
    yline2_0 = sz_full_effect(xline_full, *exampe_params, mu=mu2_0)
    yline2_1 = sz_full_effect(xline_full, *exampe_params, mu=mu2_1)

    # Pictures
    ax.plot(xline_full, yline2_0, label=f"SZ full-effect (mu = {mu2_0})")
    ax.plot(xline_full, yline2_1, label=f"SZ full-effect (mu = {mu2_1})")

    # print
    yline_origin = sz_effect(xline_full, *exampe_params, corrs=True)
    show(yline2_0, yline_origin, comment=f"SZ full-effect (mu = {mu2_0})")
    show(yline2_1, yline_origin, comment=f"SZ full-effect (mu = {mu2_1})")

ax.set_xlabel(r"$\nu$, GHz")
ax.set_ylabel(r"SZ signal $\mu$K")
ax.legend(frameon=False)

fig.savefig('pics/sz_effect.png', bbox_inches='tight')
plt.show()
