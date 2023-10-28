#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Show: models
"""

# Imports
import datetime
import numpy as np
import matplotlib.pyplot as plt 
from scipy import integrate

import paths # manage path imports by paths.py
from core_sz.trans_function import filtrs
from core_sz.model import model_112, model_100, model_110, model_111


# example params: T0, Te, beta, Tau
exampe_params = [2.7255, 6.9, 0.0, 1.4]


# models analisys
sz_100 = model_100(*exampe_params) # without corrs
sz_110 = model_110(*exampe_params) # clear
sz_111 = model_111(*exampe_params) # row trans
sz_112 = model_112(*exampe_params) # gauss trans
print(f"Waves = {filtrs} GHz")
print("100: SZ clear (no corrs): ", *list(map(" {:5.0f},".format, sz_100)))
print("110: SZ clear:            ", *list(map(" {:5.0f},".format, sz_110)))
print("111: SZ row trans-func:   ", *list(map(" {:5.0f},".format, sz_111)))
print("112: SZ gauss trans-func: ", *list(map(" {:5.0f},".format, sz_112)))

# Picture:
fig, ax = plt.subplots(figsize=(8, 8))

ax.axhline(y=0, color='k')
ax.plot(filtrs, sz_100, 'o-', label="SZ clear (no corrs)")
ax.plot(filtrs, sz_110, 'o-', label="SZ clear")
ax.plot(filtrs, sz_111, 'o-', label="SZ row trans-func")
ax.plot(filtrs, sz_112, 'o-', label="SZ gauss trans-func")

ax.set_xlabel(r"$\nu$, GHz")
ax.set_ylabel(r"SZ signal $\mu$K")
ax.legend(frameon=False)

fig.savefig('pics/models.png', bbox_inches='tight')
plt.show()
