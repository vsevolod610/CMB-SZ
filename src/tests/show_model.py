#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Show: models
"""

# Imports: built-in
import datetime
import numpy as np
import matplotlib.pyplot as plt 

import paths # manage path imports by paths.py
from core_sz.trans.trans_data import filtrs
from core_sz.model import model_112, model_100, model_110, model_111


# example params: T0, Te, beta, Tau
exampe_params = [2.7255, 6.9, 0.0, 1.4]


def timer(func, args, N):
    start = datetime.datetime.now()
    for i in range(N):
        res = func(*args)
        del res
    finish = datetime.datetime.now()
    return (finish - start).total_seconds() / N


def see(model, params, comment):
    # models analisys
    sz = model(*params)
    sz_prnt = list(map(' {:5.0f},'.format, sz))

    # timer
    N = 100 # times in timers
    t = timer(model, exampe_params, N)

    # do
    print(f"{comment:<20}", *sz_prnt, f" <- {t:.6f} sec")
    ax.plot(filtrs, sz, 'o-', label=comment)


# Picture:
print(f"Waves = {filtrs} GHz")
fig, ax = plt.subplots(figsize=(8, 8))

ax.axhline(y=0, color='k')
see(model_100, exampe_params, comment="100: SZ (no corrs)")
see(model_110, exampe_params, comment="110: SZ           ")
see(model_111, exampe_params, comment="111: SZ + row     ")
see(model_112, exampe_params, comment="112: SZ + gauss   ")

ax.set_xlabel(r"$\nu$, GHz")
ax.set_ylabel(r"SZ signal $\mu$K")
ax.legend(frameon=False)

fig.savefig('pics/models.png', bbox_inches='tight')
plt.show()
