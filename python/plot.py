#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""""
set of T0_i temperature --> Estimate T0

require:    ../Result/Result.txt
"""

import numpy as np
import matplotlib.pyplot as plt
import emcee
from multiprocessing import Pool
from chainconsumer import ChainConsumer


# data
res = []
with open('Result_1.txt') as f:
    for k, line in enumerate(f, start=1):
        if 'None' in line.split():
            continue
        res.append(line.split())
    res = np.array(res, dtype=float)
x1 = res[:, 1]
y1 = res[:, 2] * (1 + x1)
yerr1 = res[:, [3, 4]]
yerr1 = [err * (1 + x1[i]) for i, err in enumerate(yerr1)]

# data
res = []
with open('Result_2.txt') as f:
    for k, line in enumerate(f, start=1):
        if 'None' in line.split():
            continue
        res.append(line.split())
    res = np.array(res, dtype=float)
x2 = res[:, 1]
y2 = res[:, 2]
yerr2 = res[:, [3, 4]]

def model(T0, z):
    return T0 * (1 + z)


# data pic
x1 = list(x1)
x2 = list(x2)
x = list(set(x1) & set(x2))
y = [(y2[x2.index(i)] - y1[x1.index(i)]) * 2 / (yerr1[x1.index(i)][0] + yerr1[x1.index(i)][1]) for i in x]
xline = np.linspace(min(x1 + x2) - 0.2, max(x1 + x2) + 0.2, 100)

fig = plt.figure(figsize=(10, 8))
xlim = (min(x1 + x2) - 0.01, max(x1 + x2) + 0.01)
ax1 = fig.add_axes([0.1, 0.3, 0.8, 0.6], xticklabels=[], xlim=xlim)
ax2 = fig.add_axes([0.1, 0.1, 0.8, 0.2], ylim=(-3, 3), xlim=xlim)
ax1.errorbar(x1, y1, np.transpose(yerr1), capsize=3.5, mew=1.5, fmt='.k', alpha=0.3, label='1) T0 - ?')             
ax1.errorbar(x2, y2, np.transpose(yerr2), capsize=3.5, mew=1.5, fmt='.r', alpha=0.3, label='2) Tz - ?')             
ax2.plot(xline, xline * 0, 'k')
ax2.plot(x, y, '.b')
mean = np.mean(y) 
ax2.plot(xline, mean + xline * 0, '--b')
sigma = np.mean((y - mean) ** 2) ** 0.5
ax2.fill_between(xline, xline * 0 + mean + sigma, xline * 0 + mean - sigma , alpha=0.5, color='b', linewidth=1, linestyle='-')
ax2.fill_between(xline, xline * 0 + 1, xline * 0 - 1, alpha=0.3, color='g', linewidth=2, linestyle='-')
ax2.fill_between(xline, xline * 0 + 2, xline * 0 - 2, alpha=0.2, color='g', linewidth=2, linestyle='-')
ax2.set_xlabel(r'z', fontsize=12)
ax1.set_ylabel(r'$T_z$', fontsize=12)
ax1.legend()
plt.show()
