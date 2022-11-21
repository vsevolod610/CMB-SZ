# -*- coding: utf-8 -*-
"""
set of T0_i --> Estimate T0
comment: нодо бы перенести этот файл в корень
"""

import numpy as np
import matplotlib.pyplot as plt
from chainconsumer import ChainConsumer

from mcmc_kern import mcmc_kern, pic_chain, pic_fit
from start_analyze import result_read, path_result, method


x, y, yerr = result_read(path_result)


def modelT0(T0, z):
    return T0 + (0 * z)


def modelTz(T0, z):
    return T0 * (1 + z)


if method == 'T0':
    model = modelT0
if method == 'Tz':
    model = modelTz

init = np.array([[2.7255], [0.2]]).T
nwalkers = 100
nsteps = 100
amputete = int(0.2 * nsteps)
ndim = len(init)

# mcmc for estimation
sampler = mcmc_kern(model, nwalkers, nsteps, init, x, y, yerr)
flat_sample = sampler.chain[:, amputete : , :].reshape((-1, ndim))
c = ChainConsumer()
c.add_chain(flat_sample)
summary = c.analysis.get_summary()
print("\nMCMC results:")
print(*[" {:>4}: {}".format(k, summary[k]) for k in summary.keys()], sep='\n')

# Pics
fig, ax = pic_chain(sampler)
fig = c.plotter.plot(display=False, legend=False, figsize=(6, 6))
fig, ax = pic_fit(sampler, model, x, y, yerr)

plt.show()



