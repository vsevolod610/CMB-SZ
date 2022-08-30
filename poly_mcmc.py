# -*- coding: utf-8 -*-
"""  
mcmc realyze N clusters
"""

import gc
import numpy as np
import matplotlib.pyplot as plt
from chainconsumer import ChainConsumer

from data import read_SZ_data, read_prior
from data import params_names, ndim, nsteps, nwalkers, init, prior_data
from data_create import N, path_to_NSZ_data, path_to_Npriors
from model import SZmodel, gauss_model
from mcmc import mcmc_kern, pic_chain, pic_fit


# Path to result.txt
path_result = './data/N/result.txt'
path_pic_chain = './data/N/pic_chain/pic_chain{}.png'
path_pic_consumer = './data/N/pic_consumer/pic_consumer{}.png'
path_pic_fit = './data/N/pic_fit/pic_fit{}.png'


start = 1
stop = N

def read_result(path_result):
    res = []
    with open(path_result) as file:
        for k, line in enumerate(file, start=1):
            if 'None' in line.split():
                continue
            res.append(line.split())
        res = np.array(res, dtype=float)
    n = res[:, 0]
    T0 = res[:, 1]
    T0_err = res[:, [2, 3]]

    return n, T0, T0_err


if __name__ == "__main__":
    # N cluster analysis
    for i in range(start, stop + 1):
        # data
        x, y, yerr, z = read_SZ_data(path_to_NSZ_data.format(i))
        prior_data['gauss'] = read_prior(path_to_Npriors.format(i))

        # mcmc
        print(i)
        sampler = mcmc_kern(
                gauss_model, nwalkers, nsteps, ndim, init, x, y, yerr, prior_data)
        amputete = int(0.5 * nsteps)

        flat_sample = sampler.chain[:, amputete:, :].reshape((-1, ndim))
        c = ChainConsumer()
        c.add_chain(flat_sample, parameters=params_names)
        summary = c.analysis.get_summary(parameters=params_names)
        
        # write the result
        if not None in summary['T0']:
            T0m, T0, T0p = summary['T0']
            t0 = [T0, T0p - T0, T0 - T0m]
            s = "{:>4}     {:<10.6} {:<10.6} {:<10.6}".format(i, *t0)
        else:
            s = "{:>4}     {:>30}".format(i, 'None')

        with open(path_result, 'a') as file:
            file.write(s + "\n")

        # pic's
        fig = pic_chain(sampler)
        fig.savefig(path_pic_chain.format(i))

        fig = c.plotter.plot(display=False, legend=False, figsize=(6, 6))
        fig.savefig(path_pic_consumer.format(i))

        fig = pic_fit(sampler, SZmodel, x, y, yerr, prior_data)
        fig.savefig(path_pic_fit.format(i))

        # garved collector
        plt.clf()
        plt.close()
        plt.close(fig)
        plt.close('all')
        del x
        del y
        del yerr
        del sampler
        del flat_sample
        del c
        del fig
        gc.collect()
