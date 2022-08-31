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
from model import SZmodel, gauss_model, SZalternative, gauss_alternative
from mcmc import mcmc_kern, pic_chain, pic_fit


# Path to result.txt
path_result = './data/N/result{}.txt'
path_pic_chain = './data/N/pic_chain/pic_chain{}{}.png'
path_pic_consumer = './data/N/pic_consumer/pic_consumer{}{}.png'
path_pic_fit = './data/N/pic_fit/pic_fit{}{}.png'


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
    z = res[:, 1]
    T0 = res[:, 2]
    T0_err = res[:, [3, 4]]

    return n, z, T0, T0_err


def poly_mcmc_run(method, pic=False):
    # data
    x, y, yerr, z = read_SZ_data(path_to_NSZ_data.format(i))
    prior_data['gauss'] = read_prior(path_to_Npriors.format(i))

    if method == 'T0':
        params_names[0] = r'T0'
        model_function = SZmodel
        gauss_funciton = gauss_model
    elif method == 'Tz':
        params_names[0] = r'Tz'
        model_function = SZalternative
        gauss_funciton = gauss_alternative
        prior_data['const'] = {'0': z}
    else:
        print("error: method must be 'T0' or 'Tz'")

    # mcmc
    print(method, i)
    sampler = mcmc_kern(
            gauss_funciton, nwalkers, nsteps, ndim, init, x, y, yerr, prior_data)
    amputete = int(0.5 * nsteps)

    flat_sample = sampler.chain[:, amputete:, :].reshape((-1, ndim))
    c = ChainConsumer()
    c.add_chain(flat_sample, parameters=params_names)
    summary = c.analysis.get_summary(parameters=params_names)

    # write the result
    if not None in summary[method]:
        Tm, T, Tp = summary[method]
        t = [T, Tp - T, T - Tm]
        s = "{:>4}    {:<8}   {:<10.6} {:<10.6} {:<10.6}".format(i, z, *t)
    else:
        s = "{:>4}    {:<8}   {:>30}".format(i, z, 'None')

    with open(path_result.format(method), 'a') as file:
        file.write(s + "\n")

    if pic:
        # pic's
        fig = pic_chain(sampler)
        fig.savefig(path_pic_chain.format(method, i))

        fig = c.plotter.plot(display=False, legend=False, figsize=(6, 6))
        fig.savefig(path_pic_consumer.format(method, i))

        fig = pic_fit(sampler, SZfunction, x, y, yerr, prior_data)
        fig.savefig(path_pic_fit.format(method, i))

        plt.clf()
        plt.close()
        plt.close(fig)
        plt.close('all')
        del fig

    # garvege collector
    del x
    del y
    del yerr
    del sampler
    del flat_sample
    del c
    gc.collect()


if __name__ == "__main__":
    # N cluster analysis
    if 1:
        for i in range(start, stop + 1):
            poly_mcmc_run('T0')

    if 1:
        for i in range(start, stop + 1):
            poly_mcmc_run('Tz')

    
