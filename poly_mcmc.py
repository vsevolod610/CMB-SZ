# -*- coding: utf-8 -*-
"""  
mcmc realyze N clusters
"""

import numpy as np
from chainconsumer import ChainConsumer

from data import read_SZ_data, read_prior
from data import params_names, ndim, nsteps, nwalkers, init, prior_data
from mcmc import SZmcmc, pic_chain, pic_fit
from data_create import N, path_to_NSZ_data, path_to_Npriors


# Path to result.txt
path_result = './data/N/result.txt'
path_pic_chain = './data/N/pic_chain/pic_chain{}.png'
path_pic_consumer = './data/N/pic_consumer/pic_consumer{}.png'
path_pic_fit = './data/N/pic_fit/pic_fit{}.png'


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
    T0_res = []
    for i in range(N):
        # data
        x, y, yerr = read_SZ_data(path_to_NSZ_data.format(i + 1))
        prior_data['gauss'] = read_prior(path_to_Npriors.format(i + 1))

        # mcmc
        print(i + 1)
        sampler = SZmcmc(x, y, yerr, nwalkers, nsteps, ndim, init, prior_data)
        amputete = int(0.5 * nsteps)

        flat_sample = sampler.chain[:, amputete:, :].reshape((-1, ndim))
        c = ChainConsumer()
        c.add_chain(flat_sample, parameters=params_names)
        summary = c.analysis.get_summary(parameters=params_names)
        
        if not None in summary['T0']:
            T0m, T0, T0p = summary['T0']
            T0_res.append([T0, T0p - T0, T0 - T0m])
        else:
            T0_res.append(None)

        #pic's
        fig = pic_chain(sampler)
        fig.savefig(path_pic_chain.format(i + 1))

        fig = c.plotter.plot(display=False, legend=False, figsize=(6, 6))
        fig.savefig(path_pic_consumer.format(i + 1))

        fig = pic_fit(sampler, x, y, yerr, prior_data)
        fig.savefig(path_pic_fit.format(i + 1))


    # write the resutls
    print("MCMC:      T0      sigma T0   sigma T0")
    with open(path_result, 'w') as file:
        for k, t0 in enumerate(T0_res, start=1):
            if t0 is not None:
                s = "{:>2}     {:<10.6} {:<10.6} {:<10.6}".format(k, *t0)
            else:
                s = "{:>2}     {:>30}".format(k, 'None')

            file.write(s + "\n")
            print(s)
