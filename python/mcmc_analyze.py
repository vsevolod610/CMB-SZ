# -*- coding: utf-8 -*-
""""
set MCMC Analyze cluster
"""


import numpy as np
import matplotlib.pyplot as plt
from chainconsumer import ChainConsumer


parameters = [r'$T_0$', '$kT_e$', r'$\beta$', r'$\mathcal{T}$']
num_par = len(parameters)

def MCcain_analyze(filename, picname=False):
    # Read MCcain
    MCchain = []
    with open(filename) as f:
        for k, line in enumerate(f, start=1):
            MCchain.append(line.split())
    MCchain = np.array(MCchain, dtype=float)
    # Plot MCchain
    fig, ax = plt.subplots(nrows=num_par, figsize=(8, 10))
    for i, row in enumerate(ax, start=0):
            x = MCchain[:, 0]
            y = MCchain[:, 1 + i]
            y_errors = np.transpose(MCchain[:, num_par + i + 1])
            row.errorbar(x, y, yerr=y_errors)
            row.set_ylabel(parameters[i], fontsize=12)
    if picname != False:
        fig.savefig(str(picname))



def chainconsum_analize(filename, amputate=0 , picname=False):
    #Read cainconsum
    chain = []
    with open(filename) as f:
        for k, line in enumerate(f):
            chain.append(line.split())
    chain = np.array(chain, dtype=float)
    # amputate chain
    if amputate:
        chain = chain[chain[:, 0] > int(amputate)]
    # Chainconsumer Analysis
    c = ChainConsumer()
    array = chain[:, 2:] # chain : num chi2 paramns...
    c.add_chain(array, parameters=parameters, name='sz')
    # Plat CainConsumer
    fig = c.plotter.plot(display=False, legend=False, figsize=(6, 6))
    if picname != False:
        fig.savefig(str(picname))
    return c.analysis
