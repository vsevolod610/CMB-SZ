# -*- coding: utf-8 -*-
""""
mcmc analyze cluster

require:    mcmc_analyze.py
            ../fortran/MCchai.dat, ../fortran/chainconsum.dat
"""

import matplotlib.pyplot as plt
from mcmc_analyze import MCcain_analyze, chainconsum_analize


filename = '../fortran/MCchain.dat'
MCcain_analyze(filename)
filename = '../fortran/chainconsum.dat'
estim = chainconsum_analize(filename, amputate=300)
estim = estim.get_summary()
print(estim)

plt.show()
