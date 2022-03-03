# -*- coding: utf-8 -*-
""""
set clusters --> set of T0_i temperature

require:    mcmc_analyze.py
            ../OUT/MCchain{1-77}.dat, ../OUT/chainconsum{1-77}.dat
make:       ../Result/MCcain_analyze{1-77}.png, ../Result/chainconsum_analyze{1-77}.png
            ../Result/Result.txt

coments:
"""

import matplotlib.pyplot as plt
from mcmc_analyze import MCcain_analyze, chainconsum_analize

filename = '../fortran/MCchain.dat'
MCcain_analyze(filename)
filename = '../fortran/chainconsum.dat'
estim = chainconsum_analize(filename, amputate=500)
estim = estim.get_summary()
print(estim)

plt.show()