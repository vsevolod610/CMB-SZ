""""
MCMC Analysis with ChainConsumer

require:    MCchain.dat, chainconsum.dat
make:       none

coments:
"""


import numpy as np
import matplotlib.pyplot as plt
from chainconsumer import ChainConsumer


parameters = [r'$T_0$', '$kT_e$', r'$\beta$', r'$\mathcal{T}$']
num_par = len(parameters)


# Read MCcain
MCchain = []
with open('../fortran/MCchain.dat') as f:
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

#Read cainconsum
chain = []
with open('../fortran/chainconsum.dat') as f:
    for k, line in enumerate(f):
        chain.append(line.split())
chain = np.array(chain, dtype=float)
# amputate chain
amputate = 0
if amputate:
    chain = chain[chain[: 0] > 400]
# Chainconsumer Analysis
c = ChainConsumer()
array = chain[:, 2:]
c.add_chain(array, parameters=parameters, name='sz')
print(c.analysis.get_latex_table(parameters=parameters))
fig = c.plotter.plot(display=False, legend=False, figsize=(6, 6))

plt.show()
