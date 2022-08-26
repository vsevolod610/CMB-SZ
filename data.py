# -*-import matplotlib.pyplot as plt coding: utf-8 -*-
"""
Read Data
"""

import numpy as np
import matplotlib.pyplot as plt

def read_SZ_data(path):
    # DATA: import data from SZ_data.txt
    data = []
    with open(path) as file:
        for k, line in enumerate(file):
            if k == 1: z = float(line) # redshift
            if k < 3: continue # skip first 3 lines
            data.append(line.split())
        data = np.array(data, dtype=float)

    x = np.array([70.0, 100.0, 143.0, 217.0, 353.0])
    y = data[:, 0]
    yerr = data[:, [1, 2]]

    return x, y, yerr


def read_prior(path):
    # DATA: impot prior from prior.dat
    data = dict()
    with open(path) as file:
        for k, line in enumerate(file):
            if k < 1: continue # skip first line
            param = line.split()[2]
            data[param] = np.array(line.split()[3:6], dtype=float)

    return data


def read_startSZ(path):
    # DATA: import start position from startSZ.sss
    with open(path) as file:
        text = file.readlines()
        nsteps = int(text[3].split()[0])
        nwalkers = int(text[4].split()[0])
        data = []
        for k, line in enumerate(text):
            if k < 6: continue # skip first 3 line
            data.append(line.split()[1:5])
        data = np.array(data, dtype=float)

    return nsteps, nwalkers, data 

# Paths to files
path_to_SZ_data = "./data/SZ_data.txt"
path_to_prior = "./data/prior.dat"
path_to_startSZ = "./data/startSZ.sss"

# Names
params_names = [r'T0', r'Te', r'beta', r'Tau']
params_names_latex = [r'$T_0$', r'$T_e$', r'$\beta$', r'$\mathfrac{T}']
ndim = len(params_names)

# Set data
x, y, yerr = read_SZ_data(path_to_SZ_data)
prior_data = dict()
prior_data['gauss'] = read_prior(path_to_prior)
nsteps, nwalkers, data = read_startSZ(path_to_startSZ)
prior_data['box'] = data[:, [1,2]]
init = data[:, [0, 3]]

if __name__ == "__main__":
    # SZ data
    print("\nSZ data:")
    N = len(x)
    s = " {:>5}: {:>8}  {:>+8} {:>8}"
    print(*[s.format(x[i], y[i], *[1,-1]*yerr[i]) for i in range(N)], sep = '\n')

    # Prior
    print("\nPriors")
    print("  box:")
    s = "    {:<4}: [{:^6}, {:>5}]"
    l = [s.format(params_names[i], *prior_data['box'][i]) for i in range(ndim)]
    print(*l, sep='\n')

    print("  gauss:")
    names = list(prior_data['gauss'].keys())
    s = "    {:<4}: {} {:>+3} -{:>3}"
    l = [s.format(name, *prior_data['gauss'][name]) for name in names]
    print(*l, sep='\n')

    # Start
    print("\nStart:")
    s = "  {:<4}: {:>4} +- {:<5}"
    l = [s.format(params_names[i], *init[i]) for i in range(ndim)]
    print(*l, sep='\n')

    # Pic
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.errorbar(x, y, yerr.T, label='data',
            capsize=3.5, mew=1.5, fmt='.k', alpha=0.5)
    ax.plot(x, 0 * y, 'k')

    ax.set_xlabel(r"$\nu$, GHz")
    ax.set_ylabel(r"SZ signal $\mu$K")
    ax.legend(frameon=False)

    plt.show()

