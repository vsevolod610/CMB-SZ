# -*- coding: utf-8 -*-
"""
Spectral transmission functions
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits


#Paths
path_to_LFI = './data/LFI_RIMO_R3.31.fits'
path_to_HFI = './data/HFI_RIMO_R3.00.fits'

filtrs = [70, 100, 143, 217, 353]
spec_trans = dict.fromkeys(filtrs)
spec_data = dict.fromkeys(filtrs)


# Read fits files
#LFI: wavenubers in GHz, normalized, amputate is not need
#HFI: wavenubers in 1/cm, not normalized, amputate is need
with fits.open(path_to_LFI, memmap=False) as hdulist:
    spec_data[70] = np.array(hdulist[27].data)

with fits.open(path_to_HFI, memmap=False) as hdulist:
    spec_data[100] = np.array(hdulist[3].data)
    spec_data[143] = np.array(hdulist[4].data)
    spec_data[217] = np.array(hdulist[5].data)
    spec_data[353] = np.array(hdulist[6].data)


# Conver to union type
for wave in filtrs:
    nu = np.array([el[0] for el in spec_data[wave]])
    # 1 GHz = 30 * 1/cm
    nu = nu * 30 if (wave > 70) else nu
    tr = np.array([el[1] for el in spec_data[wave]])

    # sp. trans. functions have long tails of <(1e-34) values,
    # speed of calculus need amputation
    amputate = 1
    if amputate:
        tr = tr / max(tr)
        nu = nu[tr > 1e-3]
        tr = tr[tr > 1e-3]
    S = np.trapz(tr, nu)
    tr = tr / S

    spec_trans[wave] = {'nu': nu, 'tr': tr}


# Gauss approximation of trans-function
def gauss(x, mu, sigma2):
    return np.exp(-(x - mu) ** 2 / (2 * sigma2)) / np.sqrt(2 * np.pi * sigma2)


gauss_trans = dict.fromkeys(filtrs)
for wave in filtrs:
    x = spec_trans[wave]['nu']
    y = spec_trans[wave]['tr']
    mu = np.trapz(y * x, x) 
    sigma2 = np.trapz(y * (x - mu) ** 2, x) 

    gauss_trans[wave] = {'mu': mu, 'sigma2': sigma2}


if __name__ == "__main__":
    # trans function analysis
    print("\n trans funtions: nu line")
    for wave in filtrs:
        nu = spec_trans[wave]['nu']
        print(
            "Wave {:>3} GHz: nu = [{:>6.2f} ... {:>6.2f}] - {}"
            .format(wave, min(nu), max(nu), len(nu)))

    # gauss approximation analysis
    print("\n gauss approximation")
    for wave in filtrs:
        mu = gauss_trans[wave]['mu']
        sigma = np.sqrt(gauss_trans[wave]['sigma2'])
        print(
            "Wave {:>3} GHz: mu = {:>6.2f}, sigma = {:>5.2f}"
            .format(wave, mu, sigma))

    # pictures and gauss approximation
    fig, ax = plt.subplots(figsize=(8, 8))
    for wave in filtrs:
        x = spec_trans[wave]['nu']
        y = spec_trans[wave]['tr']
        mu = gauss_trans[wave]['mu']
        sigma2 = gauss_trans[wave]['sigma2']

        ax.plot(x, y, label='{} GHz'.format(wave))
        ax.plot(x, gauss(x, mu, sigma2), 'k')

        ax.set_xlabel(r'$\nu$ GHz')
        ax.legend(frameon=False)

    plt.show()
