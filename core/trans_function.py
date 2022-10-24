# -*- coding: utf-8 -*-
"""
Set: Spectral transmission functions
"""

import numpy as np
from astropy.io import fits


#Paths
path_to_LFI = '../data/LFI_RIMO_R3.31.fits'
path_to_HFI = '../data/HFI_RIMO_R3.00.fits'

filtrs = [70, 100, 143, 217, 353]
trans_data = dict.fromkeys(filtrs)


# Read fits files
#LFI: wavenubers in GHz, normalized, amputate is not need
#HFI: wavenubers in 1/cm, not normalized, amputate is need
with fits.open(path_to_LFI, memmap=False) as hdulist:
    trans_data[70] = np.array(hdulist[27].data)

with fits.open(path_to_HFI, memmap=False) as hdulist:
    trans_data[100] = np.array(hdulist[3].data)
    trans_data[143] = np.array(hdulist[4].data)
    trans_data[217] = np.array(hdulist[5].data)
    trans_data[353] = np.array(hdulist[6].data)


# Convert trans-function data to union type
for wave in filtrs:
    nu = np.array([el[0] for el in trans_data[wave]])
    # 1 GHz = 30 * 1/cm
    nu = nu * 30 if (wave > 70) else nu
    tr = np.array([el[1] for el in trans_data[wave]])

    # sp. trans. functions have long tails of <(1e-34) values,
    # speed of calculus need amputation
    amputate = 1
    if amputate:
        tr = tr / max(tr)
        nu = nu[tr > 1e-3]
        tr = tr[tr > 1e-3]
    S = np.trapz(tr, nu)
    tr = tr / S

    trans_data[wave] = {'nu': nu, 'tr': tr}


# Gauss approximation of trans-function data
def gauss(x, mu, sigma2):
    return np.exp(-(x - mu) ** 2 / (2 * sigma2)) / np.sqrt(2 * np.pi * sigma2)


trans_gauss = dict.fromkeys(filtrs)
for wave in filtrs:
    x = trans_data[wave]['nu']
    y = trans_data[wave]['tr']
    mu = np.trapz(y * x, x) 
    sigma2 = np.trapz(y * (x - mu) ** 2, x) 

    trans_gauss[wave] = {'mu': mu, 'sigma2': sigma2}

