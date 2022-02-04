""""
Show data, SZfuntion, consider spectral transmission

require:   data.py, SZfunction.py
           LFI_RIMO_R3.31.fits, HFI_RIMO_R3.00.fits, trans_function.py
make:      none

coments:
"""


import numpy as np
import matplotlib.pyplot as plt

from SZfunction import SZfunction
from trans_function import spec_trans, filtrs
from data import nu_data, SZ_data, SZ_data_errors

T0 = 2.72
kTe = 6.9
beta = 0.0
Tau = 1.12

def SZ_tr(T0, kTe, beta, Tau):
    l = dict.fromkeys(filtrs)
    for wave in filtrs:
        nu = spec_trans[wave]['nu']
        tr = spec_trans[wave]['tr']
        l[wave] = np.trapz(SZfunction(T0, kTe, beta, Tau, nu) * tr, nu)
    return l

fig, ax = plt.subplots(figsize=(12, 8))
nu = np.linspace(nu_data[0], nu_data[-1], 100)
ax.plot(nu, 0 * nu, 'k')
ax.errorbar(nu_data, SZ_data, SZ_data_errors, capsize=3.5, mew=1.5, fmt='.k')
ax.plot(nu, SZfunction(T0, kTe, beta, Tau, nu))
SZ = np.array([SZ_tr(T0, kTe, beta, Tau)[wave] for wave in filtrs])
ax.plot(nu_data, SZ, 'ob')

plt.show()