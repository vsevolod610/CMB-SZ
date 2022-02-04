""""
Show spectral transmission functions

require:   LFI_RIMO_R3.31.fits, HFI_RIMO_R3.00.fits, trans_function.py
make:      none

coments:   spectral transmission functions can be normalized with 1 or max
"""


import matplotlib.pyplot as plt
from trans_function import spec_trans, filtrs

fig, ax = plt.subplots(figsize=(12, 8))
for wave in filtrs:
    nu = spec_trans[wave]['nu']
    tr = spec_trans[wave]['tr']

    # if norm_max=0 : spectral transmission functions is normalized with 1
    # if norm_max=1 : spectral transmission functions is normalized with max
    norm_max = 0
    if norm_max == 1:
        tr = tr / max(tr)

    ax.plot(nu, tr, label=r'${}$ GHz'.format(wave))
ax.set_xlabel(r'$\nu$ GHz')
ax.set_ylabel(r'$t_{\nu_0}(\nu)$')
ax.legend(frameon=False, loc='upper right')

plt.show()