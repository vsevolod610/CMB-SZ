""""
Read spectral transmission functions of Plank telescope for 70, 100, 143, 217, 353 GHz

require:    LFI_RIMO_R3.31.fits, HFI_RIMO_R3.00.fits

coments:    LFI: wavenubers in GHz, normalized, amputate is not need
            HFI: wavenubers in 1/cm, not normalized, amputate is need
"""

import numpy as np
from astropy.io import fits


filtrs = [70, 100, 143, 217, 353]
spec_trans = dict.fromkeys(filtrs)
spec_data = dict.fromkeys(filtrs)

with fits.open('LFI_RIMO_R3.31.fits', memmap=False) as hdulist:
    spec_data[70] = np.array(hdulist[27].data)

with fits.open('HFI_RIMO_R3.00.fits', memmap=False) as hdulist:
    spec_data[100] = np.array(hdulist[3].data)
    spec_data[143] = np.array(hdulist[4].data)
    spec_data[217] = np.array(hdulist[5].data)
    spec_data[353] = np.array(hdulist[6].data)

for wave in filtrs:
    nu = np.array([el[0] for el in spec_data[wave]])
    # 1 GHz = 30 * 1/cm
    nu = nu * 30 if (wave > 70) else nu
    tr = np.array([el[1] for el in spec_data[wave]])

    # sp. trans. functions have long tails of zero(1e-34) values, speed of calculus need amputation
    amputate = 1
    if amputate:
        tr = tr / max(tr)
        nu = nu[tr > 1e-3]
        tr = tr[tr > 1e-3]

    S = np.trapz(tr, nu)
    tr = tr / S

    spec_trans[wave] = {'nu': nu, 'tr': tr}
