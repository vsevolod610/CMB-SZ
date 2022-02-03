import numpy as np
import matplotlib.pyplot as plt

def xi2_fit(y, sigma):
    sigma2fit = (np.sum(1 / sigma ** 2)) ** -1
    yfit = sigma2fit * np.sum(y / sigma ** 2)
    return yfit,  np.sqrt(sigma2fit)

N = 77
z = np.abs(np.random.normal(0, 0.3, N))
T0 = 2.729
sigmaT0 = 0.14
T = np.random.normal(T0, sigmaT0, N)
sigmaT = np.random.normal(sigmaT0, 0.2 * sigmaT0,N)
T0fit, sigmaTfit = xi2_fit(T, sigmaT)


fig, ax = plt.subplots(figsize=(12, 8))
x = z
y = T * (1 + z)
sy = sigmaT * (1 + z)
ax.errorbar(x, y, sy, capsize=3.5, mew=1.5, fmt='.k')
ax.plot(x, T0fit * (1 + x), label=r'fit $T_0 = {0} \pm {1}$'.format(round(T0fit, 3), round(sigmaTfit, 3)))
ax.plot(x, T0 * (1 + x), label=r'teor $T_0 = {0} \pm {1}$'.format(T0, sigmaT0))
ax.set_xlabel(r'$z$')
ax.set_title(r'$T(z)$'.format(N))
ax.legend(frameon=False, loc='lower right')
plt.show()