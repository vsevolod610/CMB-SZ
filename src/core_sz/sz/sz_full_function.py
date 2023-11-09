'''
Shimon, Rephaeli - 2004
'''

import numpy as np
from scipy import integrate

def sz_full_scalar(theta, beta_c, mu_c, x):
    '''
    mu = cos
    '''
    #print(f"full ... {x = }")

    # definitions
    gamma = lambda beta: (1 / np.sqrt(1 - beta ** 2))
    z     = lambda mu0, mu0h, beta: ((1 + beta * mu0h) / (1 + beta * mu0))
    mu    = lambda mu0, beta: ((mu0 + beta) / (1 + beta * mu0))

    f     = lambda mu0, mu0h: 3 / 8 * (
            1 + (mu0 * mu0h) ** 2 + 1 / 2 * (1 - mu0 ** 2) * (1 - mu0h ** 2))

    bbc   = lambda phi, mu0, beta: (
            beta_c * beta * (mu(mu0, beta) * mu_c + 
            np.cos(phi) * np.sqrt((1 - mu(mu0, beta) ** 2) * (1 - mu_c ** 2))))

    rho   = lambda phi, mu0, beta: (
            np.exp(gamma(beta_c) * gamma(beta) * 
            (bbc(phi, mu0, beta) - 1) / theta))

    delta = lambda mu0, mu0h, beta: (
            1 / (np.exp(z(mu0, mu0h, beta) * x) - 1) - 1 / (np.exp(x) - 1))


    # integrate
    sz_func = lambda phi, mu0, mu0h, beta: (
            gamma(beta_c) * (1 - bbc(phi, mu0, beta)) * f(mu0, mu0h) * 
            delta(mu0, mu0h, beta) / (gamma(beta) ** 4 * (1 + beta * mu0) ** 3) * 
            rho(phi, mu0, beta) *  gamma(beta) ** 5 * beta ** 2)

    sz, _ = integrate.nquad(
            sz_func, [[0, 2 * np.pi], [-1, 1], [-1, 1], [0, 1]], 
            opts={'epsrel': 1e-3, 'limit': 4})

    sz = sz / (2 * np.pi)

    # norm
    D_func  = lambda phi, mu0, beta: (
            gamma(beta_c) * (1 - bbc(phi, mu0, beta)) * 
            rho(phi, mu0, beta) * gamma(beta) ** 5 * beta ** 2)

    D, _ = integrate.nquad(
            D_func, [[0, 2 * np.pi], [0, 1], [0, 1]], 
            opts={'epsrel': 1e-3, 'limit': 4})

    D = D / (4 * np.pi)

    return sz / D * (1 * np.sinh(x / 2) ** 2 / x)

# NB! vectorize funciton
sz_full_function = np.vectorize(sz_full_scalar)

