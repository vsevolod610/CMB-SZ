#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Intergrate sz-effect + spectral-transmission function
"""

# Imports: from project
from .trans.trans_intergrate import trans_integrate
from .sz.sz_function import sz_function
from .sz.sz_full_function import sz_full_function


# Constants
h = 2 * 3.141 * 1.05e-27
c = 3e10
k = 1.38e-16
m = 0.91e-27
coef1 = (h * 1e9) / k
coef2 = (1e3 * 1.6e-12) / (m * c ** 2)


# FIXME: z = None, z = 0
def sz_effect(nu, T0, Te, beta, Tau, Tz, z=None, corrs=True):
    # if z is None Tz not needed
    theta = coef2 * Te
    x = (coef1 * nu * (1 + z) / Tz) if (z is not None) else (coef1 * nu / T0)
    return 1e+6 * T0 * Tau * sz_function(theta, beta, x, corrs)


def sz_full_effect(nu, T0, Te, beta, Tau, Tz, z=None, mu=1):
    # if z is None Tz not needed
    theta = coef2 * Te
    x = (coef1 * nu * (1 + z) / Tz) if (z is not None) else (coef1 * nu / T0)
    #print(f"{mu = }")
    return 1e+6 * T0 * Tau * sz_full_function(theta, beta, mu, x)


# sz-effect + spec. transmission
def sz_signal(T0, Te, beta, Tau, Tz, z=0, corrs=True, mode='guass'):
    args = (T0, Te, beta, Tau, Tz, z, corrs)
    sz_signal = trans_integrate(sz_effect, args=args, mode=mode)
    return sz_signal
