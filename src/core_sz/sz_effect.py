#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Set fuction of SZ-effect
"""

import numpy as np


# Constants
h = 2 * 3.141 * 1.05e-27
c = 3e10
k = 1.38e-16
m = 0.91e-27
coef1 = (h * 1e9) / k
coef2 = (1e3 * 1.6e-12) / (m * c ** 2)


def sz_function(theta, beta, x, corrs=True):
    X = x * (np.exp(x) + 1) / (np.exp(x) - 1)
    S = x / np.sinh(x / 2)

    # Withou relativistic corrections Y0
    Y0 = X - 4
    if corrs is False:
        return theta * Y0

    # Relativistic correctionss Y1, Y2, Y3, Y4, C1, C2, D0, D1
    Y1 = (- 10 + 47 / 2 * X - 42 / 5 * X ** 2 + 7 / 10 * X ** 3 
          + S ** 2 * (- 21 / 5 + 7 / 5 * X))

    Y2 = (- 15 / 2 + 1023 / 8 * X - 868 / 5 * X ** 2 + 329 / 5 * X ** 3 
          - 44 / 5 * X ** 4 + 11 / 30 * X ** 5 
          + S ** 2 * (- 434 / 5 + 658 / 5 * X - 242 / 5 * X ** 2 
                      + 143 / 30 * X ** 3) 
          + S ** 4 * (- 44 / 5 + 187 / 60 * X))

    Y3 = (15 / 2 + 2505 / 8 * X - 7098 / 5 * X ** 2 + 14253 / 10 * X ** 3 
          - 18594 / 35 * X ** 4 + 12059 / 140 * X ** 5 - 128 / 21 * X ** 6 
          + 16 / 105 * X ** 7 
          + S ** 2 * (- 7098 / 10 + 14253 / 5 * X - 102267 / 35 * X ** 2 
                      + 156767 / 140 * X ** 3 - 1216 / 7 * X ** 4 
                      + 64 / 7 * X ** 5) 
          + S ** 4 * (- 18594 / 35 + 205003 / 280 * X - 1920 / 7 * X ** 2 
                      + 1024 / 35 * X ** 3) 
            + S ** 6 * (- 544 / 21 + 992 / 105 * X))

    Y4 = (- 135 / 32 + 30375 / 128 * X - 62391 / 10 * X ** 2 
          + 614727 / 40 * X ** 3 - 124389 / 10 * X ** 4 
          + 355703 / 80 * X ** 5 - 16568 / 21 * X ** 6 + 7516 / 105 * X ** 7 
          - 22 / 7 * X ** 8 + 11 / 210 * X ** 9 
          + S ** 2 * (- 62391 / 20 + 614727 / 20 * X - 1368279 / 20 * X ** 2 
                      + 4624139 / 80 * X ** 3 - 157396 / 7 * X ** 4 
                      + 30064 / 7 * X ** 5 - 2717 / 7 * X ** 6 
                      + 2761 / 210 * X ** 7) 
          + S ** 4 * (- 124389 / 10 + 6046951 / 160 * X - 248520 / 7 * X ** 2 
                      + 481024 / 35 * X ** 3 - 15972 / 7 * X ** 4 
                      + 18689 / 140 * X ** 5) 
          + S ** 6 * (- 70414 / 21 + 465992 / 105 * X - 11792 / 7 * X ** 2 
                      + 19778 / 105 * X ** 3) 
          + S ** 8 * (- 682 / 7 + 7601 / 210 * X))

    C1 = (10 - 47 / 5 * X + 7 / 5 * X ** 2 + 7 / 10 * S ** 2)

    C2 = (25 - 1117 / 10 * X + 847 / 10 * X ** 2 - 183 / 10 * X ** 3 
          + 11 / 10 * X ** 4 
          + S ** 2 * (847 / 20 - 183 / 5 * X + 121 / 20 * X ** 2) 
          + 11 / 10 * S ** 4)

    D0 = (- 2 / 3 + 11 / 30 * X)

    D1 = (- 4 + 12 * X - 6 * X ** 2 + 19 / 30 * X ** 3 
          + S ** 2 * (- 3 + 19 / 15 * X))

    R = (theta ** 2 * Y1 + theta ** 3 * Y2 + theta ** 4 * Y3 + theta ** 5 * Y4
         + beta ** 2 * (1 / 3 * Y0 + theta * (5 / 6 * Y0 + 2 / 3 * Y1)) 
         + beta * (1 + theta * C1 + theta ** 2 * C2) 
         + beta ** 2 * (D0 + theta * D1))

    return theta * Y0 + R


# FIXME: if z is None Tz not needed
def sz_effect(T0, Te, beta, Tau, Tz, nu, z=None, corrs=True):
    theta = coef2 * Te

    if z is None:
        x = coef1 * nu / T0 # if z is not passed: Tz = T0 * (1 + z)
    else:
        x = coef1 * nu * (1 + z) / Tz

    return 1e+6 * T0 * Tau * sz_function(theta, beta, x, corrs)

