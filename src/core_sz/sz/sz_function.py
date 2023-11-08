#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Set function of SZ-effect
"""

# Imports: built-in
import numpy as np

# Imports: from project
from .sz_terms import *


# x can be float or array
def sz_function(theta, beta, x, corrs=True):
    X, S = XS(x)

    # Withou relativistic corrections Y0
    Y0 = X - 4
    if corrs is False: return theta * Y0

    # Relativistic correctionss Y1, Y2, Y3, Y4, C1, C2, D0, D1
    Y1 = term_Y1(X, S)
    Y2 = term_Y2(X, S)
    Y3 = term_Y3(X, S)
    Y4 = term_Y4(X, S)
    C1 = term_C1(X, S)
    C2 = term_C2(X, S)
    D0 = term_D0(X, S)
    D1 = term_D1(X, S)
    
    R = (theta ** 2 * Y1 + theta ** 3 * Y2 + theta ** 4 * Y3 + theta ** 5 * Y4 
         + beta ** 2 * (1 / 3 * Y0 + theta * (5 / 6 * Y0 + 2 / 3 * Y1)) 
         + beta * (1 + theta * C1 + theta ** 2 * C2) 
         + beta ** 2 * (D0 + theta * D1))

    return theta * Y0 + R

