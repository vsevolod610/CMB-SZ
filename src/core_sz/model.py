#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Define models to analyze # model(*params, *const, x)
"""

# Imports: built-in
import numpy as np
from scipy import integrate

# Imports: from project
from .sz_signal import sz_signal


### Models
# ------------- test only -----------------


def model_000(T0, Te, z=0): # very light test model
    beta = 0.00
    Tau = 1.0
    Tz = T0 * (1 + z)
    return sz_signal(T0, Te, beta, Tau, Tz, corrs=False, mode='off')


# ------------- default -------------------


def model_100(T0, Te, beta, Tau, z=0, x=None): # light test model
    Tz = T0 * (1 + z)
    return sz_signal(T0, Te, beta, Tau, T0, corrs=False, mode='off')


def model_110(T0, Te, beta, Tau, z=0, x=None): # test model
    Tz = T0 * (1 + z)
    return sz_signal(T0, Te, beta, Tau, Tz, corrs=True, mode='off')


def model_111(T0, Te, beta, Tau, z=0, x=None): # heavy model
    Tz = T0 * (1 + z)
    return sz_signal(T0, Te, beta, Tau, Tz, corrs=True, mode='row')


def model_112(T0, Te, beta, Tau, z=0, x=None): # default model
    Tz = T0 * (1 + z)
    return sz_signal(T0, Te, beta, Tau, Tz, corrs=True, mode='gauss')


# ----------- alternative -----------------


def model_212(Tz, Te, beta, Tau, z, x=None): # default alternative model
    T0 = 2.7255
    return sz_signal(T0, Te, beta, Tau, Tz, z=z, corrs=True, mode='gauss')


def model_210(Tz, Te, beta, Tau, z, x=None): # test alternative model
    T0 = 2.7255
    return sz_signal(T0, Te, beta, Tau, Tz, z=z, corrs=True, mode='off')


