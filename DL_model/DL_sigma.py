#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
conductividad de Ag
"""
import numpy as np
graficar = 0

#%%

def sigma_DL(omega,omega_bulk,gamma_in): 
    """
    Parameters
    ----------
    omega : frecuencia en Hz
    omega_bulk : frecuencia de bulk Hz
    gamma_in : frecuencia de colision en eV
    d : espesor del plano Ag en micrometros
    Returns
        conductividad del Ag 
    -------

    """
    hb = 6.58211899*10**(-16)     ### Planck constant hbar in eV*s
    den = omega*(omega + 1j*gamma_in/hb)
    num = omega_bulk**2

    return (num/den)*1j*omega/(4*np.pi)

#%%


#%%
