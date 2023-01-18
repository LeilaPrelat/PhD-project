#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
constantes pi,hb,c,alfac,mu1,mu2 que se mantienen 
a lo largo del proyecto
"""
import numpy as np

#%%

def constantes():
    """
    Returns
    -------
    pi : cte universal pi
    hb : cte universal hbar en eV*s
    c : cte universal vel de la luz en micron/seg
    alfac : cte universal de estructura fina
    hbargama : freq de colision del grafeno en eV
    mu1 : permeabilidad magnetica del medio 1 (medio arriba del plano)
    mu2 : permeabilidad magnetica del medio 2 (medio abajo del plano) 
    epsi1 : permeabilidad electrica del medio 1 (medio arriba del plano) 
    epsi2 : permeabilidad electrica del medio 2 (medio abajo del plano) 
    """

    global mu1,mu2,c ### variables globales ---> no pueden cambiarse
    pi = np.pi
    hb = 6.58211899*10**(-16)     ### Planck constant hbar in eV*s
    c = 3*10**(14)                ### light velocity in micron/seg
    alfac = 1/137.0359            ### fine structure
#    hbargama = 0.0001             ### collision frequency in eV para el grafeno
    
    mu1, mu2 = 1, 1
    
    return pi,hb,c,alfac,mu1,mu2

#%%
