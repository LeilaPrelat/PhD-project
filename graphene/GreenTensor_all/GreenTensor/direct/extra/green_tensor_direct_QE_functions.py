#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila

funciones dentro del green tensor directo (antes de integrar)
"""
from scipy import special
import numpy as np
import sys
import os 

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('GreenTensor/' + 'direct/extra','')
#print('Importar modulos necesarios para este codigo')

try:
    sys.path.insert(1, path_constants)
    from constants import constantes
except ModuleNotFoundError:
    print('constants.py no se encuentra en ' + path_constants)

pi,hb,c,alfac,mu1,mu2 = constantes()
aux = c*hb

#%%


def functions_QE_GreenTensor(omegac,epsi1,alpha_parallel,x,y,z,xD,yD,zD):
    """
    Parameters
    ----------
    omegac : omega/c = k0 en 1/micrometros    
    epsi1 : epsilon del medio de arriba del plano
    alpha_parallel : alpha_parallel (k_parallel/k1)
    x : coordenada x 
    y : coordenada y 
    z : coordenada z
    xD : coordenada x del dipolo 
    yD : coordenada y del dipolo 
    zD : coordenada z del dipolo 
    Returns
    -------
    Gxx direct antes de integrar 
    con z hacia abajo (convencion del paper)
    con QE
    """

    k0 = omegac #=omega/c
    # x_y = ky/k0
    n1 = epsi1*mu1
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1
    k1_3 = (k0*cte1)**3
    Rbarra =  k1*np.sqrt((x-xD)**2 + (y-yD)**2)
    z_dip_barra = k1*np.abs(z-zD)
    phi = np.arctan2(np.abs(y-yD),np.abs(x-xD))

    exp = np.exp(-alpha_parallel*z_dip_barra)
    
    J0 = special.jn(0,alpha_parallel*Rbarra) 
    J2 = special.jn(2,alpha_parallel*Rbarra) 

############ todo junto ###################################
    int_re = np.real((J0 + np.cos(2*phi)*J2)*exp*alpha_parallel*(alpha_parallel - 1/alpha_parallel))
    int_im = np.imag((J0 + np.cos(2*phi)*J2)*exp*alpha_parallel*(alpha_parallel - 1/alpha_parallel))

    rta = int_re + 1j*int_im
    
    cte = -0.5*k1_3 
    
    return rta*cte

#%%
