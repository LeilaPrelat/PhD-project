#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila

green tensor 
"""
from scipy import special
import numpy as np
import sys
import os 

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
# path_graphene =  path_basic.replace('/' + 'plane','')
#print('Importar modulos necesarios para este codigo')

try:
    sys.path.insert(1, path_basic)
    from constants import constantes
except ModuleNotFoundError:
    print('constants.py no se encuentra en ' + path_basic)

pi,hb,c,alfac,mu1,mu2 = constantes()
aux = c*hb

#%%

def green_tensor_functionDIR_2(omegac,epsi1,x,y,z,xD,yD,zD):
    """
    Parameters
    ----------
    omegac : omega/c = k0 en 1/micrometros    
    epsi1 : epsilon del medio de arriba del plano
    x : coordenada x 
    y : coordenada y 
    z : coordenada z
    xD : coordenada x del dipolo 
    yD : coordenada y del dipolo 
    zD : coordenada z del dipolo 
    Returns
    -------
    Gxx direct (self interaction of the dipole)
    con z hacia abajo (convencion del paper)
    """
    
    k0 = omegac #=omega/c
    # x_y = ky/k0
    n1 = epsi1*mu1
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1
    k1_3 = (k0*cte1)**3
    Rbarra =  np.sqrt((x-xD)**2 + (y-yD)**2)/k1
    z_dip_barra = np.abs(z-zD)/k1
    phi = np.arctan2(np.abs(y-yD),np.abs(x-xD))
    cte = 0.5*1j*k1_3
    
    arg1 = np.sqrt(z_dip_barra + 1j*Rbarra)
    arg2 = np.sqrt(z_dip_barra - 1j*Rbarra)
    aux = (z_dip_barra**2 + Rbarra**2)**(-3/2)
    aux0_5 = 0.25*k1_3/np.sqrt(Rbarra)
    
    I0_5 = aux0_5*(-(arg1 + arg2) + 2*np.sin(np.arctan(0.5*Rbarra/z_dip_barra))*(aux**(1/4)))
    I0_6 = -cte*z_dip_barra*aux

    aux2_5 = k1_3/(Rbarra**2)
    
    aux_term = (np.sqrt(z_dip_barra**2+Rbarra**2) - z_dip_barra)**2
    I2_5 = 0.25*aux2_5*np.cos(2*phi)*aux_term
    I2_6 = -0.5*aux2_5*np.cos(2*phi)*aux*aux_term*(2*np.sqrt(z_dip_barra**2+Rbarra**2) + z_dip_barra)
    
    return I0_5, I2_5, I0_6, I2_6 

#%%
