#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila

green tensor : integrales resueltas analiticamente
luego de aplicar la aprox QE 
"""
import numpy as np
import sys
import os 

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('GreenTensor/' + 'direct','')
#print('Importar modulos necesarios para este codigo')

try:
    sys.path.insert(1, path_constants)
    from constants import constantes
except ModuleNotFoundError:
    print('constants.py no se encuentra en ' + path_constants)

pi,hb,c,alfac,mu1,mu2 = constantes()
aux = c*hb

#%%

def green_tensor_ANA_QE(omegac,epsi1,x,y,z,xD,yD,zD):
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
    analitico luego de aplicar QE
    """
    
    k0 = omegac #=omega/c
    # x_y = ky/k0
    n1 = epsi1*mu1
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1
    k1_3 = (k0*cte1)**3
    
    Rbarra = k1*np.sqrt((x-xD)**2 + (y-yD)**2)  #adimensional
    Rbarra_2 = Rbarra**2

    z_dip_barra = k1*np.abs(z-zD) #adimensional
    z_dip_barra_2 = z_dip_barra**2
    
    phi = np.arctan2(np.abs(y-yD),np.abs(x-xD))

    aux0 = z_dip_barra**2 + Rbarra**2
    aux1 = aux0**(-1/2)
    aux2 = aux0**(1/2) - z_dip_barra
    
    I0_5 = 0.5*k1_3*aux1
    I2_5 = 0.5*(k1_3/Rbarra_2)*np.cos(2*phi)*aux1*(aux2**2) 
    I0_6 = 0.5*k1_3*(3*z_dip_barra_2/(aux0**(5/2)) - aux0**(-3/2))

    # derivada del I6 orden 2 ### overleaf
    term1 = 3*z_dip_barra_2*(aux2**2)/(aux0**(5/2))
    term2 = 2*z_dip_barra*aux2*(2*z_dip_barra*aux1-2)/(aux0**(3/2))
    term3 = aux2**2/(aux0**(3/2))
    term4 = aux2*(-2*z_dip_barra_2*(aux1**3) + 2*aux1)*aux1
    term5 = 2*((z_dip_barra*aux1-1)**2)*aux1
    
    I2_6 = 0.5*(k1_3/Rbarra_2)*np.cos(2*phi)*(term1 - term2 - term3 + term4 + term5)
    
    return I0_5, I2_5, -I0_6, -I2_6 

#%%
