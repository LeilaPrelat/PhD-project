#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila

funciones dentro del green tensor reflejado (antes de integrar)
"""
from scipy import special
import numpy as np
import sys
import os 

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('GreenTensor/' + 'reflected/extra','')
#print('Importar modulos necesarios para este codigo')

try:
    sys.path.insert(1, path_basic)
    from graphene_sigma import sigma_DL
except ModuleNotFoundError:
    print('graphene_sigma.py no se encuentra en ' + path_basic)
    
try:
    sys.path.insert(1, path_constants)
    from constants import constantes
except ModuleNotFoundError:
    print('constants.py no se encuentra en ' + path_constants)

pi,hb,c,alfac,mu1,mu2 = constantes()
aux = c*hb

#%%


def functions_GreenTensor(omegac,epsi1,epsi2,hbmu,hbgama,alpha_parallel,x,y,z,xD,yD,zD,zp):
    """
    Parameters
    ----------
    omegac : omega/c = k0 en 1/micrometros    
    epsi1 : epsilon del medio de arriba del plano
    epsi2 : epsilon del medio de abajo del plano
    hbmu : chemical potential in eV  
    hbgama : collision frequency in eV
    alpha_parallel : alpha_parallel (k_parallel/k1)
    x : coordenada x 
    y : coordenada y 
    z : coordenada z
    xD : coordenada x del dipolo 
    yD : coordenada y del dipolo 
    zD : coordenada z del dipolo 
    zp : coordenada zp del plano
    Returns
    -------
    Gxx direct antes de integrar 
    con z hacia abajo (convencion del paper)
    y sin QE
    """
    E = omegac*aux
    k0 = omegac #=omega/c
    # x_y = ky/k0
    n1 = epsi1*mu1
    n2 = epsi2*mu2
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1
    k1_3 = (k0*cte1)**3
    Rbarra =  k1*np.sqrt((x+xD)**2 + (y+yD)**2)
    z_dip_barra = k1*(np.abs(z) + 2*zp + np.abs(zD))
    phi = np.arctan2(np.abs(y+yD),np.abs(x+xD))

    aux2 = n2/n1
    alpha_z1 = np.sqrt(1-alpha_parallel**2) if alpha_parallel<1 else 1j*np.sqrt(alpha_parallel**2-1)
    alpha_z2 = np.sqrt(aux2-alpha_parallel**2) if alpha_parallel<aux2 else 1j*np.sqrt(alpha_parallel**2-aux2)
    
    expB = np.exp(1j*alpha_z1*z_dip_barra) 
    
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    
    rp_num = epsi2*alpha_z1 - epsi1*alpha_z2 + cte1*cond*alpha_z1*alpha_z2
    rp_den = epsi2*alpha_z1 + epsi1*alpha_z2 + cte1*cond*alpha_z1*alpha_z2
    rp = rp_num/rp_den
    
    rs_num = alpha_z1 - alpha_z2 - cond/cte1
    rs_den = alpha_z1 + alpha_z2 + cond/cte1
    rs = rs_num/rs_den    

    J0 = special.jn(0,alpha_parallel*Rbarra) 
    J2 = special.jn(2,alpha_parallel*Rbarra) 

############ todo junto ###################################


    int_re = np.real((J0 + np.cos(2*phi)*J2)*expB*alpha_parallel*(alpha_z1*rp + rs/alpha_z1))
    int_im = np.imag((J0 + np.cos(2*phi)*J2)*expB*alpha_parallel*(alpha_z1*rp + rs/alpha_z1))
 

    rta = int_re + 1j*int_im
    
    cte = 0.5*1j*k1_3 
    
    return rta*cte

#%%
