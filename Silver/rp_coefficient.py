#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila
"""
import numpy as np
import sys
import os 

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
#print('Importar modulos necesarios para este codigo')

try:
    sys.path.insert(1, path_basic)
    from Silver_PP import Silver_lambda_p, Silver_Rp, epsilon_m
except ModuleNotFoundError:
    print('graphene_sigma.py no se encuentra en ' + path_basic)


try:
    sys.path.insert(1, path_basic)
    from constants import constantes
except ModuleNotFoundError:
    print('constants.py no se encuentra en ' + path_basic)

pi,hb,c,alfac,mu1,mu2 = constantes()
aux = c*hb

#%%

def rp_pole_aprox(omegac,epsi1,epsi3,d_nano,k_parallel):     # k_parallel in 1/mu
    """    
    Parameters
    ----------
    omegac : omega/c = k0 en 1/micrometros    
    epsi1 : epsilon del medio de arriba del plano
    epsi2 : epsilon del medio de abajo del plano
    hbmu : chemical potential in eV  
    hbgama : collision frequency in eV
    z : coordenada z
    xD : coordenada x del dipolo 
    yD : coordenada y del dipolo
    zD : coordenada z del dipolo 
    zp : posicion del plano (>0)
    px : coordenada x del dipolo 
    py : coordenada y del dipolo
    pz : coordenada z del dipolo
    Returns
    -------
    px,py,pz en unidades de k*alfa_eff
    """

    E = omegac*aux
    omegac_nano = omegac*1e-3
    d_micros = d_nano*1e-3
    Rp = Silver_Rp(E,epsi1,epsi3)
    lambda_p_v = Silver_lambda_p(E,epsi1,epsi3)*d_micros
    kp = 2*np.pi/lambda_p_v
   
    rp = Rp*kp/(k_parallel - kp)
      
    return rp

#%%

def rp_fresnel_num(omegac,epsi1,epsi3,d_nano,u):     # u in 1/mu
    """    
    Parameters
    ----------
    omegac : omega/c = k0 en 1/micrometros    
    epsi1 : epsilon del medio de arriba del plano
    epsi2 : epsilon del medio de abajo del plano
    hbmu : chemical potential in eV  
    hbgama : collision frequency in eV
    z : coordenada z
    xD : coordenada x del dipolo 
    yD : coordenada y del dipolo
    zD : coordenada z del dipolo 
    zp : posicion del plano (>0)
    px : coordenada x del dipolo 
    py : coordenada y del dipolo
    pz : coordenada z del dipolo
    Returns
    -------
    formula del potencial electric con QE approximation, rp con 
    aproximacion del polo y con aprox de principal value para las integrales
    con rp
    """

    E = omegac*aux
    omegac_nano = omegac*1e-3
    
    d_micros = d_nano*1e-3
    
    epsi_2 = epsilon_m(E)

    kz1 =  np.sqrt(epsi1 - u**2) if (u**2 <= epsi1) else 1j*np.sqrt(u**2 - epsi1)
    kz2 =  np.sqrt(epsi_2 - u**2) if (u**2 <= epsi_2) else 1j*np.sqrt(u**2 - epsi_2)
    kz3 =  np.sqrt(epsi3 - u**2) if (u**2 <= epsi3) else 1j*np.sqrt(u**2 - epsi3)


    if np.imag(kz1) > 0:
        kz1 = kz1
    else:
        kz1 = - kz1
        
    if np.imag(kz2) > 0:
        kz2 = kz2
    else:
        kz2 = - kz2


    if np.imag(kz3) > 0:
        kz3 = kz3
    else:
        kz3 = - kz3


    r12 =  (kz1*epsi_2 - kz2*epsi1)/(kz1*epsi_2 + kz2*epsi1)
    r21 = (kz2*epsi1 - kz1*epsi_2)/(kz2*epsi1 + kz1*epsi_2)
    r23 =  (kz2*epsi3 - kz3*epsi_2)/(kz2*epsi3 + kz3*epsi_2)
    

    exp_fresnel =  np.exp(1j*2*kz2*d_micros)
    
    cte_t = np.sqrt(epsi1*epsi_2)
    t12 =  2*kz1*cte_t/(kz1*epsi_2 + kz2*epsi1)
    t21 =  2*kz2*cte_t/(kz2*epsi1 + kz1*epsi_2)

    rp_num = t12*t21*r23*exp_fresnel
    rp_den = 1 - r21*r23*exp_fresnel
    rp = r12 +  rp_num/rp_den   
 
    return rp

#%%
