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
    from hBn_PP import hBn_lambda_p, hBn_Rp, epsilon_x, epsilon_z
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

def rp_pole_aprox(omegac,epsi1,epsi3,d_nano,k_parallel_nano):     
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
    
#    d_micros = d_nano*1e-3
    Rp = hBn_Rp(E,epsi1,epsi3)
    lambda_p_v = hBn_lambda_p(E,epsi1,epsi3)*d_nano
    kp = 2*np.pi/lambda_p_v
   
    rp = Rp*kp/(k_parallel_nano - kp)
      
    return rp

def rp_pole_aprox_v2(omegac,epsi1,epsi3,d_nano,k_parallel_nano):     
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
    epsi_x = epsilon_x(E)
    epsi_z = epsilon_z(E)
    
    epsi_m = np.sqrt(epsi_x*epsi_z)
#    if np.imag(epsi_m) > 0:
#        epsi_m = epsi_m
#    else:
#        epsi_m = -epsi_m
        
#    d_micros = d_nano*1e-3
    Rp = hBn_Rp(E,epsi1,epsi3)
#    lambda_p_v = hBn_lambda_p(E,epsi1,epsi3)*d_nano

    kp = 1j*(epsi1 + epsi3)*omegac/(4*np.pi*np.imag(epsi_m))
   
    rp = Rp*kp/(k_parallel_nano - kp)
      
    return rp

#%%

def rp_fresnel_num(omegac,epsi1,epsi3,d_nano,k_parallel_nano):     
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
    
    epsi_x = epsilon_x(E)
    epsi_z = epsilon_z(E)
    
    epsi_HBN_par = epsi_x
    epsi_HBN_perp = epsi_z

#    d_micros = d_nano*1e-3

    kz1 = np.sqrt(epsi1*omegac_nano**2 - k_parallel_nano**2) if (k_parallel_nano**2 <= epsi1*omegac_nano**2) else 1j*np.sqrt(k_parallel_nano**2 - epsi1*omegac_nano**2)
    kz2 = np.sqrt(epsi_HBN_par*omegac_nano**2 - (epsi_HBN_par/epsi_HBN_perp)*k_parallel_nano**2) 
    kz3 =  np.sqrt(epsi3*omegac_nano**2 - k_parallel_nano**2) if (k_parallel_nano**2 <= epsi3*omegac_nano**2) else 1j*np.sqrt(k_parallel_nano**2 - epsi3*omegac_nano**2)

###############################################################################################
    if np.imag(kz1) > 0 :  ### el green self image funciona mejor sin este cambio de signo ###
        kz1 = kz1
    else:
        kz1 = -kz1
        
    if np.imag(kz2) > 0 :
        kz2 = kz2
    else:
        kz2 = -kz2
    

    if np.imag(kz3) > 0 :
        kz3 = kz3
    else:
        kz3 = -kz3
################################################################################################
    

    r12 =  (kz1*epsi_HBN_par - kz2*epsi1)/(kz1*epsi_HBN_par + kz2*epsi1)
    r21 =  (kz2*epsi1 - kz1*epsi_HBN_par)/(kz2*epsi1 + kz1*epsi_HBN_par)
    r23 =  (kz2*epsi3 - kz3*epsi_HBN_par)/(kz2*epsi3 + kz3*epsi_HBN_par)
    

    exp_fresnel =  np.exp(1j*2*kz2*d_nano)
    
    cte_t = np.sqrt(epsi1*epsi_HBN_par)
    t12 =  2*kz1*cte_t/(kz1*epsi_HBN_par + kz2*epsi1)
    t21 =  2*kz2*cte_t/(kz2*epsi1 + kz1*epsi_HBN_par)

    rp_num =  t12*t21*r23*exp_fresnel
    rp_den =  1 - r21*r23*exp_fresnel
    rp = r12 +  rp_num/rp_den   
 
    return rp

#%%
