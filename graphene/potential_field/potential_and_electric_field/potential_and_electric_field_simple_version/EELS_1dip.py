#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila
"""
import numpy as np
import sys
import os 
from scipy import special
from scipy import integrate

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/potential_field/potential_and_electric_field/potential_and_electric_field_simple_version','')
#print('Importar modulos necesarios para este codigo')



try:
    sys.path.insert(1, path_constants)
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

limitt = 100

#%%

def EELS_ana_f(omegac,epsi1,epsi2,hbmu,hbgama,px,py,pz,int_v,b,zp):     
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
#    k0 = omegac #=omega/c
    # x_y = ky/k0
    n1 = epsi1*mu1
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1


    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    Rp = 2*epsi1/(epsi1 + epsi2)
    alfa_p = 1j*(epsi1 + epsi2)/(cond*cte1)
    kp = alfa_p*k1
    kp_3 = kp**3

    if int_v < np.real(alfa_p):
        aux0 = (kp**2 - (omegac*int_v)**2)**(-1/2)
        aux1 = np.cos(np.arcsin(omegac*int_v/kp))/aux0

    else:
        aux0 = ((omegac*int_v)**2 - kp**2 )**(-1/2)
        aux1 = -kp/(aux0*(omegac*int_v +  aux0))
    
    expo = np.exp(-kp*(-b + 2*zp))

    return 2*Rp*int_v*kp_3*(pz*np.sign(pz)*aux0 + px*aux1  )*expo




#%%

def EELS_ana_f_dir(omegac,epsi1,epsi2,px,py,pz,int_v,b):     
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
    
    arg = omegac*int_v*np.abs(b)
    K0 = special.kn(0,arg)
    K1 = special.kn(1,arg)
    
    
    omega = omegac*c
    cte_aux = (omegac*int_v)**2
    
    return np.real((px*K0 - 1j*pz*np.sign(b)*K1)*cte_aux)

#%%


def EELS_num_f(omegac,epsi1,epsi2,hbmu,hbgama,px,py,pz,int_v,b,zp):     
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
#    k0 = omegac #=omega/c
    # x_y = ky/k0
    n1 = epsi1*mu1
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1


    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    Rp = 2*epsi1/(epsi1 + epsi2)
    alfa_p = 1j*(epsi1 + epsi2)/(cond*cte1)
#    alfa_p = np.real(alfa_p)
 
    
    kp = alfa_p*k1
    kp_3 = kp**3
    


    H0 = lambda u : special.hankel1(0,u*alfa_p)
    H1 = lambda u : special.hankel1(1,u*alfa_p)

    cota_inf, cota_sup = 0.01/k1, 600/k1

    term0_re_f = lambda u: np.real(H0(u)*np.cos(int_v*u))
    term0_im_f = lambda u: np.imag(H0(u)*np.cos(int_v*u))
    
    term0_int_re,err = integrate.quad(term0_re_f, cota_inf, cota_sup, limit = limitt) 
    term0_int_im,err = integrate.quad(term0_im_f, cota_inf, cota_sup, limit = limitt) 



    term1_re_f = lambda u: np.real(H1(u)*np.cos(int_v*u))
    term1_im_f = lambda u: np.imag(H1(u)*np.cos(int_v*u))
    
    term1_int_re,err = integrate.quad(term1_re_f, cota_inf, cota_sup, limit = limitt) 
    term1_int_im,err = integrate.quad(term1_im_f, cota_inf, cota_sup, limit = limitt) 



   
    aux0 = term0_int_re + 1j*term0_int_im
    aux1 = term1_int_re + 1j*term1_int_im
    
    expo = np.exp(-kp*(-b + 2*zp))

    return 2*Rp*int_v*kp_3*(pz*np.sign(pz)*aux0 + px*aux1  )*expo


#%%