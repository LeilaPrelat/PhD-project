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
path_constants =  path_basic.replace('/potential_field/potential_and_electric_field/potential_and_electric_field_with_dipole_moment_formula','')
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

#%%

def dipole_moment_ana(omegac,epsi1,epsi2,hbmu,hbgama,int_v,px,py,pz,b,zp):     
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
#    k1_2 = (k0*cte1)**2
 #   n_v1 = int_v/cte1

    # alffa_eff = alpha_function_ana(omegac,epsi1,epsi2,hbmu,hbgama,px,py,pz,zp,omega0,kappa_factor_omega0,kappa_r_factor) 

    charge_electron = 4.806e-10*int_v/(2*np.pi)
    charge_electron = int_v/(2*np.pi)
    

    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    Rp = 2*epsi1/(epsi1 + epsi2)
    alfa_p = 1j*(epsi1 + epsi2)/(cond*cte1)
    kp = alfa_p*k1

      
    arg = np.abs(b)*omegac*int_v
    K1 = special.kn(1,arg)
    
    # kx = omegac*int_v
    # expo = np.exp(-np.sqrt(kx**2 + kp**2)*(-np.abs(b) + 2*zp))
    expo = np.exp(-kp*(-np.abs(b) + 2*zp))
    
#    den = np.sqrt(kx**2 + kp**2) - kp
    
    px = charge_electron*1j*(omegac*int_v*K1 - 2*np.pi*1j*Rp*kp*expo)
    
    py = px
    
    pz = charge_electron*(-omegac*int_v*K1 - 2*np.pi*1j*Rp*kp*expo)
    
    return px,py,pz


#%%

def dipole_moment_num(omegac,epsi1,epsi2,hbmu,hbgama,int_v,px,py,pz,b,zp):     
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
    k0 = omegac #=omega/c
    # x_y = ky/k0
#    n1 = epsi1*mu1
#    cte1 = np.sqrt(n1)
    # k1 = omegac*cte1
#    k1_2 = (k0*cte1)**2
 #   n_v1 = int_v/cte1
    
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama) #no es necesario dividir por c 

    # alffa_eff = alpha_function_num(omegac,epsi1,epsi2,hbmu,hbgama,px,py,pz,zp,omega0,kappa_factor_omega0,kappa_r_factor)
    
    charge_electron = 4.806e-10*int_v/(2*np.pi)
    charge_electron = int_v/(2*np.pi)
    
    cota_inf = 1
    cota_sup = 1001*k0 
    kx = int_v*omegac    
    
#    term6 = np.sign(z)*pz*Rp*kp_2*J0*exp_electron
    kparallel = lambda u: np.sqrt(kx**2 + u**2)
    
    
    rp_num = lambda u: epsi2*1j*kparallel(u) - epsi1*1j*kparallel(u) - cond*kparallel(u)**2/k0
    rp_den = lambda u: epsi2*1j*kparallel(u) + epsi1*1j*kparallel(u) - cond*kparallel(u)**2/k0
    rp = lambda u: rp_num(u)/rp_den(u)
 
      
    arg = np.abs(b)*omegac*int_v
    K1 = special.kn(1,arg)
    



    expo = lambda u: np.exp(-np.sqrt(kx**2 + u**2)*(-np.abs(b) + 2*zp))
    
    int_f_re = lambda u: np.real(rp(u)*expo(u)/kparallel(u))
    int_f_im = lambda u: np.imag(rp(u)*expo(u)/kparallel(u))
    
    INT_re,err = integrate.quad(int_f_re, cota_inf, cota_sup) 
    INT_im,err = integrate.quad(int_f_im, cota_inf, cota_sup) 
    
    INT = INT_re + 1j*INT_im
    
    px = charge_electron*1j*(omegac*int_v*K1 - INT)
    
    py = px
    
    pz = charge_electron*(-omegac*int_v*K1 - INT)   
    
    return px,py,pz

#%%