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
    from dipole_moment import dipole_moment_anav2
except ModuleNotFoundError:
    print('dipole_moment.py no se encuentra en ' + path_basic)
    
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

def EELS_ana_f_ind1(omegac,epsi1,epsi2,hbmu,hbgama,int_v,b,zp,omega0,kappa_factor_omega0,kappa_r_factor):     
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
    k_3 = omegac**3

    
    px_v,py_v,pz_v = dipole_moment_anav2(omegac,epsi1,epsi2,hbmu,hbgama,int_v,b,zp,omega0,kappa_factor_omega0,kappa_r_factor)

    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    Rp = 2*epsi1/(epsi1 + epsi2)
    alfa_p = 1j*(epsi1 + epsi2)/(cond*cte1)
    kp = alfa_p*k1
    kp_2 = kp**2


    aux0 = kp*( kp**2 - (omegac*int_v)**2 )**(-1/2)
    aux1 = -kp/(-b + 2*zp) + np.exp(-int_v*omegac*(-b + 2*zp))

#        print('opcion 2')
    expo = np.exp(-kp*(-b + 2*zp))


    rta = Rp*int_v*(pz_v*np.sign(pz_v)*aux1 + px_v*aux0*2*np.pi*(int_v**2)  )*expo*k_3/E
    
    
    return rta/(np.pi*E)





def EELS_PP_f_ind(omegac,epsi1,epsi2,hbmu,hbgama,int_v,b,zp,omega0,kappa_factor_omega0,kappa_r_factor):     
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
    k_3 = omegac**3


    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    Rp = 2*epsi1/(epsi1 + epsi2)
    alfa_p = 1j*(epsi1 + epsi2)/(cond*cte1)
    
#    alfa_p = np.real(alfa_p)
    px_v,py,pz_v = dipole_moment_anav2(omegac,epsi1,epsi2,hbmu,hbgama,int_v,b,zp,omega0,kappa_factor_omega0,kappa_r_factor) 
    
    kp = alfa_p*k1
    kp_3 = kp**3
    
#    rp_num = lambda u: epsi2*1j*u - epsi1*1j*u - u**2
#    rp_den = lambda u: epsi2*1j*u + epsi1*1j*u - u**2
#    rp = lambda u: rp_num(u)/rp_den(u)
    rp = lambda alpha_parallel : Rp*alfa_p/(alpha_parallel-alfa_p)
    exp_electron = lambda alpha_parallel: np.exp(-alpha_parallel*omegac*(np.abs(b) + 2*zp))



    cota_inf, cota_sup = int_v + 0.001, 1000
    
    term_sqr = lambda u : 1/np.sqrt(u**2 - int_v**2)

    term0_re_f = lambda u: np.real(term_sqr(u)*rp(u)*exp_electron(u))
    term0_im_f = lambda u: np.imag(term_sqr(u)*rp(u)*exp_electron(u))
    
    term0_int_re,err = integrate.quad(term0_re_f, int_v + 0.001 , 400, limit = limitt) 
    term0_int_im,err = integrate.quad(term0_im_f, int_v+ 0.001 , 400, limit = limitt) 

    term_sqr2 = lambda u : 1/np.sqrt(int_v**2 - u**2 )

    term1_re_f = lambda u: np.real((1 + int_v*term_sqr2(u))*rp(u)*exp_electron(u)*u)
    term1_im_f = lambda u: np.imag((1 + int_v*term_sqr2(u))*rp(u)*exp_electron(u)*u)
    
    term1_int_re,err = integrate.quad(term1_re_f, 0, int_v - 0.001, limit = limitt) 
    term1_int_im,err = integrate.quad(term1_im_f, 0, int_v - 0.001, limit = limitt) 


   
    aux0 = term0_int_re + 1j*term0_int_im
    aux1 = term1_int_re + 1j*term1_int_im
    

    return 2*int_v*(pz_v*np.sign(pz_v)*aux1 + px_v*aux0*(int_v**2))*k_3/E


#%%
    

def EELS_ana_f_dir(omegac,epsi1,epsi2,hbmu,hbgama,int_v,b,zp,omega0,kappa_factor_omega0,kappa_r_factor):     
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
    arg = omegac*int_v*np.abs(b)
    K0 = special.kn(0,arg)
    K1 = special.kn(1,arg)
    px_v,py_v,pz_v = dipole_moment_anav2(omegac,epsi1,epsi2,hbmu,hbgama,int_v,b,zp,omega0,kappa_factor_omega0,kappa_r_factor)
    
    omega = omegac*c
    cte_aux = (omegac*int_v)**2
    
    rta =  np.real((px_v*K0 - 1j*pz_v*np.sign(b)*K1)*cte_aux)*np.pi
    return rta/(E)




