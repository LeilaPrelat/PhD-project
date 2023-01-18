#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila
"""
import numpy as np
import sys
import os 
from scipy import integrate

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/potential_field/infinite_dipoles','')
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

def phi_many_dipoles_ana_n(omegac,epsi1,epsi2,hbmu,hbgama,int_v,zp,x,y,z,px,py,pz,a,n):     
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
    n1 = epsi1*mu1
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1
   
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    Rp = 2*epsi1/(epsi1 + epsi2)
    alfa_p = 1j*(epsi1 + epsi2)/(cond*cte1)
    kp = alfa_p*k1

   
    kp_re = np.real(kp)
    kp_im = np.imag(kp)
    kp_2 = kp_re**2 - kp_im**2 + 2*1j*kp_re*kp_im 
    

#    list_dipoles = np.linspace(-Nmax,Nmax,2*Nmax + 1)
#            
    kx = omegac*int_v + 2*np.pi*n/a     
    expo_kx = np.exp(1j*kx*x)
    
    if np.real(kp) > kx: 
        termy = np.sqrt(kp_2 - kx**2)    
#        print('kp>kx: todo ok')
    else:
#        print('warning kp<kx:' )
        termy = 1j*np.sqrt(kx**2 - kp_2)    

    alfa_pv2 = np.sqrt(alfa_p**2)
    kp_v2 = np.sqrt(kp**2)
    
    term_alfa = (1 + alfa_p/alfa_pv2)
    
    term_den = np.sqrt(kp**2 - kx**2)
    
    exp_electron = np.exp(-alfa_pv2*k1*(2*zp - z))*np.exp(1j*termy*y) 

    term1 = Rp*np.pi*1j*kp*exp_electron*term_alfa/term_den

    term2 = Rp*np.pi*1j*kp*exp_electron*term_alfa ## multiplicado por k
    
    
    term3 = Rp*np.pi*1j*kp*exp_electron*(kp_v2 + kp)/term_den  ## multiplicado por k

    final =  1j*px*kx*term1 + 1j*py*term2 - pz*term3


    cte = 4*(np.pi**2)*a
    return final*expo_kx/cte

#%%

def phi_many_dipoles_num_n(omegac,epsi1,epsi2,hbmu,hbgama,int_v,zp,x,y,z,px,py,pz,a,n):     
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
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)  

    
    cota_inf = 1*k1
    cota_sup = 600*k1



#    list_dipoles = np.linspace(-Nmax,Nmax,2*Nmax + 1)
#    
#    tot0 = 0
#    for n in list_dipoles:

    alpha_x = int_v + 2*np.pi*n/(a*omegac)
    
    kx = alpha_x*omegac
    
    expo_kx = np.exp(1j*alpha_x*k1*x)
    
    alpha_parallel = lambda alpha_y: np.sqrt(alpha_x**2 + alpha_y**2)
    
    exp_electron_f = lambda alpha_y: np.exp(-alpha_parallel(alpha_y)*k1*(2*zp - z))*np.exp(1j*alpha_y*k1*np.abs(y)) 



    rp_num = lambda alpha_y: epsi2*1j*alpha_parallel(alpha_y) - epsi1*1j*alpha_parallel(alpha_y) - cte1*cond*(alpha_parallel(alpha_y))**2
    rp_den = lambda alpha_y: epsi2*1j*alpha_parallel(alpha_y) + epsi1*1j*alpha_parallel(alpha_y) - cte1*cond*(alpha_parallel(alpha_y))**2
    rp = lambda alpha_y: rp_num(alpha_y)/rp_den(alpha_y)   
    
    term1_re_f = lambda alpha_y: np.real(rp(alpha_y)*exp_electron_f(alpha_y)/alpha_parallel(alpha_y))
    term1_im_f = lambda alpha_y: np.imag(rp(alpha_y)*exp_electron_f(alpha_y)/alpha_parallel(alpha_y))
    
    term1_int_re,err = integrate.quad(term1_re_f, cota_inf, cota_sup) 
    term1_int_im,err = integrate.quad(term1_im_f, cota_inf, cota_sup) 

    term1_int = (term1_int_re + 1j*term1_int_im)*1j*px*kx
    
    term2_re_f = lambda alpha_y: np.real(alpha_y*rp(alpha_y)*exp_electron_f(alpha_y)/alpha_parallel(alpha_y))
    term2_im_f = lambda alpha_y: np.imag(alpha_y*rp(alpha_y)*exp_electron_f(alpha_y)/alpha_parallel(alpha_y))
    
    term2_int_re,err = integrate.quad(term2_re_f, cota_inf, cota_sup) 
    term2_int_im,err = integrate.quad(term2_im_f, cota_inf, cota_sup) 
        
    term2_int = (term2_int_re + 1j*term2_int_im)*1j*py*omegac
    
    
    term3_re_f = lambda alpha_y: np.real(rp(alpha_y)*exp_electron_f(alpha_y))
    term3_im_f = lambda alpha_y: np.imag(rp(alpha_y)*exp_electron_f(alpha_y))
    
    term3_int_re,err = integrate.quad(term3_re_f, cota_inf, cota_sup) 
    term3_int_im,err = integrate.quad(term3_im_f, cota_inf, cota_sup) 
    
    term3_int = -(term3_int_re + 1j*term3_int_im)*pz*omegac
    
    cte = 4*(np.pi**2)*a    
    return expo_kx*(term1_int + term2_int + term3_int)/cte

#%%

def phi_many_dipoles_pole_aprox_n(omegac,epsi1,epsi2,hbmu,hbgama,int_v,zp,x,y,z,px,py,pz,a,n):     
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
#    k1_2 = k1**2
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    Rp = 2*epsi1/(epsi1 + epsi2)
    alfa_p = 1j*(epsi1 + epsi2)/(cond*cte1)

    
    cota_inf = 1*k1
    cota_sup = 600*k1
        
#
#    list_dipoles = np.linspace(-Nmax,Nmax,2*Nmax + 1)
#    
#    tot0 = 0
#    for n in list_dipoles:
        
    alpha_x = int_v + 2*np.pi*n/(a*omegac)
 
    kx = alpha_x*omegac
    expo_kx = np.exp(1j*alpha_x*k1*x)
    
    alpha_parallel = lambda alpha_y: np.sqrt(alpha_x**2 + alpha_y**2)

    rp = lambda alpha_y: Rp*alfa_p/(alpha_parallel(alpha_y) - alfa_p)  
    
    exp_electron_f = lambda alpha_y: np.exp(-alpha_parallel(alpha_y)*k1*(2*zp - z))*np.exp(1j*alpha_y*k1*np.abs(y))  
    
    term1_re_f = lambda alpha_y: np.real(rp(alpha_y)*exp_electron_f(alpha_y)/alpha_parallel(alpha_y))
    term1_im_f = lambda alpha_y: np.imag(rp(alpha_y)*exp_electron_f(alpha_y)/alpha_parallel(alpha_y))
    
    term1_int_re,err = integrate.quad(term1_re_f, cota_inf, cota_sup) 
    term1_int_im,err = integrate.quad(term1_im_f, cota_inf, cota_sup) 

    term1_int = (term1_int_re + 1j*term1_int_im)*1j*px*kx
    
    term2_re_f = lambda alpha_y: np.real(alpha_y*rp(alpha_y)*exp_electron_f(alpha_y)/alpha_parallel(alpha_y))
    term2_im_f = lambda alpha_y: np.imag(alpha_y*rp(alpha_y)*exp_electron_f(alpha_y)/alpha_parallel(alpha_y))
    
    term2_int_re,err = integrate.quad(term2_re_f, cota_inf, cota_sup) 
    term2_int_im,err = integrate.quad(term2_im_f, cota_inf, cota_sup) 
        
    term2_int = (term2_int_re + 1j*term2_int_im)*1j*py*omegac
    
    
    term3_re_f = lambda alpha_y: np.real(rp(alpha_y)*exp_electron_f(alpha_y))
    term3_im_f = lambda alpha_y: np.imag(rp(alpha_y)*exp_electron_f(alpha_y))
    
    term3_int_re,err = integrate.quad(term3_re_f, cota_inf, cota_sup) 
    term3_int_im,err = integrate.quad(term3_im_f, cota_inf, cota_sup) 
    
    term3_int = -(term3_int_re + 1j*term3_int_im)*pz*omegac
    
    cte = 4*(np.pi**2)*a    
    return expo_kx*(term1_int + term2_int + term3_int)/cte


#%%
    
    