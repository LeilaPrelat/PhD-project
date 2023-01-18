#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila
"""
import numpy as np
import sys
import os 
from mpmath import *

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

def G1_ana(omegac,epsi1,epsi2,hbmu,hbgama,int_v,zp,y):     
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
 #   n_v1 = int_v/cte1
#    k1_2 = k1**2
   
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    Rp = 2*epsi1/(epsi1 + epsi2)
    alfa_p = 1j*(epsi1 + epsi2)/(cond*cte1)
    kp = alfa_p*k1
    
    
#    nmax =  (np.real(kp) - omega/v)*a/(2*np.pi)    #maximum order 
#    nmax = omegac*(np.real(alfa_p) - int_v)*a/(2*np.pi)
#    nmax = int(nmax)
##    print('k: f2',omegac)
##    
##    print('alfa_p: f2', alfa_p - int_v)
##    print('a/2pi: f2', a/(2*np.pi))
##    print(nmax)
##    print(nmax)
#
#    if nmax > 0:
#        list_dipolos = np.linspace(0,nmax,nmax+1)
#    else: 
#        list_dipolos = np.linspace(nmax,0,-nmax+1)
    a = 0.05
    kx = omegac*int_v
    
#    a = 0.05
#    kx = omegac*int_v 
    
    if np.real(kp) > kx: 
        termy = np.sqrt(kp**2 - kx**2)    
        print('kp>kx: todo ok')
    else:
        print('warning kp<kx:' )
        termy = 1j*np.sqrt(kx**2 - kp**2)    
            
    
    
    exp_electron = np.exp(-alfa_p*2*k1*zp)*np.cos(termy*y) 


    term1 = Rp*4*np.pi*1j*kp*exp_electron/termy

    
###################################################################################################

    return term1 

#%%

def G1_ana_v2(omegac,epsi1,epsi2,hbmu,hbgama,int_v,zp,y):     
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
 #   n_v1 = int_v/cte1
#    k1_2 = k1**2
   
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    Rp = 2*epsi1/(epsi1 + epsi2)
    alfa_p = 1j*(epsi1 + epsi2)/(cond*cte1)
    kp = alfa_p*k1
    
    
#    nmax =  (np.real(kp) - omega/v)*a/(2*np.pi)    #maximum order 
#    nmax = omegac*(np.real(alfa_p) - int_v)*a/(2*np.pi)
#    nmax = int(nmax)
##    print('k: f2',omegac)
##    
##    print('alfa_p: f2', alfa_p - int_v)
##    print('a/2pi: f2', a/(2*np.pi))
##    print(nmax)
##    print(nmax)
#
#    if nmax > 0:
#        list_dipolos = np.linspace(0,nmax,nmax+1)
#    else: 
#        list_dipolos = np.linspace(nmax,0,-nmax+1)

    kx = omegac*int_v 
    
    termy = kx*np.sinh(np.arccosh(kp/kx))
    
    
    exp_electron = np.exp(-alfa_p*2*k1*zp)*np.exp(1j*termy*y) 


    term1 = Rp*2*np.pi*1j*kp*exp_electron

    
###################################################################################################

    return term1 


#%%

def G1_num(omegac,epsi1,epsi2,hbmu,hbgama,int_v,zp,y):     
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

    
    cota_inf = 0.01
    cota_sup = 40
    
    a = 0.05
    alpha_x = int_v 
    
    alpha_parallel = lambda alpha_y: sqrt(alpha_x**2 + alpha_y**2)
    
    exp_electron_f = lambda alpha_y: exp(-alpha_parallel(alpha_y)*2*k1*zp)*cos(alpha_y*k1*y) 



    rp_num = lambda alpha_y: epsi2*1j*alpha_parallel(alpha_y) - epsi1*1j*alpha_parallel(alpha_y) - cte1*cond*(alpha_parallel(alpha_y))**2
    rp_den = lambda alpha_y: epsi2*1j*alpha_parallel(alpha_y) + epsi1*1j*alpha_parallel(alpha_y) - cte1*cond*(alpha_parallel(alpha_y))**2
    rp = lambda alpha_y: rp_num(alpha_y)/rp_den(alpha_y)   
    
    term2_re_f = lambda alpha_y: np.real(rp(alpha_y)*exp_electron_f(alpha_y)/alpha_parallel(alpha_y))
    term2_im_f = lambda alpha_y: np.imag(rp(alpha_y)*exp_electron_f(alpha_y)/alpha_parallel(alpha_y))
    
    term2_int_re = quad(term2_re_f, [cota_inf, cota_sup]) 
    term2_int_im = quad(term2_im_f, [cota_inf, cota_sup]) 

    term2_int = term2_int_re + 1j*term2_int_im
    
    return 2*term2_int  

#%%

def G1_pole_aprox(omegac,epsi1,epsi2,hbmu,hbgama,int_v,zp,y):     
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

    
    cota_inf = 0.01
    cota_sup = 40
        
    a = 0.05
    alpha_x = int_v 
    alpha_parallel = lambda alpha_y: sqrt(alpha_x**2 + alpha_y**2)
    
    exp_electron_f = lambda alpha_y: exp(-alpha_parallel(alpha_y)*2*k1*zp)*cos(alpha_y*k1*y) 

    rp = lambda alpha_y: Rp*alfa_p/(alpha_parallel(alpha_y) - alfa_p)   
    
    term2_re_f = lambda alpha_y: np.real(rp(alpha_y)*exp_electron_f(alpha_y)/alpha_parallel(alpha_y))
    term2_im_f = lambda alpha_y: np.imag(rp(alpha_y)*exp_electron_f(alpha_y)/alpha_parallel(alpha_y))
    
    term2_int_re = quad(term2_re_f, [cota_inf, cota_sup]) 
    term2_int_im = quad(term2_im_f, [cota_inf, cota_sup]) 

    term2_int = term2_int_re + 1j*term2_int_im
    
    return 2*term2_int  


#%%
    

def G2_ana(omegac,epsi1,epsi2,hbmu,hbgama,int_v,zp,y):     
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
 #   n_v1 = int_v/cte1
#    k1_2 = k1**2
   
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    Rp = 2*epsi1/(epsi1 + epsi2)
    alfa_p = 1j*(epsi1 + epsi2)/(cond*cte1)
    kp = alfa_p*k1
    
    
#    nmax =  (np.real(kp) - omega/v)*a/(2*np.pi)    #maximum order 
#    nmax = omegac*(np.real(alfa_p) - int_v)*a/(2*np.pi)
#    nmax = int(nmax)
##    print('k: f2',omegac)
##    
##    print('alfa_p: f2', alfa_p - int_v)
##    print('a/2pi: f2', a/(2*np.pi))
##    print(nmax)
##    print(nmax)
#
#    if nmax > 0:
#        list_dipolos = np.linspace(0,nmax,nmax+1)
#    else: 
#        list_dipolos = np.linspace(nmax,0,-nmax+1)

    a = 0.05
    kx = omegac*int_v 
    
    kp_re = np.real(kp)
    kp_im = np.imag(kp)
    kp_2 = kp_re**2 - kp_im**2 + 2*1j*kp_re*kp_im 
    
    if np.real(kp) > kx: 
        termy = np.sqrt(kp_2 - kx**2)    
        print('kp>kx: todo ok')
    else:
        print('warning kp<kx:' )
        termy = 1j*np.sqrt(kx**2 - kp**2)    
        
            
    termy = np.sqrt(kp_2 - kx**2)    
#    termy = np.sqrt(kp**2 - kx**2)
    exp_electron = np.exp(-alfa_p*2*k1*zp)*np.sin(termy*y) 


    term1 = -Rp*4*np.pi*kp*exp_electron/k1

    
###################################################################################################

    return term1 


#%%

def G2_num(omegac,epsi1,epsi2,hbmu,hbgama,int_v,zp,y):     
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

    
    cota_inf = 0.01
    cota_sup = 40
    

    a = 0.05   
    alpha_x = int_v 
    alpha_parallel = lambda alpha_y: sqrt(alpha_x**2 + alpha_y**2)
    
    exp_electron_f = lambda alpha_y: exp(-alpha_parallel(alpha_y)*2*k1*zp)*1j*sin(alpha_y*k1*y) 



    rp_num = lambda alpha_y: epsi2*1j*alpha_parallel(alpha_y) - epsi1*1j*alpha_parallel(alpha_y) - cte1*cond*(alpha_parallel(alpha_y))**2
    rp_den = lambda alpha_y: epsi2*1j*alpha_parallel(alpha_y) + epsi1*1j*alpha_parallel(alpha_y) - cte1*cond*(alpha_parallel(alpha_y))**2
    rp = lambda alpha_y: rp_num(alpha_y)/rp_den(alpha_y)   
    
    term2_re_f = lambda alpha_y: np.real(2*alpha_y*rp(alpha_y)*exp_electron_f(alpha_y)/alpha_parallel(alpha_y))
    term2_im_f = lambda alpha_y: np.imag(2*alpha_y*rp(alpha_y)*exp_electron_f(alpha_y)/alpha_parallel(alpha_y))
    
    term2_int_re = quad(term2_re_f, [cota_inf, cota_sup]) 
    term2_int_im = quad(term2_im_f, [cota_inf, cota_sup]) 

    term2_int = term2_int_re + 1j*term2_int_im

    
    return term2_int  

#%%

def G2_pole_aprox(omegac,epsi1,epsi2,hbmu,hbgama,int_v,zp,y):     
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

    
    cota_inf = 0.01
    cota_sup = 40
        
    a = 0.05   
    alpha_x = int_v 
    alpha_parallel = lambda alpha_y: sqrt(alpha_x**2 + alpha_y**2)
    
    exp_electron_f = lambda alpha_y: exp(-alpha_parallel(alpha_y)*2*k1*zp)*1j*sin(alpha_y*k1*y) 

    rp = lambda alpha_y: Rp*alfa_p/(alpha_parallel(alpha_y) - alfa_p)   
    
    term2_re_f = lambda alpha_y: np.real(alpha_y*2*rp(alpha_y)*exp_electron_f(alpha_y)/alpha_parallel(alpha_y))
    term2_im_f = lambda alpha_y: np.imag(alpha_y*2*rp(alpha_y)*exp_electron_f(alpha_y)/alpha_parallel(alpha_y))
    
    term2_int_re = quad(term2_re_f, [cota_inf, cota_sup]) 
    term2_int_im = quad(term2_im_f, [cota_inf, cota_sup]) 

    term2_int = term2_int_re + 1j*term2_int_im
    
    return term2_int  

#%%
    
def G3_num(omegac,epsi1,epsi2,hbmu,hbgama,int_v,zp,y):     
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

    
    cota_inf = 0.01
    cota_sup = 40
    

    a = 0.05   
    alpha_x = int_v 
    alpha_parallel = lambda alpha_y: sqrt(alpha_x**2 + alpha_y**2)
    
    exp_electron_f = lambda alpha_y: exp(-alpha_parallel(alpha_y)*2*k1*zp)*cos(alpha_y*k1*y) 



    rp_num = lambda alpha_y: epsi2*1j*alpha_parallel(alpha_y) - epsi1*1j*alpha_parallel(alpha_y) - cte1*cond*(alpha_parallel(alpha_y))**2
    rp_den = lambda alpha_y: epsi2*1j*alpha_parallel(alpha_y) + epsi1*1j*alpha_parallel(alpha_y) - cte1*cond*(alpha_parallel(alpha_y))**2
    rp = lambda alpha_y: rp_num(alpha_y)/rp_den(alpha_y)   
    
    term2_re_f = lambda alpha_y: np.real(alpha_y*rp(alpha_y)*exp_electron_f(alpha_y))
    term2_im_f = lambda alpha_y: np.imag(alpha_y*rp(alpha_y)*exp_electron_f(alpha_y))
    
    term2_int_re = quad(term2_re_f, [cota_inf, cota_sup]) 
    term2_int_im = quad(term2_im_f, [cota_inf, cota_sup]) 

    term2_int = term2_int_re + 1j*term2_int_im
    
    return 2*term2_int  

#%%
    

def G3_pole_aprox(omegac,epsi1,epsi2,hbmu,hbgama,int_v,zp,y):     
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

    
    cota_inf = 0.01
    cota_sup = 40
        
    a = 0.05   
    alpha_x = int_v 
    alpha_parallel = lambda alpha_y: sqrt(alpha_x**2 + alpha_y**2)
    
    exp_electron_f = lambda alpha_y: exp(-alpha_parallel(alpha_y)*2*k1*zp)*cos(alpha_y*k1*y) 

    rp = lambda alpha_y: Rp*alfa_p/(alpha_parallel(alpha_y) - alfa_p)   
    
    term2_re_f = lambda alpha_y: np.real(alpha_y*rp(alpha_y)*exp_electron_f(alpha_y))
    term2_im_f = lambda alpha_y: np.imag(alpha_y*rp(alpha_y)*exp_electron_f(alpha_y))
    
    term2_int_re = quad(term2_re_f, [cota_inf, cota_sup]) 
    term2_int_im = quad(term2_im_f, [cota_inf, cota_sup]) 

    term2_int = term2_int_re + 1j*term2_int_im
    
    return 2*term2_int  


#%%
    

def G3_ana(omegac,epsi1,epsi2,hbmu,hbgama,int_v,zp,y):     
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
 #   n_v1 = int_v/cte1
#    k1_2 = k1**2
   
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    Rp = 2*epsi1/(epsi1 + epsi2)
    alfa_p = 1j*(epsi1 + epsi2)/(cond*cte1)
    kp = alfa_p*k1
    
    
#    nmax =  (np.real(kp) - omega/v)*a/(2*np.pi)    #maximum order 
#    nmax = omegac*(np.real(alfa_p) - int_v)*a/(2*np.pi)
#    nmax = int(nmax)
##    print('k: f2',omegac)
##    
##    print('alfa_p: f2', alfa_p - int_v)
##    print('a/2pi: f2', a/(2*np.pi))
##    print(nmax)
##    print(nmax)
#
#    if nmax > 0:
#        list_dipolos = np.linspace(0,nmax,nmax+1)
#    else: 
#        list_dipolos = np.linspace(nmax,0,-nmax+1)

    a = 0.05
    kx = omegac*int_v 
    
    if np.real(kp) > kx: 
        termy = np.sqrt(kp**2 - kx**2)    
        print('kp>kx: todo ok')
    else:
        print('warning kp<kx:' )
        termy = 1j*np.sqrt(kx**2 - kp**2)    
        
            
    
#    termy = np.sqrt(np.real(kp)**2 - kx**2)
    exp_electron = np.exp(-alfa_p*2*k1*zp)*np.cos(termy*y) 


    term1 = Rp*4*np.pi*1j*kp*alfa_p*exp_electron/termy
    
###################################################################################################

    return term1 


#%%