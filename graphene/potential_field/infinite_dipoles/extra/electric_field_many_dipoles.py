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
from scipy import special

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

def G1_ana(omegac,epsi1,epsi2,hbmu,hbgama,int_v,zp,a,y):     
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


#    nmax =  (np.real(kp) - omegac*int_v)*a/(2*np.pi)    #maximum order 
    nmax = omegac*(np.real(alfa_p) - int_v)*a/(2*np.pi)
    nmax = int(nmax)    
#    a = 0.05
#    kx = omegac*int_v
    if nmax > 0:
        list_dipolos = np.linspace(0,nmax,nmax+1)
    else: 
        list_dipolos = np.linspace(0,-nmax,-nmax+1)
        
    int_tot0 = 0         
    factor = a*(2*np.pi)**2
    for n in list_dipolos:
        
        kx = omegac*int_v + 2*n*pi/a
            
        if np.real(kp) > kx: 
            termy = np.sqrt(kp**2 - kx**2)    
            print('kp>kx: todo ok')
        else:
            print('warning kp<kx:' )
            termy = 1j*np.sqrt(kx**2 - kp**2)    
            
    
    
        exp_electron = np.exp(-alfa_p*2*k1*zp)*np.cos(termy*y) 


        term1 = Rp*4*np.pi*1j*kp*exp_electron/termy
        int_tot0 = int_tot0 + term1
    
###################################################################################################

    return int_tot0/factor

#%%

def G1_num(omegac,epsi1,epsi2,hbmu,hbgama,int_v,zp,a,y):     
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

    
    cota_inf = 0.1*k1
    cota_sup = 40*k1
    
    N = 20
    list_dipolos = np.linspace(0,2*N,2*N + 1)
    int_tot0 = 0 
    
    factor = a*(2*np.pi)**2
    for n in list_dipolos:
        
        alpha_x = int_v + 2*pi*n/(a*omegac)
    
        alpha_parallel = lambda alpha_y: np.sqrt(alpha_x**2 + alpha_y**2)
        
        exp_electron_f = lambda alpha_y: np.exp(-alpha_parallel(alpha_y)*2*k1*zp)*np.cos(alpha_y*k1*y) 
    
        rp_num = lambda alpha_y: epsi2*1j*alpha_parallel(alpha_y) - epsi1*1j*alpha_parallel(alpha_y) - cte1*cond*(alpha_parallel(alpha_y))**2
        rp_den = lambda alpha_y: epsi2*1j*alpha_parallel(alpha_y) + epsi1*1j*alpha_parallel(alpha_y) - cte1*cond*(alpha_parallel(alpha_y))**2
        rp = lambda alpha_y: rp_num(alpha_y)/rp_den(alpha_y)   
        
        term2_re_f = lambda alpha_y: np.real(rp(alpha_y)*exp_electron_f(alpha_y)/alpha_parallel(alpha_y))
        term2_im_f = lambda alpha_y: np.imag(rp(alpha_y)*exp_electron_f(alpha_y)/alpha_parallel(alpha_y))
        
        term2_int_re,err = integrate.quad(term2_re_f, cota_inf, cota_sup) 
        term2_int_im,err = integrate.quad(term2_im_f, cota_inf, cota_sup) 
    
        term2_int = term2_int_re + 1j*term2_int_im
        
        int_tot0 = int_tot0 + term2_int
        
    return 2*int_tot0/factor 

#%%

def G1_pole_aprox(omegac,epsi1,epsi2,hbmu,hbgama,int_v,zp,a,y):     
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
    
    nmax = omegac*(np.real(alfa_p) - int_v)*a/(2*np.pi)
    nmax = int(nmax)    
    if nmax > 0:
        list_dipolos = np.linspace(0,nmax,nmax+1)
    else: 
        list_dipolos = np.linspace(0,-nmax,-nmax+1)
        
    cota_inf = 0.1*k1
    cota_sup = 40*k1
        
    factor = a*(2*np.pi)**2
    int_tot0 = 0
    for n in list_dipolos:
        
        alpha_x = int_v + 2*pi*n/(a*omegac)
        alpha_parallel = lambda alpha_y: np.sqrt(alpha_x**2 + alpha_y**2)
        
        exp_electron_f = lambda alpha_y: np.exp(-alpha_parallel(alpha_y)*2*k1*zp)*np.cos(alpha_y*k1*y) 
    
        rp = lambda alpha_y: Rp*alfa_p/(alpha_parallel(alpha_y) - alfa_p)   
        
        term2_re_f = lambda alpha_y: np.real(rp(alpha_y)*exp_electron_f(alpha_y)/alpha_parallel(alpha_y))
        term2_im_f = lambda alpha_y: np.imag(rp(alpha_y)*exp_electron_f(alpha_y)/alpha_parallel(alpha_y))
        
        term2_int_re,err1 = integrate.quad(term2_re_f, cota_inf, cota_sup) 
        term2_int_im,err2 = integrate.quad(term2_im_f, cota_inf, cota_sup) 
    
        term2_int = term2_int_re + 1j*term2_int_im
        
        int_tot0 = int_tot0 + term2_int   
        
        print('G1:', err1/term2_int_re, err2/term2_int_im)
    
    return 2*int_tot0/factor


#%%
    

def G2_ana(omegac,epsi1,epsi2,hbmu,hbgama,int_v,zp,a,y):     
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
    
    nmax = omegac*(np.real(alfa_p) - int_v)*a/(2*np.pi)
    nmax = int(nmax)    
    if nmax > 0:
        list_dipolos = np.linspace(0,nmax,nmax+1)
    else: 
        list_dipolos = np.linspace(0,-nmax,-nmax+1)
        
    factor = a*(2*np.pi)**2
    int_tot0 = 0   
    for n in list_dipolos:
        kx = omegac*int_v + 2*pi*n/a
        
        if np.real(kp) > kx: 
            termy = np.sqrt(kp**2 - kx**2)    
            print('kp>kx: todo ok')
        else:
            print('warning kp<kx:' )
            termy = 1j*np.sqrt(kx**2 - kp**2)    
            
                
        
    #    termy = np.sqrt(kp**2 - kx**2)
        exp_electron = np.exp(-alfa_p*2*k1*zp)*np.sin(termy*y) 
    
    
        term1 = -Rp*4*np.pi*kp*exp_electron/k1
        
        int_tot0 = int_tot0 + term1

    
###################################################################################################

    return int_tot0/factor


#%%

def G2_num(omegac,epsi1,epsi2,hbmu,hbgama,int_v,zp,a,y):     
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

    
    cota_inf = 0.1*k1
    cota_sup = 40*k1
    

    N = 20
    list_dipolos = np.linspace(0,N,2*N+1)
        
    factor = a*(2*np.pi)**2
    int_tot0 = 0   
    
    for n in list_dipolos: 
        
        alpha_x = int_v + 2*pi*n/(a*omegac)
        
        alpha_parallel = lambda alpha_y: np.sqrt(alpha_x**2 + alpha_y**2)
        
        exp_electron_f = lambda alpha_y: np.exp(-alpha_parallel(alpha_y)*2*k1*zp)*1j*np.sin(alpha_y*k1*y) 
    
    
    
        rp_num = lambda alpha_y: epsi2*1j*alpha_parallel(alpha_y) - epsi1*1j*alpha_parallel(alpha_y) - cte1*cond*(alpha_parallel(alpha_y))**2
        rp_den = lambda alpha_y: epsi2*1j*alpha_parallel(alpha_y) + epsi1*1j*alpha_parallel(alpha_y) - cte1*cond*(alpha_parallel(alpha_y))**2
        rp = lambda alpha_y: rp_num(alpha_y)/rp_den(alpha_y)   
        
        term2_re_f = lambda alpha_y: np.real(2*alpha_y*rp(alpha_y)*exp_electron_f(alpha_y)/alpha_parallel(alpha_y))
        term2_im_f = lambda alpha_y: np.imag(2*alpha_y*rp(alpha_y)*exp_electron_f(alpha_y)/alpha_parallel(alpha_y))
        
        term2_int_re,err = integrate.quad(term2_re_f, cota_inf, cota_sup) 
        term2_int_im,err = integrate.quad(term2_im_f, cota_inf, cota_sup) 
    
        term2_int = term2_int_re + 1j*term2_int_im

        int_tot0 = int_tot0 + term2_int
        
    return int_tot0/factor  

#%%

def G2_pole_aprox(omegac,epsi1,epsi2,hbmu,hbgama,int_v,zp,a,y):     
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


    nmax = omegac*(np.real(alfa_p) - int_v)*a/(2*np.pi)
    nmax = int(nmax)    
    if nmax > 0:
        list_dipolos = np.linspace(0,nmax,nmax+1)
    else: 
        list_dipolos = np.linspace(0,-nmax,-nmax+1)
   
    cota_inf = 0.1*k1
    cota_sup = 40*k1

    factor = a*(2*np.pi)**2
    int_tot0 = 0  
    for n in list_dipolos:
 
        alpha_x = int_v + 2*pi*n/(a*omegac)
        alpha_parallel = lambda alpha_y: np.sqrt(alpha_x**2 + alpha_y**2)
        
        exp_electron_f = lambda alpha_y: np.exp(-alpha_parallel(alpha_y)*2*k1*zp)*1j*np.sin(alpha_y*k1*y) 
    
        rp = lambda alpha_y: Rp*alfa_p/(alpha_parallel(alpha_y) - alfa_p)   
        
        term2_re_f = lambda alpha_y: np.real(alpha_y*2*rp(alpha_y)*exp_electron_f(alpha_y)/alpha_parallel(alpha_y))
        term2_im_f = lambda alpha_y: np.imag(alpha_y*2*rp(alpha_y)*exp_electron_f(alpha_y)/alpha_parallel(alpha_y))
        
        term2_int_re,err1 = integrate.quad(term2_re_f, cota_inf, cota_sup) 
        term2_int_im,err2 = integrate.quad(term2_im_f, cota_inf, cota_sup) 
    
        term2_int = term2_int_re + 1j*term2_int_im
        
        int_tot0 = term2_int + int_tot0

        print('G2:', err1/term2_int_re, err2/term2_int_im)
    
    return int_tot0/factor  

#%%
    
def G3_num(omegac,epsi1,epsi2,hbmu,hbgama,int_v,zp,a,y):     
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

    
    cota_inf = 0.1*k1
    cota_sup = 40*k1
    

    N = 20
    list_dipolos = np.linspace(0,N,2*N+1)

    factor = a*(2*np.pi)**2
    int_tot0 = 0  
    for n in list_dipolos:
        
        alpha_x = int_v + 2*pi*n/(a*omegac)
        alpha_parallel = lambda alpha_y: np.sqrt(alpha_x**2 + alpha_y**2)
        
        exp_electron_f = lambda alpha_y: np.exp(-alpha_parallel(alpha_y)*2*k1*zp)*np.cos(alpha_y*k1*y) 
    
    
    
        rp_num = lambda alpha_y: epsi2*1j*alpha_parallel(alpha_y) - epsi1*1j*alpha_parallel(alpha_y) - cte1*cond*(alpha_parallel(alpha_y))**2
        rp_den = lambda alpha_y: epsi2*1j*alpha_parallel(alpha_y) + epsi1*1j*alpha_parallel(alpha_y) - cte1*cond*(alpha_parallel(alpha_y))**2
        rp = lambda alpha_y: rp_num(alpha_y)/rp_den(alpha_y)   
        
        term2_re_f = lambda alpha_y: np.real(alpha_y*rp(alpha_y)*exp_electron_f(alpha_y))
        term2_im_f = lambda alpha_y: np.imag(alpha_y*rp(alpha_y)*exp_electron_f(alpha_y))
        
        term2_int_re,err = integrate.quad(term2_re_f, cota_inf, cota_sup) 
        term2_int_im,err = integrate.quad(term2_im_f, cota_inf, cota_sup) 
    
        term2_int = term2_int_re + 1j*term2_int_im
        
        int_tot0 = int_tot0 + term2_int
    
    return 2*int_tot0/factor  

#%%
    

def G3_pole_aprox(omegac,epsi1,epsi2,hbmu,hbgama,int_v,zp,a,y):     
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


    nmax = omegac*(np.real(alfa_p) - int_v)*a/(2*np.pi)
    nmax = int(nmax)    
    if nmax > 0:
        list_dipolos = np.linspace(0,nmax,nmax+1)
    else: 
        list_dipolos = np.linspace(0,-nmax,-nmax+1)
 
    cota_inf = 0.1*k1
    cota_sup = 40*k1


    factor = a*(2*np.pi)**2     
    int_tot0 = 0  
    for n in list_dipolos: 
        
        alpha_x = int_v  + 2*pi*n/(a*omegac)
        alpha_parallel = lambda alpha_y: np.sqrt(alpha_x**2 + alpha_y**2)
        
        exp_electron_f = lambda alpha_y: np.exp(-alpha_parallel(alpha_y)*2*k1*zp)*np.cos(alpha_y*k1*y) 
    
        rp = lambda alpha_y: Rp*alfa_p/(alpha_parallel(alpha_y) - alfa_p)   
        
        term2_re_f = lambda alpha_y: np.real(alpha_y*rp(alpha_y)*exp_electron_f(alpha_y))
        term2_im_f = lambda alpha_y: np.imag(alpha_y*rp(alpha_y)*exp_electron_f(alpha_y))
        
        term2_int_re,err1 = integrate.quad(term2_re_f, cota_inf, cota_sup) 
        term2_int_im,err2 = integrate.quad(term2_im_f, cota_inf, cota_sup) 
    
        term2_int = term2_int_re + 1j*term2_int_im
        
        int_tot0 = int_tot0 + term2_int
        
        
        print('G3:', err1/term2_int_re, err2/term2_int_im)
    
    return 2*int_tot0/factor  


#%%
    

def G3_ana(omegac,epsi1,epsi2,hbmu,hbgama,int_v,zp,a,y):     
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


    nmax = omegac*(np.real(alfa_p) - int_v)*a/(2*np.pi)
    nmax = int(nmax)    
    if nmax > 0:
        list_dipolos = np.linspace(0,nmax,nmax+1)
    else: 
        list_dipolos = np.linspace(0,-nmax,-nmax+1)

    factor = a*(2*np.pi)**2

    int_tot0 = 0  
    for n in list_dipolos: 
        
        kx = omegac*int_v + 2*pi*n/(a*omegac)
        
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
        
        int_tot0 = int_tot0 + term1
        
    return int_tot0/factor 


#%%
    

def Ex_many_dipoles_num(omegac,epsi1,epsi2,hbmu,hbgama,x,y,z,a,zp,int_v,px,py,pz):
    """
    Parameters
    ----------
    omegac : omega/c = k0 en 1/micrometros    
    epsi1 : epsilon del medio de arriba del plano
    epsi2 : epsilon del medio de abajo del plano
    hbmu : chemical potential in eV  
    hbgama : collision frequency in eV
    x : coordenada x 
    y : coordenada y 
    z : coordenada z
    a : distancia en x entre dipolos
    N : cantidad de dipolos
    zp : coordenada zp del plano
    Returns
    -------
    potential creado por muchos dipolos 
    +
    """

#    k0 = omegac #=omega/c
    # x_y = ky/k0
    n1 = epsi1*mu1
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1    
    
    G1 = G1_num(omegac,epsi1,epsi2,hbmu,hbgama,int_v,zp,a,y)
    G2 = G2_num(omegac,epsi1,epsi2,hbmu,hbgama,int_v,zp,a,y)
    G3 = G3_num(omegac,epsi1,epsi2,hbmu,hbgama,int_v,zp,a,y)
    
            
#        rp = Rp*kp/(k_parallel(ky) - kp)
    ################### I1 ################################################
    N = 20
    list_dipolos = np.linspace(0,N,2*N+1)

    factor = a*(2*np.pi)**2

    int_tot0 = 0  
    for n in list_dipolos:
        
        kx = omegac*int_v + 2*pi*n/(a*omegac)            
        exp_x = np.exp(1j*kx*x)
        
        K0 = special.kv(0,kx*y)
        K1 = special.kv(1,kx*y)
        
        term3 = (-1j*px*K0 + 2*py*K1)*kx
        
        term4 = pz*1j*np.sign(z)/y
        
        
        term_total = term3 + term4 + G1*1j*px*kx + G2*1j*py*k1  - G3*pz*k1
        
        int_tot0 = term_total + int_tot0*exp_x
    

    return int_tot0/factor
    
    
#%%
    

def Ex_many_dipoles_pole_aprox(omegac,epsi1,epsi2,hbmu,hbgama,x,y,z,a,zp,int_v,px,py,pz):
    """
    Parameters
    ----------
    omegac : omega/c = k0 en 1/micrometros    
    epsi1 : epsilon del medio de arriba del plano
    epsi2 : epsilon del medio de abajo del plano
    hbmu : chemical potential in eV  
    hbgama : collision frequency in eV
    x : coordenada x 
    y : coordenada y 
    z : coordenada z
    a : distancia en x entre dipolos
    N : cantidad de dipolos
    zp : coordenada zp del plano
    Returns
    -------
    potential creado por muchos dipolos 
    +
    """
    E = omegac*aux
#    k0 = omegac #=omega/c
    # x_y = ky/k0
    n1 = epsi1*mu1
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1    
    
    G1 = G1_pole_aprox(omegac,epsi1,epsi2,hbmu,hbgama,int_v,zp,a,y)
    G2 = G2_pole_aprox(omegac,epsi1,epsi2,hbmu,hbgama,int_v,zp,a,y)
    G3 = G3_pole_aprox(omegac,epsi1,epsi2,hbmu,hbgama,int_v,zp,a,y)
    
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    alfa_p = 1j*(epsi1 + epsi2)/(cond*cte1)
#        rp = Rp*kp/(k_parallel(ky) - kp)
    ################### I1 ################################################
    nmax = omegac*(int_v - np.real(alfa_p))*a/(2*np.pi)
    nmax = int(nmax)    
    print(nmax)
    if nmax > 0:
        list_dipolos = np.linspace(0,nmax,nmax+1)
    else: 
        list_dipolos = np.linspace(0,-nmax,-nmax+1)

    factor = a*(2*np.pi)**2

    int_tot0 = 0  
    for n in list_dipolos:
        
        kx = omegac*int_v + 2*pi*n/(a*omegac)            
        exp_x = np.exp(1j*kx*x)
        
        K0 = special.kv(0,kx*y)
        K1 = special.kv(1,kx*y)
        
        term3 = (-1j*px*K0 + 2*py*K1)*kx
        
        term4 = pz*1j*np.sign(z)/y
        
        
        term_total = term3 + term4 + G1*1j*px*kx + G2*1j*py*k1  - G3*pz*k1
        
        int_tot0 = term_total + int_tot0*exp_x
    

    return int_tot0/factor    
    
 #%%   
    
def Ex_many_dipoles_ana(omegac,epsi1,epsi2,hbmu,hbgama,x,y,z,a,zp,int_v,px,py,pz):
    """
    Parameters
    ----------
    omegac : omega/c = k0 en 1/micrometros    
    epsi1 : epsilon del medio de arriba del plano
    epsi2 : epsilon del medio de abajo del plano
    hbmu : chemical potential in eV  
    hbgama : collision frequency in eV
    x : coordenada x 
    y : coordenada y 
    z : coordenada z
    a : distancia en x entre dipolos
    N : cantidad de dipolos
    zp : coordenada zp del plano
    Returns
    -------
    potential creado por muchos dipolos 
    +
    """
    E = omegac*aux
#    k0 = omegac #=omega/c
    # x_y = ky/k0
    n1 = epsi1*mu1
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1    
    
    G1 = G1_ana(omegac,epsi1,epsi2,hbmu,hbgama,int_v,zp,a,y)
    G2 = G2_ana(omegac,epsi1,epsi2,hbmu,hbgama,int_v,zp,a,y)
    G3 = G3_ana(omegac,epsi1,epsi2,hbmu,hbgama,int_v,zp,a,y)
    
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    alfa_p = 1j*(epsi1 + epsi2)/(cond*cte1)

#        rp = Rp*kp/(k_parallel(ky) - kp)
    ################### I1 ################################################
    nmax = omegac*(int_v - np.real(alfa_p))*a/(2*np.pi)
    nmax = int(nmax)    
    print(nmax)  
    if nmax > 0:
        list_dipolos = np.linspace(0,nmax,nmax+1)
    else: 
        list_dipolos = np.linspace(0,-nmax,-nmax+1)

    factor = a*(2*np.pi)**2

    int_tot0 = 0  
    for n in list_dipolos:
        
        kx = omegac*int_v + 2*pi*n/(a*omegac)            
        exp_x = np.exp(1j*kx*x)
        
        K0 = special.kv(0,kx*y)
        K1 = special.kv(1,kx*y)
        
        term3 = (-1j*px*K0 + 2*py*K1)*kx
        
        term4 = pz*1j*np.sign(z)/y
        
        
        term_total = term3 + term4 + G1*1j*px*kx + G2*1j*py*k1  - G3*pz*k1
        
        int_tot0 = term_total + int_tot0*exp_x
    

    return int_tot0/factor           


#%%  