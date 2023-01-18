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

def G1_num(omegac,epsi1,epsi2,hbmu,hbgama,int_v,zp,y,alpha_y):     
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

    
#    cota_inf = 0.01
#    cota_sup = 1001*k1
#    
#    a = 0.05
    alpha_x = int_v 
    
    alpha_parallel =  np.sqrt(alpha_x**2 + alpha_y**2)
    
    exp_electron_f = np.exp(-alpha_parallel*2*k1*zp)*np.exp(1j*alpha_y*k1*y) 



    rp_num =  epsi2*1j*alpha_parallel- epsi1*1j*alpha_parallel - cte1*cond*(alpha_parallel)**2
    rp_den = epsi2*1j*alpha_parallel + epsi1*1j*alpha_parallel - cte1*cond*(alpha_parallel)**2
    rp =  rp_num/rp_den
    
    term2 =  rp*exp_electron_f/alpha_parallel
    
    return term2  

#%%

def G1_pole_aprox(omegac,epsi1,epsi2,hbmu,hbgama,int_v,zp,y,alpha_y):     
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

#
#    a = 0.05
    alpha_x = int_v 
    alpha_parallel =  np.sqrt(alpha_x**2 + alpha_y**2)
    
    exp_electron_f =  np.exp(-alpha_parallel*2*k1*zp)*np.exp(1j*alpha_y*k1*y) 

    rp =  Rp*alfa_p/(alpha_parallel - alfa_p)   
    
    term2 = rp*exp_electron_f/alpha_parallel

    return term2 

#%%

def G2_num(omegac,epsi1,epsi2,hbmu,hbgama,int_v,zp,y,alpha_y):     
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
    
    alpha_x = int_v 
    alpha_parallel = np.sqrt(alpha_x**2 + alpha_y**2)
    
    exp_electron_f = np.exp(-alpha_parallel*2*k1*zp)*1j*np.sin(alpha_y*k1*y) 
    
    rp_num = epsi2*1j*alpha_parallel - epsi1*1j*alpha_parallel - cte1*cond*(alpha_parallel)**2
    rp_den = epsi2*1j*alpha_parallel + epsi1*1j*alpha_parallel - cte1*cond*(alpha_parallel)**2
    rp = rp_num/rp_den   
    
    term2_re_f = np.real(2*alpha_y*rp*exp_electron_f/alpha_parallel)
    term2_im_f = np.imag(2*alpha_y*rp*exp_electron_f/alpha_parallel)

    term2_int = term2_re_f + 1j*term2_im_f
    
    return term2_int  

#%%

def G2_pole_aprox(omegac,epsi1,epsi2,hbmu,hbgama,int_v,zp,y,alpha_y):     
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

    alpha_x = int_v 
    alpha_parallel = np.sqrt(alpha_x**2 + alpha_y**2)
    
    exp_electron_f = np.exp(-alpha_parallel*2*k1*zp)*1j*np.sin(alpha_y*k1*y) 

    rp = Rp*alfa_p/(alpha_parallel - alfa_p)   
    
    term2_re_f = np.real(alpha_y*2*rp*exp_electron_f/alpha_parallel)
    term2_im_f = np.imag(alpha_y*2*rp*exp_electron_f/alpha_parallel)
    

    term2_int = term2_re_f + 1j*term2_im_f

    
    return term2_int  

#%%
    
def G3_num(omegac,epsi1,epsi2,hbmu,hbgama,int_v,zp,y,alpha_y):     
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

    a = 0.05   
    alpha_x = int_v 
    alpha_parallel = np.sqrt(alpha_x**2 + alpha_y**2)
    
    exp_electron_f = np.exp(-alpha_parallel*2*k1*zp)*np.exp(1j*alpha_y*k1*y) 

    rp_num = epsi2*1j*alpha_parallel - epsi1*1j*alpha_parallel - cte1*cond*(alpha_parallel)**2
    rp_den = epsi2*1j*alpha_parallel + epsi1*1j*alpha_parallel - cte1*cond*(alpha_parallel)**2
    rp = rp_num/rp_den   
    
    term2_re_f = np.real(alpha_y*rp*exp_electron_f)
    term2_im_f = np.imag(alpha_y*rp*exp_electron_f)
    
    term2_int = term2_re_f + 1j*term2_im_f
    
    return term2_int  

#%%

def G3_num_odd(omegac,epsi1,epsi2,hbmu,hbgama,int_v,zp,y,alpha_y):     
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

    a = 0.05   
    alpha_x = int_v 
    alpha_parallel = np.sqrt(alpha_x**2 + alpha_y**2)
    
    exp_electron_f = np.exp(-alpha_parallel*2*k1*zp)*np.cos(alpha_y*k1*y) 

    rp_num = epsi2*1j*alpha_parallel - epsi1*1j*alpha_parallel - cte1*cond*(alpha_parallel)**2
    rp_den = epsi2*1j*alpha_parallel + epsi1*1j*alpha_parallel - cte1*cond*(alpha_parallel)**2
    rp = rp_num/rp_den   
    
    term2_re_f = np.real(alpha_y*rp*exp_electron_f)
    term2_im_f = np.imag(alpha_y*rp*exp_electron_f)
    
    term2_int = term2_re_f + 1j*term2_im_f
    
    return 2*term2_int  

#%%

def G3_pole_aprox(omegac,epsi1,epsi2,hbmu,hbgama,int_v,zp,y,alpha_y):     
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

        
    a = 0.05   
    alpha_x = int_v 
    alpha_parallel = np.sqrt(alpha_x**2 + alpha_y**2)
    
    exp_electron_f = np.exp(-alpha_parallel*2*k1*zp)*np.exp(1j*alpha_y*k1*y) 

    rp = Rp*alfa_p/(alpha_parallel - alfa_p)   
    
    term2_re_f =  np.real(alpha_y*rp*exp_electron_f)
    term2_im_f =  np.imag(alpha_y*rp*exp_electron_f)

    term2_int = term2_re_f + 1j*term2_im_f
    
    return term2_int  


#%%
    


def G3_pole_aprox_odd(omegac,epsi1,epsi2,hbmu,hbgama,int_v,zp,y,alpha_y):     
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

        
    a = 0.05   
    alpha_x = int_v 
    alpha_parallel = np.sqrt(alpha_x**2 + alpha_y**2)
    
    exp_electron_f = np.exp(-alpha_parallel*2*k1*zp)*np.cos(alpha_y*k1*y) 

    rp = Rp*alfa_p/(alpha_parallel - alfa_p)   
    
    term2_re_f =  np.real(alpha_y*rp*exp_electron_f)
    term2_im_f =  np.imag(alpha_y*rp*exp_electron_f)

    term2_int = term2_re_f + 1j*term2_im_f
    
    return 2*term2_int  