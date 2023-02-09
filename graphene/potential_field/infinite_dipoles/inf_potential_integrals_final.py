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


def F1_ana(omegac,epsi1,epsi2,hbmu,hbgama,int_v,zp,a,n,y,z):     
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
    

    kx = omegac*int_v + 2*np.pi*n/a
    
    term_ky = np.sqrt(kp**2 - kx**2)
    
    kp_pv2 = np.sqrt(kp**2)
    
    term = (1 + kp/kp_pv2)

    
    exp_z = np.exp(-kp_pv2*(2*zp-z))*np.exp(1j*term_ky*np.abs(y)) 

    term1 = np.pi*1j*Rp*kp*exp_z*term/term_ky

###################################################################################################

    return term1 


#%%

def F1_num(omegac,epsi1,epsi2,hbmu,hbgama,int_v,zp,a,n,y,z):     
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

#    k1_2 = (k0*cte1)**2
 #   n_v1 = int_v/cte1
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)  

    
    cota_inf = 1*omegac
    cota_sup = 600*omegac
    

    alpha_x = int_v + 2*np.pi*n/(omegac*a)
    
    alpha_parallel = lambda alpha_y: np.sqrt(alpha_x**2 + alpha_y**2)
    
    exp_z = lambda alpha_y: np.exp(-alpha_parallel(alpha_y)*omegac*(2*zp-z))

    exp_y = lambda alpha_y: np.exp(1j*omegac*alpha_y*np.abs(y))

    rp_num = lambda alpha_y: epsi2*1j*alpha_parallel(alpha_y) - epsi1*1j*alpha_parallel(alpha_y) - cte1*cond*(alpha_parallel(alpha_y))**2
    rp_den = lambda alpha_y: epsi2*1j*alpha_parallel(alpha_y) + epsi1*1j*alpha_parallel(alpha_y) - cte1*cond*(alpha_parallel(alpha_y))**2
    rp = lambda alpha_y: rp_num(alpha_y)/rp_den(alpha_y)   
    
    termF1_re = lambda alpha_y: np.real(rp(alpha_y)*exp_z(alpha_y)*exp_y(alpha_y)/alpha_parallel(alpha_y))
    termF1_im = lambda alpha_y: np.imag(rp(alpha_y)*exp_z(alpha_y)*exp_y(alpha_y)/alpha_parallel(alpha_y))
    
    INT_re,err = integrate.quad(termF1_re, cota_inf, cota_sup) 
    INT_im,err = integrate.quad(termF1_im, cota_inf, cota_sup) 
    
    INT = INT_re + 1j*INT_im
    
    return INT 

#%%

def F1_pole_aprox(omegac,epsi1,epsi2,hbmu,hbgama,int_v,zp,a,n,y,z):     
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

    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    Rp = 2*epsi1/(epsi1 + epsi2)
    alfa_p = 1j*(epsi1 + epsi2)/(cond*cte1)
    
    cota_inf = 1*omegac
    cota_sup = 600*omegac
    

    alpha_x = int_v + 2*np.pi*n/(omegac*a)
    
    alpha_parallel = lambda alpha_y: np.sqrt(alpha_x**2 + alpha_y**2)
    
    exp_z = lambda alpha_y: np.exp(-alpha_parallel(alpha_y)*omegac*(2*zp-z))

    exp_y = lambda alpha_y: np.exp(1j*omegac*alpha_y*np.abs(y))

    rp = lambda alpha_y: Rp*alfa_p/(alpha_parallel(alpha_y) - alfa_p)    
    
    termF1_re = lambda alpha_y: np.real(rp(alpha_y)*exp_z(alpha_y)*exp_y(alpha_y)/alpha_parallel(alpha_y))
    termF1_im = lambda alpha_y: np.imag(rp(alpha_y)*exp_z(alpha_y)*exp_y(alpha_y)/alpha_parallel(alpha_y))
    
    INT_re,err = integrate.quad(termF1_re, cota_inf, cota_sup) 
    INT_im,err = integrate.quad(termF1_im, cota_inf, cota_sup) 
    
    INT = INT_re + 1j*INT_im
    
    return INT 


#%%
    
def F2_ana(omegac,epsi1,epsi2,hbmu,hbgama,int_v,zp,a,n,y,z):     
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
    

    kx = omegac*int_v + 2*np.pi*n/a
    
    term_ky = np.sqrt(kp**2 - kx**2)

    kp_pv2 = np.sqrt(kp**2)
    
    term = (1 + kp/kp_pv2)    
    
    exp_z = np.exp(-kp_pv2*(2*zp-z))*np.exp(1j*term_ky*np.abs(y)) 

    term1 = np.pi*1j*Rp*kp*exp_z*term

###################################################################################################

    return term1 


#%%

def F2_num(omegac,epsi1,epsi2,hbmu,hbgama,int_v,zp,a,n,y,z):     
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

#    k1_2 = (k0*cte1)**2
 #   n_v1 = int_v/cte1
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)  

    
    cota_inf = 1*omegac
    cota_sup = 600*omegac
    

    alpha_x = int_v + 2*np.pi*n/(omegac*a)
    
    alpha_parallel = lambda alpha_y: np.sqrt(alpha_x**2 + alpha_y**2)
    
    exp_z = lambda alpha_y: np.exp(-alpha_parallel(alpha_y)*omegac*(2*zp-z))

    exp_y = lambda alpha_y: np.exp(1j*omegac*alpha_y*np.abs(y))

    rp_num = lambda alpha_y: epsi2*1j*alpha_parallel(alpha_y) - epsi1*1j*alpha_parallel(alpha_y) - cte1*cond*(alpha_parallel(alpha_y))**2
    rp_den = lambda alpha_y: epsi2*1j*alpha_parallel(alpha_y) + epsi1*1j*alpha_parallel(alpha_y) - cte1*cond*(alpha_parallel(alpha_y))**2
    rp = lambda alpha_y: rp_num(alpha_y)/rp_den(alpha_y)   
    
    termF1_re = lambda alpha_y: np.real(alpha_y*rp(alpha_y)*exp_z(alpha_y)*exp_y(alpha_y)/alpha_parallel(alpha_y))
    termF1_im = lambda alpha_y: np.imag(alpha_y*rp(alpha_y)*exp_z(alpha_y)*exp_y(alpha_y)/alpha_parallel(alpha_y))
    
    INT_re,err = integrate.quad(termF1_re, cota_inf, cota_sup) 
    INT_im,err = integrate.quad(termF1_im, cota_inf, cota_sup) 
    
    INT = INT_re + 1j*INT_im
    
    return INT*omegac 

#%%

def F2_pole_aprox(omegac,epsi1,epsi2,hbmu,hbgama,int_v,zp,a,n,y,z):     
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

    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    Rp = 2*epsi1/(epsi1 + epsi2)
    alfa_p = 1j*(epsi1 + epsi2)/(cond*cte1)
    
    cota_inf = 1*omegac
    cota_sup = 600*omegac

    alpha_x = int_v + 2*np.pi*n/(omegac*a)
    
    alpha_parallel = lambda alpha_y: np.sqrt(alpha_x**2 + alpha_y**2)
    
    exp_z = lambda alpha_y: np.exp(-alpha_parallel(alpha_y)*omegac*(2*zp-z))

    exp_y = lambda alpha_y: np.exp(1j*omegac*alpha_y*np.abs(y))


    rp = lambda alpha_y: Rp*alfa_p/(alpha_parallel(alpha_y) - alfa_p)    
    
    termF1_re = lambda alpha_y: np.real(alpha_y*rp(alpha_y)*exp_z(alpha_y)*exp_y(alpha_y)/alpha_parallel(alpha_y))
    termF1_im = lambda alpha_y: np.imag(alpha_y*rp(alpha_y)*exp_z(alpha_y)*exp_y(alpha_y)/alpha_parallel(alpha_y))
    
    INT_re,err = integrate.quad(termF1_re, cota_inf, cota_sup) 
    INT_im,err = integrate.quad(termF1_im, cota_inf, cota_sup) 
    
    INT = INT_re + 1j*INT_im
    
    return INT*omegac 



#%%
    

def F3_ana(omegac,epsi1,epsi2,hbmu,hbgama,int_v,zp,a,n,y,z):     
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
    

    kx = omegac*int_v + 2*np.pi*n/a
    
    term_ky = np.sqrt(kp**2 - kx**2)
    
    kp_pv2 = np.sqrt(kp**2)
    
    term_kp = (kp + kp_pv2)    
   
    exp_z = np.exp(-kp_pv2*(2*zp-z))*np.exp(1j*term_ky*np.abs(y)) 

    term1 = np.pi*1j*Rp*kp*term_kp*exp_z/term_ky

###################################################################################################

    return term1 


#%%

def F3_num(omegac,epsi1,epsi2,hbmu,hbgama,int_v,zp,a,n,y,z):     
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

#    k1_2 = (k0*cte1)**2
 #   n_v1 = int_v/cte1
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)  

    
    cota_inf = 1*omegac
    cota_sup = 600*omegac
    

    alpha_x = int_v + 2*np.pi*n/(omegac*a)
    
    alpha_parallel = lambda alpha_y: np.sqrt(alpha_x**2 + alpha_y**2)
    
    exp_z = lambda alpha_y: np.exp(-alpha_parallel(alpha_y)*omegac*(2*zp-z))

    exp_y = lambda alpha_y: np.exp(1j*omegac*alpha_y*np.abs(y))

    rp_num = lambda alpha_y: epsi2*1j*alpha_parallel(alpha_y) - epsi1*1j*alpha_parallel(alpha_y) - cte1*cond*(alpha_parallel(alpha_y))**2
    rp_den = lambda alpha_y: epsi2*1j*alpha_parallel(alpha_y) + epsi1*1j*alpha_parallel(alpha_y) - cte1*cond*(alpha_parallel(alpha_y))**2
    rp = lambda alpha_y: rp_num(alpha_y)/rp_den(alpha_y)   
    
    termF1_re = lambda alpha_y: np.real(rp(alpha_y)*exp_z(alpha_y)*exp_y(alpha_y))
    termF1_im = lambda alpha_y: np.imag(rp(alpha_y)*exp_z(alpha_y)*exp_y(alpha_y))
    
    INT_re,err = integrate.quad(termF1_re, cota_inf, cota_sup) 
    INT_im,err = integrate.quad(termF1_im, cota_inf, cota_sup) 
    
    INT = INT_re + 1j*INT_im
    
    return INT*omegac

#%%

def F3_pole_aprox(omegac,epsi1,epsi2,hbmu,hbgama,int_v,zp,a,n,y,z):     
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

    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    Rp = 2*epsi1/(epsi1 + epsi2)
    alfa_p = 1j*(epsi1 + epsi2)/(cond*cte1)
    
    cota_inf = 1*omegac
    cota_sup = 600*omegac

    alpha_x = int_v + 2*np.pi*n/(omegac*a)
    
    alpha_parallel = lambda alpha_y: np.sqrt(alpha_x**2 + alpha_y**2)
    
    exp_z = lambda alpha_y: np.exp(-alpha_parallel(alpha_y)*omegac*(2*zp-z))

    exp_y = lambda alpha_y: np.exp(1j*omegac*alpha_y*np.abs(y))


    rp = lambda alpha_y: Rp*alfa_p/(alpha_parallel(alpha_y) - alfa_p)    
    
    termF1_re = lambda alpha_y: np.real(rp(alpha_y)*exp_z(alpha_y)*exp_y(alpha_y))
    termF1_im = lambda alpha_y: np.imag(rp(alpha_y)*exp_z(alpha_y)*exp_y(alpha_y))
    
    INT_re,err = integrate.quad(termF1_re, cota_inf, cota_sup) 
    INT_im,err = integrate.quad(termF1_im, cota_inf, cota_sup) 
    
    INT = INT_re + 1j*INT_im
    
    return INT*omegac



