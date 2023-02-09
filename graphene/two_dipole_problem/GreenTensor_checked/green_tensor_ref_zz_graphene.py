#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila

green tensor : integrales resueltas numericamente
luego de aplicar la aprox QE + sin aplicar la aprox QE
"""
from scipy import integrate
from scipy import special
import numpy as np
import sys
import os 

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('GreenTensor_checked/' + name_this_py,'')
# path_graphene =  path_basic.replace('/' + 'plane','')
#print('Importar modulos necesarios para este codigo')

try:
    sys.path.insert(1, path_basic)
    from graphene_sigma import sigma_DL
except ModuleNotFoundError:
    print('graphene_sigma.py no se encuentra en ' + path_basic)

try:
    sys.path.insert(1, path_basic)
    from constants import constantes
except ModuleNotFoundError:
    print('constants.py no se encuentra en ' + path_basic)

pi,hb,c,alfac,mu1,mu2 = constantes()
aux = c*hb

limitt = 100

#%%

def green_tensor_ref_fresnel_rp(omegac,epsi1,epsi2,hbmu,hbgama,R,phi,z,ze,zp):
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
    xD : coordenada x del dipolo 
    yD : coordenada y del dipolo 
    zD : coordenada z del dipolo 
    zp : coordenada zp del plano
    Returns
    -------
    Gxx reflejado (superficie)
    con z hacia abajo (convencion del paper)
    numerico con la aprox QE (no hay alpha_z)
    """
    E = omegac*aux
    k0 = omegac #=omega/c
    # x_y = ky/k0
    n1 = epsi1*mu1
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1
    k1_3 = (k0*cte1)**3
    
    
    Rbarra =  k1*R
    z_dip_barra = k1*np.abs(np.abs(z) + 2*zp + np.abs(ze))

    
    expB = lambda u: np.exp(-u*z_dip_barra) 
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    
    rp_num = lambda u: epsi2*1j*u - epsi1*1j*u - cte1*cond*u**2
    rp_den = lambda u: epsi2*1j*u + epsi1*1j*u - cte1*cond*u**2
    rp = lambda u: rp_num(u)/rp_den(u)

    J0 = lambda u: special.jn(0,u*Rbarra) 

    
    cota_inf = 0.01
    cota_sup = 400
    
############ I0 5 ##################################
##################################################################
    Int06_B_function_re = lambda u: np.real(rp(u)*J0(u)*expB(u)*u**2)
    Int06_B_function_im = lambda u: np.imag(rp(u)*J0(u)*expB(u)*u**2)

    int06B_re,err = integrate.quad(Int06_B_function_re, cota_inf, cota_sup)
    int06B_im,err = integrate.quad(Int06_B_function_im, cota_inf, cota_sup)
##################################################################
    
    
    int0_6 = int06B_re + 1j*int06B_im
    
    cte = k1_3

    return int0_6*cte 

def green_tensor_ref_fresnel(omegac,epsi1,epsi2,hbmu,hbgama,R,phi,z,ze,zp):
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
    xD : coordenada x del dipolo 
    yD : coordenada y del dipolo 
    zD : coordenada z del dipolo 
    zp : coordenada zp del plano
    Returns
    -------
    Gxx reflejado (superficie)
    con z hacia abajo (convencion del paper)
    numerico con la aprox QE (no hay alpha_z)
    """
    E = omegac*aux
    k0 = omegac #=omega/c
    # x_y = ky/k0
    n1 = epsi1*mu1
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1
    k1_3 = (k0*cte1)**3
    
    
    Rbarra =  k1*R
    z_dip_barra = k1*np.abs(np.abs(z) + 2*zp + np.abs(ze))

    
    expB = lambda u: np.exp(-u*z_dip_barra) 
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    
    rp_num = lambda u: epsi2*1j*u - epsi1*1j*u - cte1*cond*u**2
    rp_den = lambda u: epsi2*1j*u + epsi1*1j*u - cte1*cond*u**2
    rp = lambda u: rp_num(u)/rp_den(u)


    J0 = lambda u: special.jn(0,u*Rbarra)     
    
    cota_inf = 0.01
    cota_sup = 400
    
############ I0 5 ##################################
##################################################################
    Int06_B_function_re = lambda u: np.real(rp(u)*J0(u)*expB(u)*u**2)
    Int06_B_function_im = lambda u: np.imag(rp(u)*J0(u)*expB(u)*u**2)

    int06B_re,err = integrate.quad(Int06_B_function_re, cota_inf, cota_sup,limit = limitt)
    int06B_im,err = integrate.quad(Int06_B_function_im, cota_inf, cota_sup,limit = limitt)
##################################################################
    
    
    int0_6 = int06B_re + 1j*int06B_im
    
    cte = k1_3

    return int0_6*cte 

#%%


def green_tensor_ref_pole_aprox(omegac,epsi1,epsi2,hbmu,hbgama,R,phi,z,ze,zp):
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
    xD : coordenada x del dipolo 
    yD : coordenada y del dipolo 
    zD : coordenada z del dipolo 
    zp : coordenada zp del plano
    Returns
    -------
    Gxx reflejado (superficie)
    con z hacia abajo (convencion del paper)
    numerico con la aprox QE (no hay alpha_z)
    """
    E = omegac*aux
    k0 = omegac #=omega/c
    # x_y = ky/k0
    n1 = epsi1*mu1
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1
    k1_3 = (k0*cte1)**3
    
    
    Rbarra =  k1*R
    z_dip_barra = k1*np.abs(np.abs(z) + 2*zp + np.abs(ze))

    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    
    alpha_p = 1j*(epsi1 + epsi2)/(cond)
    Rp = 2*epsi1/(epsi1 + epsi2)
    
    expB = lambda u: np.exp(-u*z_dip_barra) 
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    

    rp = lambda u: Rp*alpha_p/(u - alpha_p)

    J0 = lambda u: special.jn(0,u*Rbarra) 
 
    
    cota_inf = 0.01/k1
    cota_sup = 600/k1
    
############ I0 5 ##################################
##################################################################
    Int06_B_function_re = lambda u: np.real(rp(u)*J0(u)*expB(u)*u**2)
    Int06_B_function_im = lambda u: np.imag(rp(u)*J0(u)*expB(u)*u**2)

    int06B_re,err = integrate.quad(Int06_B_function_re, cota_inf, cota_sup)
    int06B_im,err = integrate.quad(Int06_B_function_im, cota_inf, cota_sup)
##################################################################
    
    
    int0_6 = int06B_re + 1j*int06B_im
    
    cte = k1_3

    return int0_6*cte 

#%%

#def green_tensor_PP1(omegac,epsi1,epsi2,hbmu,hbgama,R,phi,z,ze,zp):
#    """
#    Parameters
#    ----------
#    omegac : omega/c = k0 en 1/micrometros    
#    epsi1 : epsilon del medio de arriba del plano
#    epsi2 : epsilon del medio de abajo del plano
#    hbmu : chemical potential in eV  
#    hbgama : collision frequency in eV
#    x : coordenada x 
#    y : coordenada y 
#    z : coordenada z
#    xD : coordenada x del dipolo 
#    yD : coordenada y del dipolo 
#    zp : coordenada zp del plano
#    Returns
#    -------
#    Gxx reflejado (superficie)
#    con z hacia abajo (convencion del paper)
#    analitico con la aprox QE + PP
#    """
#    E = omegac*aux
##    k0 = omegac #=omega/c
#    # x_y = ky/k0
##    n1 = epsi1*mu1
##    cte1 = np.sqrt(n1)
#    k1 = omegac
#    k1_3 = k1**3
#    Rbarra =  k1*R
#
#    
#    z_dip_barra = k1*np.abs(np.abs(z) + 2*zp + np.abs(ze))
#
#    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
#    
#    alpha_p = 1j*(epsi1 + epsi2)/(cond)
#    Rp = 2*epsi1/(epsi1 + epsi2)
#    
#    expo = np.exp(-alpha_p*z_dip_barra)
#    alphap_3 = alpha_p**3
#    
#    arg = alpha_p*Rbarra
#    
#    
#    J0 = special.jn(0,arg)
#    J2 = special.jn(2,arg)
#
#
#    final = -0.5*Rp*alphap_3*k1_3*(J0 + np.cos(2*phi)*J2)*expo
#
#    
#    return final

#%%


def green_tensor_PP2(omegac,epsi1,epsi2,hbmu,hbgama,R,phi,z,ze,zp):
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
    xD : coordenada x del dipolo 
    yD : coordenada y del dipolo 
    zp : coordenada zp del plano
    Returns
    -------
    Gxx reflejado (superficie)
    con z hacia abajo (convencion del paper)
    analitico con la aprox QE + PP
    """
    E = omegac*aux
#    k0 = omegac #=omega/c
    # x_y = ky/k0
#    n1 = epsi1*mu1
#    cte1 = np.sqrt(n1)
    k1 = omegac
    k1_3 = k1**3
    Rbarra =  k1*R
    
    z_dip_barra = k1*np.abs(np.abs(z) + 2*zp + np.abs(ze))

    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    
    alpha_p = 1j*(epsi1 + epsi2)/(cond)
    Rp = 2*epsi1/(epsi1 + epsi2)
    
    expo = np.exp(-alpha_p*z_dip_barra)
    alphap_3 = alpha_p**3
    
    arg = alpha_p*Rbarra
    
    
    H0 = special.hankel1(0,arg)


    final = np.pi*1j*Rp*alphap_3*k1_3*H0*expo
    
    return final

#%%