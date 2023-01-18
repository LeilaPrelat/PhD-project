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
path_basic = path.replace('GreenTensor/reflected/' + name_this_py,'')
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

#%%

def green_tensor_NUM_QE(omegac,epsi1,epsi2,hbmu,hbgama,x,y,z,xD,yD,zD,b,zp):
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
    b: electron position en z, b < 0 (dipolo en 0 y eje de z apunta hacia abajo)
    zp : coordenada zp del plano
    Returns
    -------
    Gzz reflejado (superficie)
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
    Rbarra =  k1*np.sqrt((x + xD)**2 + (y + yD)**2)
    z_dip_barra = k1*(z + np.abs(b) + 2*np.abs(zp))  
    
    expB = lambda u: np.exp(-u*z_dip_barra) 
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    
    rp_num = lambda u: epsi2*1j*u - epsi1*1j*u - cte1*cond*u**2
    rp_den = lambda u: epsi2*1j*u + epsi1*1j*u - cte1*cond*u**2
    rp = lambda u: rp_num(u)/rp_den(u)

    J0 = lambda u: special.jn(0,u*Rbarra) 
################### I0 6 ################ ###################
    Int06_B_function_re = lambda u: np.real(rp(u)*J0(u)*expB(u)*u**2)
    Int06_B_function_im = lambda u: np.imag(rp(u)*J0(u)*expB(u)*u**2)

    int06B_re,err = integrate.quad(Int06_B_function_re, 0, np.inf)
    int06B_im,err = integrate.quad(Int06_B_function_im, 0, np.inf)
################### I2 6 ################ ###################
    
    int0_6 = int06B_re + 1j*int06B_im
    
    return int0_6*k1_3

#%%


def green_tensor_NUM(omegac,epsi1,epsi2,hbmu,hbgama,x,y,z,xD,yD,zD,b,zp):
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
    b: electron position en z, b < 0 (dipolo en 0 y eje de z apunta hacia abajo)
    zp : posicion del plano (>0)
    Returns
    -------
    Gzz reflejado (superficie)
    con z hacia abajo (convencion del paper)
    sin QE
    """
    E = omegac*aux
    k0 = omegac #=omega/c
    # x_y = ky/k0
    n1 = epsi1*mu1
    n2 = epsi2*mu2
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1
    k1_3 = (k0*cte1)**3
    Rbarra =  k1*np.sqrt((x + xD)**2 + (y + yD)**2)
    z_dip_barra = k1*(z + np.abs(b) + 2*np.abs(zp))  

    aux2 = n2/n1
    alpha_z1 = lambda u: np.sqrt(1-u**2) if u<1 else 1j*np.sqrt(u**2-1)
    alpha_z2 = lambda u: np.sqrt(aux2-u**2) if u<aux2 else 1j*np.sqrt(u**2-aux2)
    
    exp = lambda u: np.exp(1j*alpha_z1(u)*z_dip_barra) 

    J0 = lambda u: special.jn(0,u*Rbarra) 
    
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    
    rp_num = lambda u: epsi2*alpha_z1(u) - epsi1*alpha_z2(u) + cte1*cond*alpha_z1(u)*alpha_z2(u)
    rp_den = lambda u: epsi2*alpha_z1(u) + epsi1*alpha_z2(u) + cte1*cond*alpha_z1(u)*alpha_z2(u)
    rp = lambda u: rp_num(u)/rp_den(u)
    
    cota_sup1 = 0.95
    cota_inf1 = 1.05
    cota_sup2 = 1e5

    Int06_function_re = lambda u: np.real(J0(u)*rp(u)*exp(u)*u**3/alpha_z1(u))
    Int06_function_im = lambda u: np.imag(J0(u)*rp(u)*exp(u)*u**3/alpha_z2(u))

    int06A_re,err = integrate.quad(Int06_function_re, 0, cota_sup1)
    int06B_re,err = integrate.quad(Int06_function_re, cota_inf1, cota_sup2)
    int06A_im,err = integrate.quad(Int06_function_im, 0, cota_sup1)
    int06B_im,err = integrate.quad(Int06_function_im, cota_inf1, cota_sup2)
    
    int0_6 = int06A_re + int06B_re + 1j*int06A_im + 1j*int06B_im
    
    cte = 1j*k1_3
    
    return int0_6*cte

#%%