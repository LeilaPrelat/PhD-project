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
#path_constants =  path_basic.replace('GreenTensor/' + 'reflected','')
# path_graphene =  path_basic.replace('/' + 'plane','')
#print('Importar modulos necesarios para este codigo')

try:
    sys.path.insert(1, path_basic)
    from Ag_sigma import sigma_DL
except ModuleNotFoundError:
    print('Ag_sigma.py no se encuentra en ' + path_basic)

try:
    sys.path.insert(1, path_basic)
    from constants import constantes
except ModuleNotFoundError:
    print('constants.py no se encuentra en ' + path_basic)

pi,hb,c,alfac,mu1,mu2 = constantes()
aux = c*hb

#%%

def green_tensor_NUM_QE(omegac,epsi1,epsi2,omega_bulk,hbar_gamma_in,d,epsilon_b,x,y,z,xD,yD,zD,zp):
    """
    Parameters
    ----------
    omegac : omega/c = k0 en 1/micrometros    
    epsi1 : epsilon del medio de arriba del plano
    epsi2 : epsilon del medio de abajo del plano
    omega_bulk : chemical potential in eV  
    hbar_gamma_in : collision frequency in eV
    d : espesor del plano de Ag
    epsilon_b : permeabilidad electrica del plano
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
    k0 = omegac #=omega/c
    omega = k0*c
    # x_y = ky/k0
    n1 = epsi1*mu1
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1
    k1_3 = (k0*cte1)**3
    Rbarra =  k1*np.sqrt((x + xD)**2 + (y + yD)**2)
    z_dip_barra = k1*np.abs(z + 2*np.abs(zD))
    phi = np.arctan2(np.abs(y + yD),np.abs(x + xD))
    
    expB = lambda u: np.exp(-u*z_dip_barra) 
    cond = 4*np.pi*sigma_DL(omega,omega_bulk,hbar_gamma_in,d,epsilon_b)/c

    rp_num = lambda u: epsi2*1j*u - epsi1*1j*u - cte1*cond*u**2
    rp_den = lambda u: epsi2*1j*u + epsi1*1j*u - cte1*cond*u**2
    rp = lambda u: rp_num(u)/rp_den(u)
    
    rs_num = lambda u: 1j*u - 1j*u - cond/cte1
    rs_den = lambda u: 1j*u + 1j*u + cond/cte1
    rs = lambda u: rs_num(u)/rs_den(u)

    J0 = lambda u: special.jn(0,u*Rbarra) 
    J2 = lambda u: special.jn(2,u*Rbarra) 
############ I0 5 ###################################
    Int05_B_function_re = lambda u: np.real(rs(u)*J0(u)*expB(u))
    Int05_B_function_im = lambda u: np.imag(rs(u)*J0(u)*expB(u))

    int05B_re,err = integrate.quad(Int05_B_function_re, 0, np.inf)
    int05B_im,err = integrate.quad(Int05_B_function_im, 0, np.inf)
############ I2 5 ############### ###################
    Int25_B_function_re = lambda u: np.real(rs(u)*J2(u)*expB(u))
    Int25_B_function_im = lambda u: np.imag(rs(u)*J2(u)*expB(u))

    int25B_re,err = integrate.quad(Int25_B_function_re, 0, np.inf)
    int25B_im,err = integrate.quad(Int25_B_function_im, 0, np.inf)
################### I0 6 ################ ###################
    Int06_B_function_re = lambda u: np.real(rp(u)*J0(u)*expB(u)*u**2)
    Int06_B_function_im = lambda u: np.imag(rp(u)*J0(u)*expB(u)*u**2)

    int06B_re,err = integrate.quad(Int06_B_function_re, 0, np.inf)
    int06B_im,err = integrate.quad(Int06_B_function_im, 0, np.inf)
################### I2 6 ################ ###################
    Int26_B_function_re = lambda u: np.real(rp(u)*J2(u)*expB(u)*u**2)
    Int26_B_function_im = lambda u: np.imag(rp(u)*J2(u)*expB(u)*u**2)

    int26B_re,err = integrate.quad(Int26_B_function_re, 0, np.inf)
    int26B_im,err = integrate.quad(Int26_B_function_im, 0, np.inf)

    int0_5 = int05B_re + 1j*int05B_im
    int2_5 = int25B_re + 1j*int25B_im
    
    int0_6 = int06B_re + 1j*int06B_im
    int2_6 = int26B_re + 1j*int26B_im
    
    cte = 0.5*k1_3
    cte2 = np.cos(2*phi)  
    
    return int0_5*cte, int2_5*cte*cte2, int0_6*cte, int2_6*cte*cte2

#%%


def green_tensor_NUM(omegac,epsi1,epsi2,omega_bulk,hbar_gamma_in,d,epsilon_b,x,y,z,xD,yD,zD,zp):
    """
    Parameters
    ----------
    omegac : omega/c = k0 en 1/micrometros    
    epsi1 : epsilon del medio de arriba del plano
    epsi2 : epsilon del medio de abajo del plano
    omega_bulk : chemical potential in eV  
    hbar_gamma_in : collision frequency in eV
    d : espesor del plano de Ag
    epsilon_b : permeabilidad electrica del plano
    x : coordenada x 
    y : coordenada y 
    z : coordenada z
    xD : coordenada x del dipolo 
    yD : coordenada y del dipolo 
    zD : coordenada z del dipolo 
    zp : posicion del plano (>0)
    Returns
    -------
    Gxx reflejado (superficie)
    con z hacia abajo (convencion del paper)
    sin QE
    """
    k0 = omegac #=omega/c
    omega = omegac*c
    # x_y = ky/k0
    n1 = epsi1*mu1
    n2 = epsi2*mu2
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1
    k1_3 = (k0*cte1)**3
    Rbarra =  k1*np.sqrt((x + xD)**2 + (y + yD)**2)
    z_dip_barra = k1*np.abs(z + 2*np.abs(zD))
    phi = np.arctan2(np.abs(y + yD),np.abs(x + xD))

    aux2 = n2/n1
    alpha_z1 = lambda u: np.sqrt(1-u**2) if u<1 else 1j*np.sqrt(u**2-1)
    alpha_z2 = lambda u: np.sqrt(aux2-u**2) if u<aux2 else 1j*np.sqrt(u**2-aux2)
    
    exp = lambda u: np.exp(1j*alpha_z1(u)*z_dip_barra) 

    J0 = lambda u: special.jn(0,u*Rbarra) 
    J2 = lambda u: special.jn(2,u*Rbarra) 
    
    cond = 4*np.pi*sigma_DL(omega,omega_bulk,hbar_gamma_in,d,epsilon_b)/c
    
    rp_num = lambda u: epsi2*alpha_z1(u) - epsi1*alpha_z2(u) + cte1*cond*alpha_z1(u)*alpha_z2(u)
    rp_den = lambda u: epsi2*alpha_z1(u) + epsi1*alpha_z2(u) + cte1*cond*alpha_z1(u)*alpha_z2(u)
    rp = lambda u: rp_num(u)/rp_den(u)
    
    rs_num = lambda u: alpha_z1(u) - alpha_z2(u) - cond/cte1
    rs_den = lambda u: alpha_z1(u) + alpha_z2(u) + cond/cte1
    rs = lambda u: rs_num(u)/rs_den(u)
    
    cota_sup1 = 0.95
    cota_inf1 = 1.05
    cota_sup2 = 1e5

############ I0 5 ######################################################
    Int05_function_re = lambda u: np.real(J0(u)*rs(u)*exp(u)*u/alpha_z1(u))
    Int05_function_im = lambda u: np.imag(J0(u)*rs(u)*exp(u)*u/alpha_z1(u))
    
    int05A_re,err = integrate.quad(Int05_function_re, 0, cota_sup1)
    int05B_re,err = integrate.quad(Int05_function_re, cota_inf1, cota_sup2)
    int05A_im,err = integrate.quad(Int05_function_im, 0, cota_sup1)
    int05B_im,err = integrate.quad(Int05_function_im, cota_inf1, cota_sup2)

############ I2 5 ############### ######################################
    Int25_function_re = lambda u: np.real(J2(u)*rs(u)*exp(u)*u/alpha_z1(u))
    Int25_function_im = lambda u: np.imag(J2(u)*rs(u)*exp(u)*u/alpha_z2(u))

    int25A_re,err = integrate.quad(Int25_function_re, 0, cota_sup1)
    int25B_re,err = integrate.quad(Int25_function_re, cota_inf1, cota_sup2)
    int25A_im,err = integrate.quad(Int25_function_im, 0, cota_sup1)
    int25B_im,err = integrate.quad(Int25_function_im, cota_inf1, cota_sup2)
################### I0 6 ################ ############################
    Int06_function_re = lambda u: np.real(J0(u)*rp(u)*exp(u)*u*alpha_z1(u))
    Int06_function_im = lambda u: np.imag(J0(u)*rp(u)*exp(u)*u*alpha_z2(u))

    int06A_re,err = integrate.quad(Int06_function_re, 0, cota_sup1)
    int06B_re,err = integrate.quad(Int06_function_re, cota_inf1, cota_sup2)
    int06A_im,err = integrate.quad(Int06_function_im, 0, cota_sup1)
    int06B_im,err = integrate.quad(Int06_function_im, cota_inf1, cota_sup2)

################### I2 6 ################ ###########################
    Int26_function_re = lambda u: np.real(J2(u)*rp(u)*exp(u)*u*alpha_z1(u))
    Int26_function_im = lambda u: np.imag(J2(u)*rp(u)*exp(u)*u*alpha_z2(u))

    int26A_re,err = integrate.quad(Int26_function_re, 0, cota_sup1)
    int26B_re,err = integrate.quad(Int26_function_re, cota_inf1, cota_sup2)
    int26A_im,err = integrate.quad(Int26_function_im, 0, cota_sup1)
    int26B_im,err = integrate.quad(Int26_function_im, cota_inf1, cota_sup2)

    int0_5 = int05A_re + int05B_re + 1j*int05A_im + 1j*int05B_im
    int2_5 = int25A_re + int25B_re + 1j*int25A_im + 1j*int25B_im
    
    int0_6 = int06A_re + int06B_re + 1j*int06A_im + 1j*int06B_im
    int2_6 = int26A_re + int26B_re + 1j*int26A_im + 1j*int26B_im
    
    cte = 0.5*1j*k1_3
    cte2 = np.cos(2*phi)  
    
    return int0_5*cte, int2_5*cte*cte2, -int0_6*cte, -int2_6*cte*cte2  

#%%
