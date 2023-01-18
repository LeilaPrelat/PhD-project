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
path_basic = path.replace('GreenTensor/self/' + name_this_py,'')
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

def green_tensor_NUM_QE(omegac,epsi1,epsi2,hbmu,hbgama,x,y,z,xD,yD,zD,zp):
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
    Gxx self (superficie)
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
    Rbarra =  k1*np.sqrt(2)*xD
    z_dip_barra = k1*(np.abs(z) + 2*zp + np.abs(zD))
    
    
    expB = lambda u: np.exp(-u*z_dip_barra) 
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    
    rp_num = lambda u: epsi2*1j*u - epsi1*1j*u - cte1*cond*u**2
    rp_den = lambda u: epsi2*1j*u + epsi1*1j*u - cte1*cond*u**2
    rp = lambda u: rp_num(u)/rp_den(u)
    
    rs_num = lambda u: 1j*u - 1j*u - cond/cte1
    rs_den = lambda u: 1j*u + 1j*u + cond/cte1
    rs = lambda u: rs_num(u)/rs_den(u)

    J0 = lambda u: special.jn(0,u*Rbarra) 
    J2 = lambda u: special.jn(2,u*Rbarra) 
    
    
    cota_inf = 1
    cota_sup = 1001    
    
############ I0 5 ###################################
    Int05_B_function_re = lambda u: np.real(rs(u)*J0(u)*expB(u))
    Int05_B_function_im = lambda u: np.imag(rs(u)*J0(u)*expB(u))

    int05B_re,err = integrate.quad(Int05_B_function_re, cota_inf, cota_sup)
    int05B_im,err = integrate.quad(Int05_B_function_im, cota_inf, cota_sup)
############ I2 5 ############### ###################
    Int25_B_function_re = lambda u: np.real(rs(u)*J2(u)*expB(u))
    Int25_B_function_im = lambda u: np.imag(rs(u)*J2(u)*expB(u))

    int25B_re,err = integrate.quad(Int25_B_function_re, cota_inf, cota_sup)
    int25B_im,err = integrate.quad(Int25_B_function_im,cota_inf, cota_sup)
################### I0 6 ################ ###################
    Int06_B_function_re = lambda u: np.real(rp(u)*J0(u)*expB(u)*u**2)
    Int06_B_function_im = lambda u: np.imag(rp(u)*J0(u)*expB(u)*u**2)

    int06B_re,err = integrate.quad(Int06_B_function_re, cota_inf, cota_sup)
    int06B_im,err = integrate.quad(Int06_B_function_im, cota_inf, cota_sup)
################### I2 6 ################ ###################
    Int26_B_function_re = lambda u: np.real(rp(u)*J2(u)*expB(u)*u**2)
    Int26_B_function_im = lambda u: np.imag(rp(u)*J2(u)*expB(u)*u**2)

    int26B_re,err = integrate.quad(Int26_B_function_re, cota_inf, cota_sup)
    int26B_im,err = integrate.quad(Int26_B_function_im, cota_inf, cota_sup)

    int0_5 = int05B_re + 1j*int05B_im
    int2_5 = int25B_re + 1j*int25B_im
    
    int0_6 = int06B_re + 1j*int06B_im
    int2_6 = int26B_re + 1j*int26B_im
    
    cte = 0.5*k1_3
    
    return int0_5*cte, int2_5*cte, int0_6*cte, int2_6*cte


#%%


def green_tensor_NUM(omegac,epsi1,epsi2,hbmu,hbgama,x,y,z,xD,yD,zD,zp):
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
    zp : posicion del plano (>0)
    Returns
    -------
    Gxx reflejado (superficie)
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
    Rbarra =  k1*np.sqrt(2)*xD
    z_dip_barra =  k1*(np.abs(z) + 2*zp + np.abs(zD))

    aux2 = n2/n1
    alpha_z1 = lambda u: np.sqrt(1-u**2) if u<1 else 1j*np.sqrt(u**2-1)
    alpha_z2 = lambda u: np.sqrt(aux2-u**2) if u<aux2 else 1j*np.sqrt(u**2-aux2)
    
    exp = lambda u: np.exp(1j*alpha_z1(u)*z_dip_barra) 

    J0 = lambda u: special.jn(0,u*Rbarra) 
    J2 = lambda u: special.jn(2,u*Rbarra) 
    
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    
    rp_num = lambda u: epsi2*alpha_z1(u) - epsi1*alpha_z2(u) + cte1*cond*alpha_z1(u)*alpha_z2(u)
    rp_den = lambda u: epsi2*alpha_z1(u) + epsi1*alpha_z2(u) + cte1*cond*alpha_z1(u)*alpha_z2(u)
    rp = lambda u: rp_num(u)/rp_den(u)
    
    rs_num = lambda u: alpha_z1(u) - alpha_z2(u) - cond/cte1
    rs_den = lambda u: alpha_z1(u) + alpha_z2(u) + cond/cte1
    rs = lambda u: rs_num(u)/rs_den(u)
    
    cota_sup1 = 0.95
    cota_inf1 = 1.05
    cota_sup2 = 1000

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
    
    return int0_5*cte, int2_5*cte, -int0_6*cte, -int2_6*cte  

#%%


def green_tensor_NUM_PP_QE(omegac,epsi1,epsi2,hbmu,hbgama,z,xD,zD,zp):
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
    zD : coordenada z del dipolo 
    zp : posicion del plano (>0)
    Returns
    -------
    Gxx self (superficie)
    con z hacia abajo (convencion del paper)
    numerico con la aprox QE (no hay alpha_z)
    + p polarization approx (fresnel coefficient)
    solo para x,y,z = xD,yD,zD = 0,0,0
    """
    E = omegac*aux
    k0 = omegac #=omega/c
    # x_y = ky/k0
    n1 = epsi1*mu1
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1
    k1_3 = (k0*cte1)**3
    Rbarra =  k1*np.sqrt(2)*xD
    z_dip_barra =  k1*(np.abs(z) + 2*zp + np.abs(zD))
    

    
    expB = lambda u: np.exp(-u*z_dip_barra) 
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)

    alfa_p = (k0/k1)*1j*(epsi1 + epsi2)/cond
    Rp = 2*epsi1/(epsi1 + epsi2)
    
    rp = lambda u: u**2/(u-alfa_p)
    


    J0 = lambda u: special.jn(0,u*Rbarra) 
    J2 = lambda u: special.jn(2,u*Rbarra) 
    
    
    cota_inf = 1
    cota_sup = 1001    


################### I2 6 ################ ###################
    Int26_B_function_re = lambda u: np.real(rp(u)*(J2(u)+J0(u))*expB(u))
    Int26_B_function_im = lambda u: np.imag(rp(u)*(J2(u)+J0(u))*expB(u))

    int26B_re,err = integrate.quad(Int26_B_function_re, cota_inf, cota_sup)
    int26B_im,err = integrate.quad(Int26_B_function_im, cota_inf, cota_sup)

    int2_6 = int26B_re + 1j*int26B_im
    
    cte = -1j*0.5*k1_3*Rp*alfa_p
    
    return int2_6*cte


#%%


def green_tensor_ana1(omegac,epsi1,epsi2,hbmu,hbgama,z,xD,zD,zp):
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
    zD : coordenada z del dipolo 
    zp : posicion del plano (>0)
    Returns
    -------
    Gxx self (superficie)
    con z hacia abajo (convencion del paper)
    numerico con la aprox QE (no hay alpha_z)
    + p polarization approx (fresnel coefficient)
    solo para x,y,z = xD,yD,zD = 0,0,0
    Aprox analitica
    """
    E = omegac*aux
    k0 = omegac #=omega/c
    # x_y = ky/k0
    n1 = epsi1*mu1
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1
    k1_3 = (k0*cte1)**3
    Rbarra =  k1*np.sqrt(2)*xD
    z_dip_barra = k1*(np.abs(z) + 2*zp + np.abs(zD))
    

    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)

    alfa_p = (k0/k1)*1j*(epsi1 + epsi2)/cond
 
    Rp = 2*epsi1/(epsi1 + epsi2)
    
    expB = np.exp(-alfa_p*z_dip_barra) 
    
    J0 = special.jn(0,alfa_p*Rbarra) 
    J2 = special.jn(2,alfa_p*Rbarra) 
    
    
    final = np.pi*k1_3*Rp*(alfa_p**3)*(J0 + J2)
    
    return final*expB


#%%


def green_tensor_ana2(omegac,epsi1,epsi2,hbmu,hbgama,z,xD,zD,zp):
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
    zD : coordenada z del dipolo 
    zp : posicion del plano (>0)
    Returns
    -------
    Gxx self (superficie)
    con z hacia abajo (convencion del paper)
    numerico con la aprox QE (no hay alpha_z)
    + p polarization approx (fresnel coefficient)
    solo para x,y,z = xD,yD,zD = 0,0,0
    Aprox analitica
    """
    E = omegac*aux
    k0 = omegac #=omega/c
    # x_y = ky/k0
    n1 = epsi1*mu1
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1
    k1_3 = (k0*cte1)**3

    xD_tilde = k1*xD
    z_dip_barra = k1*(np.abs(z) + 2*zp + np.abs(zD))
    

    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)

    alfa_p = (k0/k1)*1j*(epsi1 + epsi2)/cond
 
    Rp = 2*epsi1/(epsi1 + epsi2)
    
    
    term0 = np.sqrt(z_dip_barra**2 + xD_tilde**2)
    term1 = -z_dip_barra + term0
    term2 = xD_tilde**2 + z_dip_barra**2
    
    rta1 = 3*(z_dip_barra*term1)**2/(term2**(2.5)) - 2*z_dip_barra*term1*(2*z_dip_barra/term0-2)/(term2**(1.5)) - term1**2/(term2**(1.5)) 

    rta2 = (2*z_dip_barra/term0-2)*(z_dip_barra/term0-1)/(term2**(0.5)) - z_dip_barra*(term1**2)/(term2**(3/2)) + term1*(2*z_dip_barra/term0-2)/(term2**(1/2))
    
    rta3 = term1*(-2*z_dip_barra**2/(term2**(3/2)) + 2/(term2**(1/2)))/(term2**(0.5))
    
    final = rta1 + rta2 + rta3
    
    term5 = alfa_p*np.sqrt(term2)
    
    aux_f = -1j*k1_3*0.5*Rp*alfa_p
    
    return aux_f*(final + term5)


#%%