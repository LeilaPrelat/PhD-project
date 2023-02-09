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
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/GreenTensor_checked','')
#print('Importar modulos necesarios para este codigo')

try:
    sys.path.insert(1, path_constants)
    from constants import constantes
except ModuleNotFoundError:
    print('constants.py no se encuentra en ' + path_constants)

pi,hb,c,alfac,mu1,mu2 = constantes()
aux = c*hb

#%%

def green_tensor_NUM_QE(omegac,epsi1,R,z,ze):
    """
    Parameters
    ----------
    omegac : omega/c = k0 en 1/micrometros    
    epsi1 : epsilon del medio de arriba del plano
    x : coordenada x 
    y : coordenada y 
    z : coordenada z
    xD : coordenada x del dipolo 
    yD : coordenada y del dipolo 
    zD : coordenada z del dipolo 
    Returns
    -------
    Gxx direct (self interaction of the dipole)
    con z hacia abajo (convencion del paper)
    numerico con la aprox QE (no hay alpha_z)
    """

    k0 = omegac #=omega/c
    # x_y = ky/k0
    n1 = epsi1*mu1
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1
    k1_3 = omegac**3
#    Rbarra =  k1*np.sqrt((x-xe)**2 + (y-ye)**2)
    z_dip_barra = k1*np.abs(z-ze)
#    phi = np.arctan2(np.abs(y-ye),np.abs(x-xe))
    Rbarra = k1*R
    expB = lambda u: np.exp(-u*z_dip_barra) 

    J0 = lambda u: special.jn(0,u*Rbarra) 

    cota_inf = 0.001/omegac
    cota_sup = 600/omegac
    
############ I0 5 ###################################
    Int05_B_function_re = lambda u: np.real(u**2*J0(u)*expB(u))
    Int05_B_function_im = lambda u: np.imag(u**2*J0(u)*expB(u))

    int05B_re,err = integrate.quad(Int05_B_function_re, cota_inf, cota_sup)
    int05B_im,err = integrate.quad(Int05_B_function_im, cota_inf, cota_sup)

############ I2 5 ############### ###################
#    Int25_B_function_re = lambda u: np.real(np.cos(2*phi)*J2(u)*expB(u))
#    Int25_B_function_im = lambda u: np.imag(np.cos(2*phi)*J2(u)*expB(u))
#
#    int25B_re,err = integrate.quad(Int25_B_function_re, cota_inf, cota_sup)
#    int25B_im,err = integrate.quad(Int25_B_function_im,cota_inf, cota_sup)
#################### I0 6 ################ ###################
#    Int06_B_function_re = lambda u: np.real(J0(u)*expB(u)*u**2)
#    Int06_B_function_im = lambda u: np.imag(J0(u)*expB(u)*u**2)
#
#    int06B_re,err = integrate.quad(Int06_B_function_re,cota_inf, cota_sup)
#    int06B_im,err = integrate.quad(Int06_B_function_im, cota_inf, cota_sup)
#
#################### I2 6 ################ ###################
#    Int26_B_function_re = lambda u: np.real(np.cos(2*phi)*J2(u)*expB(u)*u**2)
#    Int26_B_function_im = lambda u: np.imag(np.cos(2*phi)*J2(u)*expB(u)*u**2)
#
#    int26B_re,err = integrate.quad(Int26_B_function_re, cota_inf, cota_sup)
#    int26B_im,err = integrate.quad(Int26_B_function_im, cota_inf, cota_sup)

    int0_5 = int05B_re + 1j*int05B_im
#    int2_5 = int25B_re + 1j*int25B_im
#    
#    int0_6 = int06B_re + 1j*int06B_im
#    int2_6 = int26B_re + 1j*int26B_im
        
    return int0_5*k1_3

#%%



def green_tensor_ANA_QE(omegac,epsi1,R,z,ze):
    """
    Parameters
    ----------
    omegac : omega/c = k0 en 1/micrometros    
    epsi1 : epsilon del medio de arriba del plano
    x : coordenada x 
    y : coordenada y 
    z : coordenada z
    xD : coordenada x del dipolo 
    yD : coordenada y del dipolo 
    zD : coordenada z del dipolo 
    Returns
    -------
    Gxx direct (self interaction of the dipole)
    con z hacia abajo (convencion del paper)
    analitico luego de aplicar QE
    """
            
    den = np.abs(z-ze)**2 + R**2
    
    return 3*np.abs(z-ze)**2*(den**(-5/2))  - den**(-3/2)

#%%
