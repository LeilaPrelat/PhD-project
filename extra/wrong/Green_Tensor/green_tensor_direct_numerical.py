#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila

green tensor 
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
# path_graphene =  path_basic.replace('/' + 'plane','')
#print('Importar modulos necesarios para este codigo')

try:
    sys.path.insert(1, path_basic)
    from constants import constantes
except ModuleNotFoundError:
    print('constants.py no se encuentra en ' + path_basic)

pi,hb,c,alfac,mu1,mu2 = constantes()
aux = c*hb

#%%

def green_tensor_functionDIR(omegac,epsi1,x,y,z,xD,yD,zD):
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
    """

    k0 = omegac #=omega/c
    # x_y = ky/k0
    n1 = epsi1*mu1
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1
    k1_3 = (k0*cte1)**3
    Rbarra =  np.sqrt((x-xD)**2 + (y-yD)**2)/k1
    z_dip_barra = np.abs(z-zD)/k1
    phi = np.arctan2(np.abs(y-yD),np.abs(x-xD))

    alpha_z1 = lambda u: np.sqrt(1-u**2)
    alpha_z2 = lambda u: 1j*np.sqrt(u**2-1)
    
    expA = lambda u: np.exp(1j*alpha_z1(u)*z_dip_barra) 
    expB = lambda u: np.exp(1j*alpha_z2(u)*z_dip_barra)

    J0 = lambda u: special.jn(0,u*Rbarra) 
    J2 = lambda u: special.jn(2,u*Rbarra) 
############ I0 5 ###################################
    Int05_A_function_re = lambda u: np.real(J0(u)*expA(u)/alpha_z1(u))
    Int05_B_function_re = lambda u: np.real(J0(u)*expB(u)/alpha_z2(u))
    Int05_A_function_im = lambda u: np.imag(J0(u)*expA(u)/alpha_z1(u))
    Int05_B_function_im = lambda u: np.imag(J0(u)*expB(u)/alpha_z2(u))

    int05A_re,err = integrate.quad(Int05_A_function_re, 0, 1)
    int05B_re,err = integrate.quad(Int05_B_function_re, 1, np.inf)
    int05A_im,err = integrate.quad(Int05_A_function_im, 0, 1)
    int05B_im,err = integrate.quad(Int05_B_function_im, 1, np.inf)

############ I2 5 ############### ###################
    Int25_A_function_re = lambda u: np.real(np.cos(2*phi)*J2(u)*expA(u)/alpha_z1(u))
    Int25_B_function_re = lambda u: np.real(np.cos(2*phi)*J2(u)*expB(u)/alpha_z2(u))
    Int25_A_function_im = lambda u: np.imag(np.cos(2*phi)*J2(u)*expA(u)/alpha_z1(u))
    Int25_B_function_im = lambda u: np.imag(np.cos(2*phi)*J2(u)*expB(u)/alpha_z2(u))

    int25A_re,err = integrate.quad(Int25_A_function_re, 0, 1)
    int25B_re,err = integrate.quad(Int25_B_function_re, 1, np.inf)
    int25A_im,err = integrate.quad(Int25_A_function_im, 0, 1)
    int25B_im,err = integrate.quad(Int25_B_function_im, 1, np.inf)
################### I0 6 ################ ###################
    Int06_A_function_re = lambda u: np.real(J0(u)*expA(u)*alpha_z1(u))
    Int06_B_function_re = lambda u: np.real(J0(u)*expB(u)*alpha_z2(u))
    Int06_A_function_im = lambda u: np.imag(J0(u)*expA(u)*alpha_z1(u))
    Int06_B_function_im = lambda u: np.imag(J0(u)*expB(u)*alpha_z2(u))

    int06A_re,err = integrate.quad(Int06_A_function_re, 0, 1)
    int06B_re,err = integrate.quad(Int06_B_function_re, 1, np.inf)
    int06A_im,err = integrate.quad(Int06_A_function_im, 0, 1)
    int06B_im,err = integrate.quad(Int06_B_function_im, 1, np.inf)

################### I2 6 ################ ###################
    Int26_A_function_re = lambda u: np.real(np.cos(2*phi)*J2(u)*expA(u)*alpha_z1(u))
    Int26_B_function_re = lambda u: np.real(np.cos(2*phi)*J2(u)*expB(u)*alpha_z2(u))
    Int26_A_function_im = lambda u: np.imag(np.cos(2*phi)*J2(u)*expA(u)*alpha_z1(u))
    Int26_B_function_im = lambda u: np.imag(np.cos(2*phi)*J2(u)*expB(u)*alpha_z2(u))

    int26A_re,err = integrate.quad(Int26_A_function_re, 0, 1)
    int26B_re,err = integrate.quad(Int26_B_function_re, 1, np.inf)
    int26A_im,err = integrate.quad(Int26_A_function_im, 0, 1)
    int26B_im,err = integrate.quad(Int26_B_function_im, 1, np.inf)

    int0_5 = int05A_re + int05B_re + 1j*int05A_im + 1j*int05B_im
    int2_5 = int25A_re + int25B_re + 1j*int25A_im + 1j*int25B_im
    
    int0_6 = int06A_re + int06B_re + 1j*int06A_im + 1j*int06B_im
    int2_6 = int26A_re + int26B_re + 1j*int26A_im + 1j*int26B_im
    
    cte = 0.5*1j*k1_3
    
    return int0_5*cte, int2_5*cte, int0_6*cte, int2_6*cte  

#%%
