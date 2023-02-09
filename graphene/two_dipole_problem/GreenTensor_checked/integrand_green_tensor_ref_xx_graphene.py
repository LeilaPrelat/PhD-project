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

#%%

def green_tensor_NUM_p_fresnel_QE_integrand(omegac,epsi1,epsi2,hbmu,hbgama,R,phi,z,ze,zp,u):
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

    
    expB = np.exp(-u*z_dip_barra) 
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    
    rp_num = epsi2*1j*u - epsi1*1j*u - cte1*cond*u**2
    rp_den = epsi2*1j*u + epsi1*1j*u - cte1*cond*u**2
    rp = rp_num/rp_den

    J0 = special.jn(0,u*Rbarra) 
    J2 = special.jn(2,u*Rbarra) 
    

############ I0 5 ##################################
##################################################################
    Int06_B_function_re = np.real(rp*(J0 - np.cos(2*phi)*J2)*expB*u**2)
    Int06_B_function_im = np.imag(rp*(J0 - np.cos(2*phi)*J2)*expB*u**2)

    
    
    int0_6 = Int06_B_function_re + 1j*Int06_B_function_im
    
    cte = 0.5*k1_3
    
    return int0_6*cte 


#%%