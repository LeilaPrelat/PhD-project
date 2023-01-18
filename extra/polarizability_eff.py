#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila

polarizability alpha_E
efectiva : A/((1/alfa) -green tensor)
"""
import sys
import os 

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
#print('Importar modulos necesarios para este codigo')

try:
    sys.path.insert(1, path_basic)
    from green_tensor import green_tensor_function
except ModuleNotFoundError:
    print('green_tensor.py no se encuentra en ' + path_basic)

try:
    sys.path.insert(1, path_basic)
    from constants import constantes
except ModuleNotFoundError:
    print('constants_plane.py no se encuentra en ' + path_basic)

pi,hb,c,alfac,mu1,mu2 = constantes()
aux = c*hb

#%%

def alpha_function_eff(omegac,hbmu,z_plane,epsi1,epsi2,epsi_h,hbargama_sigma,A,hbargama_alfa,omega0,px,py,pz):
    """
    Parameters
    ----------
    omegac : omega/c = k0 en 1/micrometros    
    hbmu: potencial quimico del grafeno en eV
    z_plane: radio del cilindro en micrometros
    epsi1 : epsilon del medio de arriba del plano
    epsi2 : epsilon del medio de abajo del plano
    hbargama_sigma : collision frequency in eV of sigma (grafeno)
    A : constante A de polarizability 
    hbargama_alfa : collision frequency in eV of the polarizability
    omega0 : frequency of transition of the polarizability
    px : coordenada x del dipolo 
    py : coordenada y del dipolo 
    pz : coordenada z del dipolo     
    Returns
    -------
    polarizabilty efectiva  A/((1/alfa) -green tensor)
    mismo sistema de coordenadas que el paper 370 (z hacia abajo)
    """
    
    omega = omegac*c
    # k0 = omegac #=omega/c
    green_tensor_f = green_tensor_function(omegac,hbmu,z_plane,epsi1,epsi2,epsi_h,hbargama_sigma,px,py,pz)
    den = omega0 - omega - 1j*hbargama_alfa/2 - A*green_tensor_f
    
    return A/den

#%%
