#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila
"""
import numpy as np
import sys
import os 
from scipy import special

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('potential_field/potential_and_electric_field/potential_and_electric_field_simple_version/potential_1_dipolo_analytical_vs_num','')
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


def compare_fresnel_vs_pole_approx(omegac,epsi1,epsi2,hbmu,hbgama,u): # u es alfa parallel y rp es fresnel coefficient
    """
    Parameters
    ----------
    omegac : omega/c = k0 en 1/micrometros    
    epsi1 : epsilon del medio de arriba del plano
    epsi2 : epsilon del medio de abajo del plano
    hbmu : chemical potential in eV  
    hbgama : collision frequency in eV
    int_v : fraccion de la vel de la c que vale la velocidad uniforme del electron en x, si v = c/int
    zp : coordenada zp del plano, zp > 0  (dipolo en 0 y eje de z apunta hacia abajo) 
    x : punto x donde miramos el campo en micrones    
    Returns
    -------
    alfa effectivo en QE approx
    """

    E = omegac*aux 
#    k0 = omegac
#    n1 = epsi1*mu1
#    cte1 = np.sqrt(n1)
#    k1 = omegac*cte1
#    k1_2 = k1**2

    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)

    rp_num = epsi2*1j*u - epsi1*1j*u - cond*u**2
    rp_den = epsi2*1j*u + epsi1*1j*u - cond*u**2
    rp_fresnel = rp_num/rp_den  
    
    
    alfa_p = 1j*(epsi1 + epsi2)/cond
 
    Rp = 2*epsi1/(epsi1 + epsi2)   
    
    rp_pp  =  Rp*alfa_p/(u - alfa_p)

    final = 1 - rp_pp/rp_fresnel

    return np.abs(final)

#%%