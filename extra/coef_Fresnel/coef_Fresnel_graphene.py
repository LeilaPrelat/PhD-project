#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila

coeficientes de Fresnel para el grafeno

"""

import numpy as np
import sys
import os 

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_graphene =  path_basic.replace('/' + 'coef_Fresnel','')
#print('Importar modulos necesarios para este codigo')

try:
    sys.path.insert(1, path_graphene)
    from graphene_sigma import sigma
except ModuleNotFoundError:
    print('graphene_sigma.py no se encuentra en ' + path_graphene + '/' + 'sigma_graphene')

try:
    sys.path.insert(1, path_graphene)
    from constants import constantes
except ModuleNotFoundError:
    print('constants_plane.py no se encuentra en ' + path_graphene)

pi,hb,c,alfac,mu1,mu2 = constantes()
aux = c*hb

#%%

def coefTM_TE(omegac,epsi1,epsi2,hbmu,hbargama,alpha_parallel):
    """
    Parameters
    ----------
    omegac : omega/c = k0 en 1/micrometros    
    epsi1 : epsilon del medio de arriba del plano
    epsi2 : epsilon del medio de abajo del plano
    hbmu: potencial quimico del grafeno en eV
    hbgama : collision frequency in eV
    alpha_parallel:  alpha_parallel
    Returns
    -------
    2 coef: refleado TM y refleado TE.
    La posicion del plano esta en 0 (no depende de z_plane)
    """
    E = omegac*aux
    # k0 = omegac #=omega/c

    n1 = (epsi1*mu1)**(1/2)
    # k1 = k0*n1
    
    if alpha_parallel<=1:
        alphaz1 = np.sqrt(1-alpha_parallel**2)
    else:
        alphaz1 = 1j*np.sqrt(alpha_parallel**2 - 1)
            
    b = epsi2*mu2/(epsi1*mu1)
    if alpha_parallel<=b:
        alphaz2 = np.sqrt(b-alpha_parallel**2)
    else:
        alphaz2 = 1j*np.sqrt(alpha_parallel**2 - b)
        
    sigmatot = sigma(E,hbmu,hbargama)
    cond3 = sigmatot*alfac*4*pi  #sigma devuelve la conductividad teniendo que multiplicar por alfac*c ---> no hay que dividir por c
     
    num_rs = mu2*alphaz1 - mu1*alphaz2 - cond3/n1
    den_rs = mu2*alphaz1 + mu1*alphaz2 + cond3/n1
    coef_rs = num_rs/den_rs
    
    num_rp = epsi2*alphaz1 - epsi1*alphaz2 - alphaz2*alphaz1*cond3*n1
    den_rp = epsi2*alphaz1 + epsi1*alphaz2 + alphaz2*alphaz1*cond3*n1
    coef_rp = num_rp/den_rp

# FORMULAS DEL OVERLEAF : 70a (r_s) y 70b (r_p)
    
    return coef_rs, coef_rp

#%%
