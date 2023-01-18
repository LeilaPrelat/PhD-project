#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila

polarizability alpha_E
"""
from scipy import integrate
import numpy as np
import sys
import os 

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_graphene =  path_basic.replace('/' + 'plane','')
#print('Importar modulos necesarios para este codigo')

try:
    sys.path.insert(1, path_graphene + '/' + 'sigma_graphene')
    from graphene_sigma import sigma_DL
except ModuleNotFoundError:
    print('graphene_sigma.py no se encuentra en ' + path_graphene + '/' + 'sigma_graphene')

try:
    sys.path.insert(1, path_basic)
    from constants import constantes
except ModuleNotFoundError:
    print('constants_plane.py no se encuentra en ' + path_basic)

pi,hb,c,alfac,mu1,mu2 = constantes()
aux = c*hb

#%%

def alpha_function(omegac,hbmu,z_plane,epsi1,epsi2,hbargama):
    """
    Parameters
    ----------
    omegac : omega/c = k0 en 1/micrometros    
    hbmu: potencial quimico del grafeno en eV
    z_plane: radio del cilindro en micrometros
    epsi1 : epsilon del medio de arriba del plano
    epsi2 : epsilon del medio de abajo del plano
    hbargamma : collision frequency in eV del sigma (grafeno)
    Returns
    -------
    polarizabilty Eq 110
    mismo sistema de coordenadas que el paper 370 (z hacia abajo)
    """
    
    E = omegac*aux
    # k0 = omegac #=omega/c
    # x_y = ky/k0
    n1 = epsi1*mu1
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1
    k1_3 = k1**3
    b = epsi2*mu2/n1
    z_plane_barra = z_plane*k1

    sigmatot = sigma_DL(E,hbmu,hbargama)
    cond3 = sigmatot*alfac*4*pi  #sigma devuelve la conductividad teniendo que multiplicar por alfac*c ---> no hay que dividir por c

    term2 = lambda u: (b-1+u**2)**(1/2)
    expA = lambda u: np.exp(1j*(u-np.sqrt(term2(u)))*z_plane_barra) 
    expB = lambda u: np.exp(-(u-np.sqrt(term2(u)))*z_plane_barra)

    num1 = lambda u: 2*mu2*u**2/np.sqrt(1-u**2)
    den1 = lambda u: mu2*u + mu1*term2(u) - cond3*mu2*mu1/cte1 
    f4 = lambda u: num1(u)/den1(u) 

    F4A_re = lambda u: np.real(f4(u)*expA(u))
    F4B_re = lambda u: np.real(f4(u)*expB(u))

    F4A_im = lambda u: np.imag(f4(u)*expA(u))
    F4B_im = lambda u: np.imag(f4(u)*expB(u))

    int4A_re,err = integrate.quad(F4A_re, 0, 1)
    int4B_re,err = integrate.quad(F4B_re, 0, np.inf)
    int4A_im,err = integrate.quad(F4A_im, 0, 1)
    int4B_im,err = integrate.quad(F4B_im, 0, np.inf)
    
    cte = 3*k1_3/2

    final = int4A_re + 1j*int4A_im + int4B_re + 1j*int4B_im
    
    return final*cte

#%%
