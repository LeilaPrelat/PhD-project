#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila

green tensor 
"""
from scipy import integrate
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
    from graphene_sigma import sigma_DL
except ModuleNotFoundError:
    print('graphene_sigma.py no se encuentra en ' + path_basic + '/' + 'sigma_graphene')

try:
    sys.path.insert(1, path_basic)
    from constants import constantes
except ModuleNotFoundError:
    print('constants.py no se encuentra en ' + path_basic)

pi,hb,c,alfac,mu1,mu2 = constantes()
aux = c*hb

#%%

def green_tensor_function(omegac,hbmu,z_plane,epsi1,epsi2,epsi_h,hbargama_sigma,px,py,pz):
    """
    Parameters
    ----------
    omegac : omega/c = k0 en 1/micrometros    
    hbmu: potencial quimico del grafeno en eV
    z_plane: radio del cilindro en micrometros
    epsi1 : epsilon del medio de arriba del plano
    epsi2 : epsilon del medio de abajo del plano
    epsi_h = epsilon del host medium (material del plano)
    hbargama_sigma : collision frequency in eV del sigma (grafeno)
    px : coordenada x del dipolo 
    py : coordenada y del dipolo 
    pz : coordenada z del dipolo 
    Returns
    -------
    Gxx = Gyy = G_parallel, Gzz = G_perp (self interaction of the dipole)
    con z hacia abajo (convencion del paper)
    """
    E = omegac*aux
    k0 = omegac #=omega/c
    # x_y = ky/k0
    n1 = epsi1*mu1
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1
    k1_3 = (k0*cte1)**3
    b = epsi2*mu2/n1
    z_plane_barra = z_plane*k1

    sigmatot = sigma_DL(E,hbmu,hbargama_sigma)
    cond3 = sigmatot*alfac*4*pi  #sigma devuelve la conductividad teniendo que multiplicar por alfac*c ---> no hay que dividir por c

    expA = lambda u: np.exp(1j*2*u*z_plane_barra) 
    expB = lambda u: np.exp(-2*u*z_plane_barra)
    term2 = lambda u: (b-1+u**2)**(1/2)
    num1 = lambda u: u - term2(u) - cond3/cte1
    den1 = lambda u: u + term2(u) + cond3/cte1      
    f1 = lambda u: num1(u)/den1(u) 

    num2 = lambda u: u*epsi2 - term2(u)*epsi1 + cond3*term2(u)*u*cte1
    den2 = lambda u: u*epsi2 + term2(u)*epsi1 + cond3*term2(u)*u*cte1  
    
    f2 = lambda u: (u**2)*(num2(u)/den2(u))
    f3 = lambda u: (1-u**2)*(num2(u)/den2(u))

    F1A_re = lambda u: np.real(f1(u)*expA(u))
    F1B_re = lambda u: np.real(f1(u)*expB(u))
    
    F2A_re = lambda u: np.real(f2(u)*expA(u))
    F2B_re = lambda u: np.real(f2(u)*expB(u))
    
    F3A_re = lambda u: np.real(f3(u)*expA(u))
    F3B_re = lambda u: np.real(f3(u)*expB(u))
    
    F1A_im = lambda u: np.imag(f1(u)*expA(u))
    F1B_im = lambda u: np.imag(f1(u)*expB(u))
    
    F2A_im = lambda u: np.imag(f2(u)*expA(u))
    F2B_im = lambda u: np.imag(f2(u)*expB(u))
    
    F3A_im = lambda u: np.imag(f3(u)*expA(u))
    F3B_im = lambda u: np.imag(f3(u)*expB(u))

    int1A_re,err = integrate.quad(F1A_re, 0, 1)
    int1B_re,err = integrate.quad(F1B_re, 0, np.inf)
    int1A_im,err = integrate.quad(F1A_im, 0, 1)
    int1B_im,err = integrate.quad(F1B_im, 0, np.inf)

    int2A_re,err = integrate.quad(F2A_re, 0, 1)
    int2B_re,err = integrate.quad(F2B_re, 0, np.inf)
    int2A_im,err = integrate.quad(F2A_im, 0, 1)
    int2B_im,err = integrate.quad(F2B_im, 0, np.inf)
    
    int3A_re,err = integrate.quad(F3A_re, 0, 1)
    int3B_re,err = integrate.quad(F3B_re, 0, np.inf)
    int3A_im,err = integrate.quad(F3A_im, 0, 1)
    int3B_im,err = integrate.quad(F3B_im, 0, np.inf)
    
    int1A = int1A_re + 1j*int1A_im
    int1B = int1B_re + 1j*int1B_im
    
    int2A = int2A_re + 1j*int2A_im
    int2B = int2B_re + 1j*int2B_im
    
    int3A = int3A_re + 1j*int3A_im
    int3B = int3B_re + 1j*int3B_im
    
    cte = 0.5*1j*k1_3/epsi_h
    Gxx =  (int1A + int1B + int2A + int2B)*cte
    Gzz = (int3A + int3B)*cte
    
    mod_p = np.sqrt(px**2 + py**2 + pz**2)
    nx,ny,nz = px/mod_p, py/mod_p , pz/mod_p
    
    rtaf = Gxx*px*nx + Gxx*py*ny + Gzz*pz*nz 
    
    return rtaf

#%%
