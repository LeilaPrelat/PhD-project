#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on March 19 20:56:58 2020

@author: leila

EELS
el dipolo esta en x,y,z = 0,0,0

"""
import numpy as np
import os 
import sys
from scipy import integrate

sys.setrecursionlimit(10000)

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


# try:
#     sys.path.insert(1, path_basic)
#     from polarizability import alpha_function
# except ModuleNotFoundError:
#     print('polarizability.py no se encuentra en ' + path_basic)

# try:
#     sys.path.insert(1, path_basic)
#     from green_tensor_ICFO import green_tensor
# except ModuleNotFoundError:
#     print('green_tensor_ICFO.py no se encuentra en ' + path_basic)
    
try:
    sys.path.insert(1, path_basic)
    from constants import constantes
except ModuleNotFoundError:
    print('constants.py no se encuentra en ' + path_basic)

pi,hb,c,alfac,mu1,mu2 = constantes()
aux = c*hb

#%%
    
#print('Definir la funcion Hz para el medio 1 y el 2')

def fieldEELS(omegac,hbmu,y,z_electron,z_plane,v,epsi1,epsi2,hbargama,px,py,pz):
    """
    Parameters
    ----------
    omegac : omega/c = k0 en 1/micrometros    
    hbmu: potencial quimico del grafeno en eV
    y: coordenada y 
    z_electron : posicion del electron en micrometros (<0)
    z_plane: posicion del plano en micrometros (>0)
    v : velocidad del electron (eje x, paralelo al plano)
    epsi1 : epsilon del medio de arriba del plano
    epsi2 : epsilon del medio de abajo del plano
    hbargamma : collision frequency in eV del sigma (grafeno)
    px : coordenada x del dipolo 
    py : coordenada y del dipolo 
    pz : coordenada z del dipolo 
    Returns
        Ex, Ey, Ez del campo inducido por un dipolo
        el dipolo esta en x,y,z = 0,0,0
    -------
    """
    E = omegac*aux
    # k0 = omegac #=omega/c
    omega = omegac*c
    # x_y = ky/k0
    n1 = epsi1*mu1
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1
    k1_3 = k1**3
    b = epsi2*mu2/n1
    v_monio_2 = (omega/(v*k1))**2


    y_barra = y*k1
    z_plane_barra = z_plane*k1
    z_electron_barra = z_electron*k1
    
#    sigmatot = sigma_DL(E,hbmu,hbargama)
#    cond3 = sigmatot*alfac*4*pi  #sigma devuelve la conductividad teniendo que multiplicar por alfac*c ---> no hay que dividir por c
    cond3 = 0
    
    exp_fresnel = lambda u: np.exp(1j*2*u*z_plane_barra) 
    term2 = lambda u: (b-1+u**2)**(1/2)
    num_rs = lambda u: u - term2(u) - cond3/cte1
    den_rs = lambda u: u + term2(u) + cond3/cte1      
    fresnel_rs = lambda u: num_rs(u)*exp_fresnel(u)/den_rs(u) 

    num_rp = lambda u: epsi2*u - epsi1*term2(u) + cond3*term2(u)*u*cte1
    den_rp = lambda u: epsi2*u + epsi1*term2(u) + cond3*term2(u)*u*cte1
    fresnel_rp = lambda u: num_rp(u)*exp_fresnel(u)/den_rp(u) 

    alfa_y = lambda u: np.sqrt(1 - v_monio_2 - u**2)
    alfa_x = (omega/v)/k1
    # alfa_parallel = lambda u: np.sqrt(alfa_y(u)**2 + alfa_x**2)
    alfa_parallel_2 = lambda u: alfa_y(u)**2 + alfa_x**2
    
    term1 = lambda u: (px*alfa_y(u) - py*alfa_x)*(alfa_y(u))/alfa_parallel_2(u)
    term1_final = lambda u: term1(u)*(1+fresnel_rs(u))

    term2 = lambda u: (px*alfa_x + py*alfa_y(u))*(1-fresnel_rp(u))
    term2_final = lambda u: term2(u)*(u*alfa_x)**2/alfa_parallel_2(u)

    term3 = lambda u: (1+fresnel_rp(u))*pz*u*alfa_x
    
    terminos = lambda u: term1_final(u) + term2_final(u) + term3(u)
    
    exp_1A = lambda u:  np.exp(1j*u*z_electron_barra) 
    exp_1B = lambda u:  np.exp(-u*z_electron_barra) 
    exp_2 = lambda u:  np.exp(1j*alfa_y(u)*y_barra) 


    funcion_finalA = lambda u: np.real(terminos(u)*exp_1A*(u)*exp_2(u)/alfa_y(u))
    funcion_finalB = lambda u: np.real(terminos(u)*exp_1B*(u)*exp_2(u)/alfa_y(u))
    
    cte_aux = np.sqrt(1-v**2)
    
    int1A_re,err = integrate.quad(funcion_finalA, 0, cte_aux)
    int1B_re,err = integrate.quad(funcion_finalB, 0, np.inf)
    
    return (int1A_re + int1B_re)*k1_3
    
#%%        
