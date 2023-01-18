
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila

campo externo directo + reflejado numerico con la convencion de z hacia abajo
en z = 0
No hay solucion analitica porque es para cualquier x,y (solo hay sol analitica en x=y=0)
Dipolo en el 0,0
"""
import numpy as np
import sys
import os 

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/polarizability/wrong','')
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

def alpha_function(omegac,epsi1,epsi2,hbargama,hbmu,z_plane,u):
    """
    Parameters
    ----------
    omegac : omega/c = k0 en 1/micrometros    
    epsi1 : epsilon del medio de arriba del plano
    epsi2 : epsilon del medio de abajo del plano
    hbargamma : collision frequency in eV del sigma (grafeno)
    hbmu: potencial quimico del grafeno en eV
    z_plane: radio del cilindro en micrometros
    u : alpha parallel
    Returns
    -------
    funcion dentro de polarizabilty Eq 110 (antes de integrar en alpha_parallel)
    mismo sistema de coordenadas que el paper 370 (z hacia abajo)
    """
    
    E = omegac*aux
    n1 = epsi1*mu1
    n2 = epsi2*mu2
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1
    k1_3 = k1**3
    
    aux2 = n2/n1
    
    def alpha_z1(u): 
        return np.sqrt(1-u**2) if u<1 else 1j*np.sqrt(u**2-1)
    
    def alpha_z2(u):
        return np.sqrt(aux2-u**2) if u<aux2 else 1j*np.sqrt(u**2-aux2)
    
    sigmatot = sigma_DL(E,hbmu,hbargama)
    cond3 = sigmatot*alfac*4*pi  #sigma devuelve la conductividad teniendo que multiplicar por alfac*c ---> no hay que dividir por c

    def function(u):
        exp = np.exp(1j*k1*(alpha_z1(u)-alpha_z2(u))*z_plane) 

        num = 2*mu2*alpha_z1(u)*2
        den = mu2*alpha_z1(u) + mu1*alpha_z2(u) - cond3*mu2*mu1/cte1 
        func = num/den 

        F4A_re = np.real(func*exp)
        F4A_im = np.imag(func*exp)

#    cota_sup1A = 1000
#    cota_sup2A = 80*k1

#    int4_re,err = integrate.quad(F4A_re, 0, cota_sup1A)
#    int4_im,err = integrate.quad(F4A_im, 0, cota_sup1A)

        cte = 3*k1_3/2

        final = F4A_re + 1j*F4A_im
    
        return final*cte

    rta = function(u)
    return rta

#%%


def alpha_functionQE(omegac,epsi1,epsi2,hbargama,hbmu,z_plane,u):
    """
    Parameters
    ----------
    omegac : omega/c = k0 en 1/micrometros    
    epsi1 : epsilon del medio de arriba del plano
    epsi2 : epsilon del medio de abajo del plano
    hbargamma : collision frequency in eV del sigma (grafeno)
    hbmu: potencial quimico del grafeno en eV
    z_plane: radio del cilindro en micrometros
    u : alpha parallel
    Returns
    -------
    funcion dentro de polarizabilty Eq 110 (antes de integrar en alpha_parallel)
    mismo sistema de coordenadas que el paper 370 (z hacia abajo)
    bajo aprox QE
    """
    
    E = omegac*aux
    n1 = epsi1*mu1
 #   n2 = epsi2*mu2
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1
    k1_3 = k1**3
    
#    aux2 = n2/n1
    
    def alpha_z1(u): 
        return 1j*u
    
    def alpha_z2(u):
        return 1j*u
    
    sigmatot = sigma_DL(E,hbmu,hbargama)
    cond3 = sigmatot*alfac*4*pi  #sigma devuelve la conductividad teniendo que multiplicar por alfac*c ---> no hay que dividir por c

    def function(u):
        exp = np.exp(1j*k1*(alpha_z1(u)-alpha_z2(u))*z_plane) 

        num = 2*mu2*alpha_z1(u)*2
        den = mu2*alpha_z1(u) + mu1*alpha_z2(u) - cond3*mu2*mu1/cte1 
        func = num/den 

        F4A_re = np.real(func*exp)
        F4A_im = np.imag(func*exp)

#    cota_sup1A = 1000
#    cota_sup2A = 80*k1

#    int4_re,err = integrate.quad(F4A_re, 0, cota_sup1A)
#    int4_im,err = integrate.quad(F4A_im, 0, cota_sup1A)

        cte = 3*k1_3/2

        final = F4A_re + 1j*F4A_im
    
        return final*cte

    rta = function(u)
    return rta

#%%
