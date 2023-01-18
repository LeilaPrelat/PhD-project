#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila
"""
import numpy as np
import sys
import os 

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/potential_field/infinite_dipoles','')
#print('Importar modulos necesarios para este codigo')

try:
    sys.path.insert(1, path_constants)
    from dipole_moment import dipole_moment_anav2
except ModuleNotFoundError:
    print('dipole_moment.py no se encuentra en ' + path_basic)

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

def Gamma_film(omegac,epsi1,epsi2,hbmu,hbgama,int_v,zp,a,b,x,y,z,omega0,kappa_factor_omega0,kappa_r_factor,n):     
    """    
    Parameters
    ----------
    omegac : omega/c = k0 en 1/micrometros    
    epsi1 : epsilon del medio de arriba del plano
    epsi2 : epsilon del medio de abajo del plano
    hbmu : chemical potential in eV  
    hbgama : collision frequency in eV
    z : coordenada z
    xD : coordenada x del dipolo 
    yD : coordenada y del dipolo
    zD : coordenada z del dipolo 
    zp : posicion del plano (>0)
    px : coordenada x del dipolo 
    py : coordenada y del dipolo
    pz : coordenada z del dipolo
    Returns
    -------
    formula del potencial electric con QE approximation, rp con 
    aproximacion del polo y con aprox de principal value para las integrales
    con rp
    """

    E = omegac*aux
#    k0 = omegac #=omega/c
    # x_y = ky/k0
    n1 = epsi1*mu1
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1
 #   n_v1 = int_v/cte1
#    k1_2 = k1**2
   
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    Rp = 2*epsi1/(epsi1 + epsi2)
    alfa_p = 1j*(epsi1 + epsi2)/(cond*cte1)
    kp = alfa_p*k1
    

    px, py, pz = dipole_moment_anav2(omegac,epsi1,epsi2,hbmu,hbgama,int_v,b,zp,omega0,kappa_factor_omega0,kappa_r_factor)

    v = c/int_v
    cte_p = alfac*c/((2*pi*v)**2)
    px_2, py_2, pz_2 = px*cte_p, py*cte_p, pz*cte_p  # cte que falta + el hbar --> aparece alfac*c
        
    kx = omegac*int_v - 2*n*pi/a
        
    if np.real(kp) > kx: 
        termy = np.sqrt(kp**2 - kx**2)    
    else:
        termy = 1j*np.sqrt(kx**2 - kp**2)   
        
    exp_y = np.exp(1j*termy*y)
    exp_z = np.exp(-kp*(2*zp - z))
    exp_x = np.exp(1j*kx*x)    
    
    
    cte = Rp*kp/(2*pi*a)
    term = exp_x*( px_2*kx/termy + py_2 + pz_2*kp/termy )*exp_y*exp_z

    Ex = cte*1j*kx*term
    Ey = cte*1j*termy*term
    Ez = cte*kp*term
    
    rta = Ex*np.conjugate(px) + Ey*np.conjugate(py) +  Ez*np.conjugate(pz)

    return rta*2 

#%%
