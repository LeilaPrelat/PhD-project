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
path_constants =  path_basic.replace('/potential_field/infinite_dipoles','')
#print('Importar modulos necesarios para este codigo')

try:
    sys.path.insert(1, path_constants)
    from hBn_PP import hBn_lambda_p,hBn_Rp
except ModuleNotFoundError:
    print('graphene_sigma.py no se encuentra en ' + path_constants)

try:
    sys.path.insert(1, path_constants)
    from dipole_moment import dipole_moment_anav2_resonance, dipole_moment_num_resonance,dipole_moment_pole_aprox_resonance
except ModuleNotFoundError:
    print('dipole_moment.py no se encuentra en ' + path_constants)
    
try:
    sys.path.insert(1, path_constants)
    from constants import constantes
except ModuleNotFoundError:
    print('constants.py no se encuentra en ' + path_constants)

pi,hb,c,alfac,mu1,mu2 = constantes()
aux = c*hb

#%%


# normalizado con el paper 149 
def potential_inf_dipoles_num(omegac,epsi1,epsi3,d_nano,int_v,zp,a,b,n):     
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

#    x, y, z = 0,0,0
    E = omegac*aux
   
    d_micros = d_nano*1e-3
    Rp = hBn_Rp(E,epsi1,epsi3)
    lambda_p_v = hBn_lambda_p(E,epsi1,epsi3)*d_micros
    kp = 2*np.pi/lambda_p_v   
    
   
    px,py,pz  = dipole_moment_num_resonance(omegac,epsi1,epsi3,d_nano,int_v,b,zp)  
#    list_dipoles = np.linspace(-Nmax,Nmax,2*Nmax + 1)
#            
    kx = omegac*int_v + 2*np.pi*n/a     
    den = np.sqrt(kp**2 - kx**2)
   # return np.imag(final_2*cte*kp*np.cos(theta))
    phi_n = np.exp(-2*kp*zp)*Rp*kp*(px*kx/den + py + 1j*pz*kp/den )/(2*np.pi*a)

    return phi_n












# normalizado con el paper 149 
def potential_inf_dipoles_pole_aprox(omegac,epsi1,epsi3,d_nano,int_v,zp,a,b,n):     
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

#    x, y, z = 0,0,0
    E = omegac*aux
   
    d_micros = d_nano*1e-3
    Rp = hBn_Rp(E,epsi1,epsi3)
    lambda_p_v = hBn_lambda_p(E,epsi1,epsi3)*d_micros
    kp = 2*np.pi/lambda_p_v   
    
   
    px,py,pz  = dipole_moment_pole_aprox_resonance(omegac,epsi1,epsi3,d_nano,int_v,b,zp)  
#    list_dipoles = np.linspace(-Nmax,Nmax,2*Nmax + 1)
#            
    kx = omegac*int_v + 2*np.pi*n/a     
    den = np.sqrt(kp**2 - kx**2)
   # return np.imag(final_2*cte*kp*np.cos(theta))
    phi_n = np.exp(-2*kp*zp)*Rp*kp*(px*kx/den + py + 1j*pz*kp/den )/(2*np.pi*a)

    return phi_n












# normalizado con el paper 149 
def potential_inf_dipoles_ana(omegac,epsi1,epsi3,d_nano,int_v,zp,a,b,n):     
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

#    x, y, z = 0,0,0
    E = omegac*aux
   
    d_micros = d_nano*1e-3
    Rp = hBn_Rp(E,epsi1,epsi3)
    lambda_p_v = hBn_lambda_p(E,epsi1,epsi3)*d_micros
    kp = 2*np.pi/lambda_p_v   
    
   
    px,py,pz  = dipole_moment_anav2_resonance(omegac,epsi1,epsi3,d_nano,int_v,b,zp)  
#    list_dipoles = np.linspace(-Nmax,Nmax,2*Nmax + 1)
#            
    kx = omegac*int_v + 2*np.pi*n/a     
    den = np.sqrt(kp**2 - kx**2)
   # return np.imag(final_2*cte*kp*np.cos(theta))
    phi_n = np.exp(-2*kp*zp)*Rp*kp*(px*kx/den + py + 1j*pz*kp/den )/(2*np.pi*a)

    return phi_n

