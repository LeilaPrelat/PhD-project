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
    from hBn_PP import hBn_lambda_p,hBn_Rp,epsilon_x
except ModuleNotFoundError:
    print('graphene_sigma.py no se encuentra en ' + path_constants)

try:
    sys.path.insert(1, path_constants)
    from dipole_moment_dif_sign import dipole_moment_ana_resonance_v1, dipole_moment_num_resonance,dipole_moment_pole_aprox_resonance_v1, dipole_moment_pole_aprox_resonance_v1_for_decay_rate,dipole_moment_anav1_for_decay_rate_resonance
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
def potential_inf_dipoles_num(omegac,epsi_silica,d_nano_film,d_thickness_disk_nano,D_disk_nano,int_v,zp,a,b,n):     
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
   
   
    Rp = 1
    epsi_x = epsilon_x(E)
    epsi_HBN_par = epsi_x
    d_micro = d_nano_film*1e-3
    alfa_p = epsi_silica(E)*2/(omegac*d_micro*(epsi_HBN_par-1))
    kp = alfa_p*omegac
   
   
    px,py,pz  = dipole_moment_num_resonance(omegac,epsi_silica,d_nano_film,d_thickness_disk_nano,D_disk_nano,int_v,b,zp)  
#    list_dipoles = np.linspace(-Nmax,Nmax,2*Nmax + 1)
#            
    kx = omegac*int_v + 2*np.pi*n/a     
    den = np.sqrt(kp**2 - kx**2)
   # return np.imag(final_2*cte*kp*np.cos(theta))
   
   
    phi_n = np.exp(-2*kp*zp)*kp*(px*kx/den + py + 1j*pz*kp/den )/(2*np.pi*a)

    return phi_n













# normalizado con el paper 149 
def potential_inf_dipoles_pole_aprox(omegac,epsi_silica,d_nano_film,d_thickness_disk_nano,D_disk_nano,int_v,zp,a,b,n):     
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
   
    Rp = 1
    epsi_x = epsilon_x(E)
    epsi_HBN_par = epsi_x
    d_micro = d_nano_film*1e-3
    alfa_p = epsi_silica(E)*2/(omegac*d_micro*(epsi_HBN_par-1))
    kp = alfa_p*omegac
   
   
    px,py,pz  = dipole_moment_pole_aprox_resonance_v1(omegac,epsi_silica,d_nano_film,d_thickness_disk_nano,D_disk_nano,int_v,b,zp)  
#    list_dipoles = np.linspace(-Nmax,Nmax,2*Nmax + 1)
#            
    kx = omegac*int_v + 2*np.pi*n/a     
    den = np.sqrt(kp**2 - kx**2)
   # return np.imag(final_2*cte*kp*np.cos(theta))
   
   
    phi_n = np.exp(-2*kp*zp)*kp*(px*kx/den + py + 1j*pz*kp/den )/(2*np.pi*a)

    return phi_n













# normalizado con el paper 149 
def potential_inf_dipoles_ana(omegac,epsi_silica,d_nano_film,d_thickness_disk_nano,D_disk_nano,int_v,zp,a,b,n):     
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
    print(E, epsi_silica(E))
    
    Rp = 1
    epsi_x = epsilon_x(E)
    epsi_HBN_par = epsi_x
    d_micro = d_nano_film*1e-3
    alfa_p = epsi_silica(E)*2/(omegac*d_micro*(epsi_HBN_par-1))
    kp = alfa_p*omegac
   
   
    px,py,pz  = dipole_moment_anav1_for_decay_rate_resonance(omegac,epsi_silica,d_nano_film,d_thickness_disk_nano,D_disk_nano,int_v,b,zp)  
#    list_dipoles = np.linspace(-Nmax,Nmax,2*Nmax + 1)
#            
    kx = omegac*int_v + 2*np.pi*n/a     
    den = np.sqrt(kp**2 - kx**2)
   # return np.imag(final_2*cte*kp*np.cos(theta))
   
   
    phi_n = -np.exp(-2*kp*zp)*kp*(px*kx/den + py + 1j*pz*kp/den )/(2*np.pi*a)

    return phi_n

