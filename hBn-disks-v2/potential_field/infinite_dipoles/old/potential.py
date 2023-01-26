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
from scipy import integrate

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/potential_field/infinite_dipoles','')
#print('Importar modulos necesarios para este codigo')

try:
    sys.path.insert(1, path_constants)
    from hBn_PP import epsilon_x
except ModuleNotFoundError:
    print('hBn_PP.py no se encuentra en ' + path_constants)


try:
    sys.path.insert(1, path_constants)
    from dipole_moment import dipole_moment_ana_resonance_v1,dipole_moment_ana_resonance_v2,dipole_moment_pole_aprox_resonance_v1, dipole_moment_num_resonance,dipole_moment_pole_aprox_resonance_v2
except ModuleNotFoundError:
    print('green_self_image.py no se encuentra en ' + path_basic)


try:
    sys.path.insert(1, path_constants)
    from constants import constantes
except ModuleNotFoundError:
    print('constants.py no se encuentra en ' + path_constants)

pi,hb,c,alfac,mu1,mu2 = constantes()
aux = c*hb


#%%
    
#     rp = lambda u: alpha_parallel(u)/(alpha_parallel(u) - alfa_p)

def potential_ana_resonance_v1(omegac,epsi_silica,d_nano,int_v,b,zp,x,y,z,a,N):     ### usando 
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
    px,py,pz en unidades de k*alfa_eff
    """

    E = omegac*aux
#    k0 = omegac #=omega/c
    # x_y = ky/k0

#    k1_2 = (k0*cte1)**2
 #   n_v1 = int_v/cte1

    epsi_x = epsilon_x(E)
    epsi_HBN_par = epsi_x
    
    d_micro = d_nano*1e-3
    alfa_p = epsi_silica(E)*2/(omegac*d_micro*(epsi_HBN_par-1))
    kp = alfa_p*omegac


    px, py, pz  =  dipole_moment_pole_aprox_resonance_v1(omegac,epsi_silica,d_nano,int_v,b,zp)

    kx = omegac*int_v + 2*np.pi*N/a
    arg_z_y = np.sqrt(np.abs(z)**2 + np.abs(y)**2)
    K0 =  special.kv(0,kx*arg_z_y)
    K1 =  special.kv(1,kx*arg_z_y)
    
    term_kp = np.sqrt(kp**2) + kp

    term1 = -2*1j*px*kx*K0 + (py*np.abs(y) + pz*np.sign(z)*np.abs(z))*2*kx*K1/arg_z_y 
    
    expo_z = np.exp(-np.sqrt(kp**2)*(2*zp-z) )
    expo_y = np.exp(-np.sqrt(kp**2 - kx**2)*np.abs(y) )
    
    term2 = -px*kx*np.pi*term_kp*expo_z*expo_y/(np.sqrt(kp**2 - kx**2))
    
    term3 = -py*np.pi*term_kp*expo_z*expo_y
    
    term4 = -pz*1j*np.pi*term_kp*np.sqrt(kp**2)*expo_z*expo_y/(np.sqrt(kp**2 - kx**2))
    
    cte_final = np.exp(1j*kx*x)/((2*np.pi)**2*a)
    
    return (term1 + term2 + term3 + term4)*cte_final


def potential_ana_resonance_v2(omegac,epsi_silica,d_nano,int_v,b,zp,x,y,z,a,N):     ### usando 
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
    px,py,pz en unidades de k*alfa_eff
    """
    E = omegac*aux
    epsi_x = epsilon_x(E)
    epsi_HBN_par = epsi_x
    
    d_micro = d_nano*1e-3
    alfa_p = epsi_silica(E)*2/(omegac*d_micro*(epsi_HBN_par-1))
    kp = alfa_p*omegac


    px, py, pz  =  dipole_moment_pole_aprox_resonance_v1(omegac,epsi_silica,d_nano,int_v,b,zp)

    kx = omegac*int_v + 2*np.pi*N/a
    arg_z_y = np.sqrt(np.abs(z)**2 + np.abs(y)**2)
    K0 =  special.kv(0,kx*arg_z_y)
    K1 =  special.kv(1,kx*arg_z_y)
    
    term_kp = np.sqrt(kp**2) + kp

    term1 = -2*1j*px*kx*K0 + (py*np.abs(y) + pz*np.sign(z)*np.abs(z))*2*kx*K1/arg_z_y 
    
    expo_z = np.exp(-np.sqrt(kp**2)*(2*zp-z) )
    expo_y = np.exp(-np.sqrt(kp**2 - kx**2)*np.abs(y) )
    
    term2 = -px*kx*np.pi*term_kp*expo_z*expo_y/(np.sqrt(kp**2 - kx**2))
    
    term3 = -py*np.pi*term_kp*expo_z*expo_y
    
    term4 = -pz*1j*np.pi*term_kp*np.sqrt(kp**2)*expo_z*expo_y/(np.sqrt(kp**2 - kx**2))
    
    cte_final = np.exp(1j*kx*x)/((2*np.pi)**2*a)
    
    return (term1 + term2 + term3 + term4)*cte_final


#%%


def potential_num_resonance(omegac,epsi_silica,d_nano,int_v,b,zp,x,y,z,a,N):     
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
    epsi_x = epsilon_x(E)
    epsi_HBN_par = epsi_x
    
    d_micro = d_nano*1e-3
    alfa_p = epsi_silica(E)*2/(omegac*d_micro*(epsi_HBN_par-1))
    kp = alfa_p*omegac


    px, py, pz  =  dipole_moment_num_resonance(omegac,epsi_silica,d_nano,int_v,b,zp)

    kx = omegac*int_v + 2*np.pi*N/a
    arg_z_y = np.sqrt(np.abs(z)**2 + np.abs(y)**2)
    K0 =  special.kv(0,kx*arg_z_y)
    K1 =  special.kv(1,kx*arg_z_y)
    
#    term_kp = np.sqrt(kp**2) + kp

    term1 = -2*1j*px*kx*K0 + (py*np.abs(y) + pz*np.sign(z)*np.abs(z))*2*kx*K1/arg_z_y 
    
    alfa_x = int_v +  2*np.pi*N/(a*omegac)
    
    expo_z = lambda u : np.exp(-np.sqrt(alfa_x**2 + u**2)*omegac*(2*zp-z) )
    expo_y = lambda u : np.exp(1j*u*omegac*np.abs(y))
    
    term2_re = lambda u: expo_z(u)*expo_y(u)*np.real(1/(np.sqrt(alfa_x**2 + u**2) - alfa_p ))
    term2_im = lambda u: expo_z(u)*expo_y(u)*np.imag(1/(np.sqrt(alfa_x**2 + u**2) - alfa_p ))
    
    term3_re = lambda u: u*expo_z(u)*expo_y(u)*np.real(1/(np.sqrt(alfa_x**2 + u**2) - alfa_p ))
    term3_im = lambda u: u*expo_z(u)*expo_y(u)*np.imag(1/(np.sqrt(alfa_x**2 + u**2) - alfa_p ))
    
    term4_re = lambda u: expo_z(u)*expo_y(u)*np.real(np.sqrt(alfa_x**2 + u**2)/(np.sqrt(alfa_x**2 + u**2) - alfa_p ))
    term4_im = lambda u: expo_z(u)*expo_y(u)*np.imag(np.sqrt(alfa_x**2 + u**2)/(np.sqrt(alfa_x**2 + u**2) - alfa_p ))
    
    
    cte_final = np.exp(1j*kx*x)/((2*np.pi)**2*a)
    
    cota_inf = 0.01/omegac
    cota_sup = 600/omegac   
    
    I2_re, err = integrate.quad(term2_re, cota_inf, cota_sup)
    I3_re, err = integrate.quad(term3_re, cota_inf, cota_sup)
    I4_re, err = integrate.quad(term4_re, cota_inf, cota_sup)
    
    
    I2_im, err = integrate.quad(term2_im, cota_inf, cota_sup)
    I3_im, err = integrate.quad(term3_im, cota_inf, cota_sup)
    I4_im, err = integrate.quad(term4_im, cota_inf, cota_sup)
    
    term2_final = I2_re + 1j*I2_im
    term3_final = I3_re + 1j*I3_im
    term4_final = I4_re + 1j*I4_im
    
    return (term1 + 1j*px*kx*term2_final + 1j*py*omegac*term3_final - pz*omegac*term4_final)*cte_final


#%%
    


