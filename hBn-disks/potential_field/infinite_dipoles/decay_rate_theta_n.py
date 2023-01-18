#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila
"""
import numpy as np
import sys
import os 
#from scipy import special

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/potential_field/infinite_dipoles/potential_and_electric_field_with_dipole_moment_formula/decay_rate_second_try','')
#print('Importar modulos necesarios para este codigo')


try:
    sys.path.insert(1, path_constants)
    from hBn_PP import hBn_lambda_p, hBn_Rp, epsilon_x, epsilon_z, hBn_lambda_p_Gself_image, hBn_Rp_Gself_image
except ModuleNotFoundError:
    print('hBn_PP.py no se encuentra en ' + path_constants)

try:
    sys.path.insert(1, path_basic)
    from dipole_moment import dipole_moment_sin_integrar_en_y, dipole_moment_sin_integrar_en_y_resonance
except ModuleNotFoundError:
    print('dipole_moment.py no se encuentra en ' + path_basic)
    
try:
    sys.path.insert(1, path_constants)
    from constants import constantes
except ModuleNotFoundError:
    print('constants.py no se encuentra en ' + path_constants)

pi,hb,c,alfac,mu1,mu2 = constantes()
aux = c*hb

#%%


def decay_rate_theta_inf_dipoles_ana(omegac,epsi1,epsi3,d_nano,int_v,zp,a,b,n,omega0,kappa_factor_omega0,kappa_r_factor,theta):     
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

    x, y, z = 0,0,zp 
    E = omegac*aux
    n1 = epsi1*mu1
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1
   
    d_micros = d_nano*1e-3
    Rp = hBn_Rp_Gself_image(E,epsi1,epsi3)
    lambda_p_v = hBn_lambda_p_Gself_image(E,epsi1,epsi3)*d_micros
    kp = 2*np.pi/lambda_p_v
    alfa_p = kp/omegac
    
    
    px,py,pz  = dipole_moment_sin_integrar_en_y(omegac,epsi1,epsi3,d_nano,int_v,b,zp,theta,omega0,kappa_factor_omega0,kappa_r_factor)  
#    list_dipoles = np.linspace(-Nmax,Nmax,2*Nmax + 1)
#            
    kx = omegac*int_v + 2*np.pi*n/a     
    expo_kx = np.exp(1j*kx*x)
    
    
    ky = kp*np.sin(theta)
    term_den = np.sqrt(ky**2 + kx**2)
#    den_dir = np.sqrt(np.abs(z)**2 + np.abs(y)**2)
#    K0 = special.kn(0,kx*den_dir)
#    K1 = special.kn(1,kx*den_dir)
    
    exp_electron = np.exp(-term_den*(2*zp - z))*np.exp(1j*ky*np.abs(y))*expo_kx

#    term1 =  -2*1j*px*kx*K0 + 2*py*kx*np.abs(y)*K1/den_dir  + pz*np.sign(z)*2*kx*np.abs(z)*K1/den_dir


    rp = Rp*kp/(term_den - kp)
    term2 = 1j*px*kx*rp/term_den
    
    term3 = 1j*py*ky*rp/term_den
    
    term4 = -pz*kp*rp
    
 

    final =  (term2 + term3 + term4)*exp_electron

    Ex = 1j*kx*final
    Ey = 1j*ky*final
    Ez = term_den*final

    final_2 = np.conjugate(px)*Ex + np.conjugate(py)*Ey + np.conjugate(pz)*Ez


    cte = 1/((2*np.pi)**2*a)
    
   # return np.imag(final_2*cte*kp*np.cos(theta))

    cte2 = alfac*int_v*1e15/(np.pi) ## cambio de unidades + agregar lo que faltaba en el momento dipolar
    
    return np.imag(final_2*cte*cte2)


#%%



def decay_rate_theta_inf_dipoles_ana_res(omegac,epsi1,epsi3,d_nano,int_v,zp,a,b,n,theta):     
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

    x, y, z = 0,0,0
    E = omegac*aux
    n1 = epsi1*mu1
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1
   
    d_micros = d_nano*1e-3
    Rp = hBn_Rp_Gself_image(E,epsi1,epsi3)
    lambda_p_v = hBn_lambda_p_Gself_image(E,epsi1,epsi3)*d_micros
    kp = 2*np.pi/lambda_p_v
    alfa_p = kp/omegac
    
    
    
    px,py,pz  = dipole_moment_sin_integrar_en_y_resonance(omegac,epsi1,epsi3,d_nano,int_v,b,zp,theta)  
#    list_dipoles = np.linspace(-Nmax,Nmax,2*Nmax + 1)
#            
    kx = omegac*int_v + 2*np.pi*n/a     
    expo_kx = np.exp(1j*kx*x)
    
    
    ky = kp*np.sin(theta)
    term_den = np.sqrt(ky**2 + kx**2)
#    den_dir = np.sqrt(np.abs(z)**2 + np.abs(y)**2)
#    K0 = special.kn(0,kx*den_dir)
#    K1 = special.kn(1,kx*den_dir)
    
    exp_electron = np.exp(-term_den*(2*zp - z))*np.exp(1j*ky*np.abs(y))*expo_kx

#    term1 =  -2*1j*px*kx*K0 + 2*py*kx*np.abs(y)*K1/den_dir  + pz*np.sign(z)*2*kx*np.abs(z)*K1/den_dir


    rp = Rp*kp/(term_den - kp)
    term2 = 1j*px*kx*rp/term_den
    
    term3 = 1j*py*ky*rp/term_den
    
    term4 = -pz*kp*rp
    
 

    final =  (term2 + term3 + term4)*exp_electron

    Ex = 1j*kx*final
    Ey = 1j*ky*final
    Ez = term_den*final

    final_2 = np.conjugate(px)*Ex + np.conjugate(py)*Ey + np.conjugate(pz)*Ez


    cte = 1/((2*np.pi)**2*a)
    
    cte2 = alfac*int_v*1e15/(np.pi) ## cambio de unidades + agregar lo que faltaba en el momento dipolar
    
   # return np.imag(final_2*cte*kp*np.cos(theta))

    
    return np.imag(final_2*cte*cte2)


#%%


