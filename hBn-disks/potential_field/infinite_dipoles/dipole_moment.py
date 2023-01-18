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
    from hBn_PP import hBn_lambda_p, hBn_Rp, epsilon_x, epsilon_z, hBn_lambda_p_Gself_image, hBn_Rp_Gself_image
except ModuleNotFoundError:
    print('graphene_sigma.py no se encuentra en ' + path_constants)


try:
    sys.path.insert(1, path_basic)
    from green_self_image import green_self_num, green_self_ana2
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

def polarizability_parallel(hbw,D_nano,epsilon):
    
    a_zeta, b_zeta, c_zeta = -0.01267, -45.34, 0.8635
    a_eta, b_eta, c_eta = 0.03801, -8.569, -0.1108

    D = D_nano*1e-3 ## tiene que estar en 1/micrones³
    
    t = D
    x = t/D
    
    zeta1 = a_zeta*np.exp(b_zeta*x)  + c_zeta    
    eta1 = a_eta*np.exp(b_eta*x)  + c_eta
    
#    omegac = hbw/(c*hb)
#    
    eta_parallel = 1j*np.imag(epsilon_x(hbw))/(D*epsilon)
    
    num = zeta1**2
    den = 1/eta_parallel - 1/eta1
    
    D_3 = D**3
    
    return D_3*epsilon*num/den
    

def polarizability_perp(hbw,D_nano,epsilon):
    
    a_zeta, b_zeta, c_zeta = -0.01267, -45.34, 0.8635
    a_eta, b_eta, c_eta = 0.03801, -8.569, -0.1108

    D = D_nano*1e-3 ## tiene que estar en 1/micrones³
   
    t = D
    x = t/D
    
    zeta1 = a_zeta*np.exp(b_zeta*x)  + c_zeta
    eta1 = a_eta*np.exp(b_eta*x)  + c_eta
    
#    omegac = hbw/(c*hb)
#    
    eta_perp = 1j*np.imag(epsilon_x(hbw))/(D*epsilon)
    
    num = zeta1**2
    den = 1/eta_perp - 1/eta1
    
    D_3 = D**3
    
    return D_3*epsilon*num/den


#%%

def dipole_moment_anav2(omegac,epsi1,epsi3,d_nano,int_v,b,zp,omega0,kappa_factor_omega0,kappa_r_factor):     
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
    n1 = epsi1*mu1
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1
#    k1_2 = (k0*cte1)**2
 #   n_v1 = int_v/cte1

    epsilon = 1
    D_disk_nano = 120
    alffa_x =  polarizability_parallel(E,D_disk_nano,epsilon)
    alffa_z =  polarizability_perp(E,D_disk_nano,epsilon)
    
    rtaself_x, rtaself_y, rtaself_z  =  green_self_ana2(omegac,epsi1,epsi3,d_nano,zp)
    alffa_eff_x = (1/alffa_x -  rtaself_x)**(-1)
    alffa_eff_y = (1/alffa_x -  rtaself_y)**(-1)
    alffa_eff_z = (1/alffa_z -  rtaself_z)**(-1)

    charge_electron = 4.806e-10/c
    cte_uni = int_v/(2*np.pi*c)
    

    d_micros = d_nano*1e-3
    Rp = hBn_Rp_Gself_image(E,epsi1,epsi3)
    lambda_p_v = hBn_lambda_p_Gself_image(E,epsi1,epsi3)*d_micros
    kp = 2*np.pi/lambda_p_v
    alfa_p = kp/omegac

      
    arg = np.abs(b)*omegac*int_v
    K1 = special.kn(1,arg)
    K0 = special.kn(0,arg)
    
    kx = omegac*int_v
#    expo = np.exp(-np.sqrt(kx**2 + kp**2)*(np.abs(b) + 2*zp))
    expo = np.exp(-kp*(np.abs(b) + 2*zp))
    
#    den = np.sqrt(kx**2 + kp**2) - kp
    ky = np.sqrt(kp**2 - kx**2)
    kp_2 = np.sqrt(kp**2)
    term_kp = 1 + kp/kp_2
    
#    term_extra = 2*np.pi*1j*Rp*kp*np.abs(kp)*expo/ky
    
    
    px = alffa_eff_x*1j*omegac*int_v*(K0 - 2*np.pi*1j*Rp*kp*expo/ky)
    
    py = alffa_eff_y*1j*(2*1j*omegac*int_v*K1 - 2*np.pi*1j*Rp*kp*expo)
    
    pz = alffa_eff_z*(-omegac*int_v*K1 + 2*np.pi*1j*Rp*(kp**2)*expo/ky )
    
    return px, py, pz

#%%


def dipole_moment_anav2_res(omegac,epsi1,epsi3,d_nano,int_v,b,zp,theta):     
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
    n1 = epsi1*mu1
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1
#    k1_2 = (k0*cte1)**2
 #   n_v1 = int_v/cte1


    rtaself_x, rtaself_y, rtaself_z  = green_self_ana2(omegac,epsi1,epsi3,d_nano,zp)
    alffa_eff_x = (-2*omegac**3/(3*epsi1) -  np.imag(rtaself_x))**(-1)
    alffa_eff_y = (-2*omegac**3/(3*epsi1) -  np.imag(rtaself_y))**(-1)
    alffa_eff_z = (-2*omegac**3/(3*epsi1) -  np.imag(rtaself_z))**(-1)

    charge_electron = 4.806e-10/c
    cte_uni = int_v/(2*np.pi*c)
    

    d_micros = d_nano*1e-3
    Rp = hBn_Rp_Gself_image(E,epsi1,epsi3)
    lambda_p_v = hBn_lambda_p_Gself_image(E,epsi1,epsi3)*d_micros
    kp = 2*np.pi/lambda_p_v
    alfa_p = kp/omegac
      
    arg = np.abs(b)*omegac*int_v
    K1 = special.kn(1,arg)
    K0 = special.kn(0,arg)
    
    kx = omegac*int_v
#    expo = np.exp(-np.sqrt(kx**2 + kp**2)*(np.abs(b) + 2*zp))
    expo = np.exp(-kp*(np.abs(b) + 2*zp))
    
#    den = np.sqrt(kx**2 + kp**2) - kp
    ky = np.sqrt(kp**2 - kx**2)
    kp_2 = np.sqrt(kp**2)
    term_kp = 1 + kp/kp_2
    
#    term_extra = 2*np.pi*1j*Rp*kp*np.abs(kp)*expo/ky
    
    
    px = alffa_eff_x*1j*omegac*int_v*(K0 - 2*np.pi*1j*Rp*kp*expo/ky)
    
    py = alffa_eff_y*1j*(2*1j*omegac*int_v*K1 - 2*np.pi*1j*Rp*kp*expo)
    
    pz = alffa_eff_z*(-omegac*int_v*K1 + 2*np.pi*1j*Rp*(kp**2)*expo/ky )
    
    return px, py, pz

#%%
    



def dipole_moment_sin_integrar_en_y(omegac,epsi1,epsi3,d_nano,int_v,b,zp,theta,omega0,kappa_factor_omega0,kappa_r_factor):     
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
    k0 = omegac #=omega/c
    # x_y = ky/k0
#    n1 = epsi1*mu1
#    cte1 = np.sqrt(n1)
    # k1 = omegac*cte1
#    k1_2 = (k0*cte1)**2
 #   n_v1 = int_v/cte1
    
    d_micros = d_nano*1e-3
    Rp = hBn_Rp_Gself_image(E,epsi1,epsi3)
    lambda_p_v = hBn_lambda_p_Gself_image(E,epsi1,epsi3)*d_micros
    kp = 2*np.pi/lambda_p_v
    alfa_p = kp/omegac

    epsilon = 1
    D_disk_nano = 120
    alffa_x =  polarizability_parallel(E,D_disk_nano,epsilon)
    alffa_z =  polarizability_perp(E,D_disk_nano,epsilon)
    
    rtaself_x, rtaself_y, rtaself_z  =  green_self_ana2(omegac,epsi1,epsi3,d_nano,zp)
    alffa_eff_x = (1/alffa_x -  rtaself_x)**(-1)
    alffa_eff_y = (1/alffa_x -  rtaself_y)**(-1)
    alffa_eff_z = (1/alffa_z -  rtaself_z)**(-1)    

 
    alpha_x = int_v    
    alpha_y = alfa_p*np.sin(theta)
    
#    term6 = np.sign(z)*pz*Rp*kp_2*J0*exp_electron
    alpha_parallel = np.sqrt(alpha_x**2 + alpha_y**2)

    rp = Rp*alfa_p/(alpha_parallel - alfa_p)
      
    expo_2 =  np.exp(-np.sqrt(alpha_x**2 + alpha_y**2)*k0*(2*zp + np.abs(b)))
    expo_1 =  np.exp(-np.sqrt(alpha_x**2 + alpha_y**2)*k0*np.abs(b))
    
    INT_x = (expo_1 - rp*expo_2)/alpha_parallel
    
    INT_y = (expo_1 - rp*expo_2)*alpha_y/alpha_parallel
            
    INT_z = -expo_1 + rp*expo_2


    
    px = alffa_eff_x*1j*omegac*alpha_x*INT_x
    
    py = alffa_eff_y*1j*omegac*INT_y
    
    pz = alffa_eff_z*omegac*INT_z
    
    return px*alfa_p*np.cos(theta) , py*alfa_p*np.cos(theta) , pz*alfa_p*np.cos(theta)

#%%


def dipole_moment_sin_integrar_en_y_resonance(omegac,epsi1,epsi3,d_nano,int_v,b,zp,theta):     
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
    k0 = omegac #=omega/c
    # x_y = ky/k0
#    n1 = epsi1*mu1
#    cte1 = np.sqrt(n1)
    # k1 = omegac*cte1
#    k1_2 = (k0*cte1)**2
 #   n_v1 = int_v/cte1
    
    d_micros = d_nano*1e-3
    Rp = hBn_Rp_Gself_image(E,epsi1,epsi3)
    lambda_p_v = hBn_lambda_p_Gself_image(E,epsi1,epsi3)*d_micros
    kp = 2*np.pi/lambda_p_v
    alfa_p = kp/omegac
#    kp = alfa_p*k0

    rtaself_x, rtaself_y, rtaself_z  = green_self_ana2(omegac,epsi1,epsi3,d_nano,zp)
    alffa_eff_x = (-2*omegac**3/(3*epsi1) -  np.imag(rtaself_x))**(-1)
    alffa_eff_y = (-2*omegac**3/(3*epsi1) -  np.imag(rtaself_y))**(-1)
    alffa_eff_z = (-2*omegac**3/(3*epsi1) -  np.imag(rtaself_z))**(-1)

    
    alpha_x = int_v    
    alpha_y = alfa_p*np.sin(theta)
    
#    term6 = np.sign(z)*pz*Rp*kp_2*J0*exp_electron
    alpha_parallel = np.sqrt(alpha_x**2 + alpha_y**2)

    rp = Rp*alfa_p/(alpha_parallel - alfa_p)
      
    expo_2 =  np.exp(-np.sqrt(alpha_x**2 + alpha_y**2)*k0*(2*zp + np.abs(b)))
    expo_1 =  np.exp(-np.sqrt(alpha_x**2 + alpha_y**2)*k0*np.abs(b))
    
    INT_x = (expo_1 - rp*expo_2)/alpha_parallel
    
    INT_y = (expo_1 - rp*expo_2)*alpha_y/alpha_parallel
            
    INT_z = -expo_1 + rp*expo_2


    
    px = alffa_eff_x*1j*omegac*alpha_x*INT_x
    
    py = alffa_eff_y*1j*omegac*INT_y
    
    pz = alffa_eff_z*omegac*INT_z
    
    return px*alfa_p*np.cos(theta) , py*alfa_p*np.cos(theta) , pz*alfa_p*np.cos(theta)


#%%
    
