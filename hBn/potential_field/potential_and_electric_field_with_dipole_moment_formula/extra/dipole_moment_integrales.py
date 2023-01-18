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
path_constants =  path_basic.replace('/potential_field/potential_and_electric_field/potential_and_electric_field_with_dipole_moment_formula','')
#print('Importar modulos necesarios para este codigo')

try:
    sys.path.insert(1, path_constants)
    from hBn_PP import hBn_lambda_p, hBn_Rp, epsilon_x, epsilon_z
except ModuleNotFoundError:
    print('graphene_sigma.py no se encuentra en ' + path_basic)


try:
    sys.path.insert(1, path_constants)
    from green_self_image import green_self_pole_aprox,green_self_ana2
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

def alpha_function(epsi1,omegac,omega0,kappa_factor_omega0,kappa_r_factor):
    """
    Parameters
    ----------
    epsilon1 : permeabilidad electrica del medio 1
    omegac : frequencia in units of 1/micrometers 
    omega0 : resonant frequency 
    kappa_factor_omega0 : kappa (collision frequency) = kappa_factor*omega0
    kappa_r_factor : kappa_r = kappa_r_factor*kappa with kappa_r_factor<1 
    Returns
    -------
    lorentzian model for polarizabilty 
    """
    omega = omegac*c
    n1 = epsi1*mu1
    cte1 = np.sqrt(n1)
    k1 = cte1
    k1_3 = k1**3
    
    kappa = kappa_factor_omega0*omega0
    kappa_r = kappa_r_factor*kappa
    A = 3*kappa_r*0.25/k1_3
    
    den = (omega0 - omega)**2 + (kappa/2)**2
    num = omega0 - omega + 1j*kappa/2

    rta = A*num/den
    return rta

#%%
    

def green_self_num_mix(omegac,epsi1,epsi3,d_nano,zp) :
    
    real_part = np.real(green_self_pole_aprox(omegac,epsi1,epsi3,d_nano,zp))
    imaginary_part = np.imag(green_self_pole_aprox(omegac,epsi1,epsi3,d_nano,zp))
    
    return real_part + 1j*imaginary_part


#%%

def INT1_ana(omegac,epsi1,epsi3,d_nano,int_v,b,zp,omega0,kappa_factor_omega0,kappa_r_factor):     
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

    alffa = alpha_function(epsi1,omegac,omega0,kappa_factor_omega0,kappa_r_factor)
    rtaself_x, rtaself_y, rtaself_z  =  green_self_num_mix(omegac,epsi1,epsi3,d_nano,zp)
    alffa_eff_x = (1/alffa -  rtaself_x)**(-1)
    alffa_eff_y = (1/alffa -  rtaself_y)**(-1)
    alffa_eff_z = (1/alffa -  rtaself_z)**(-1)

#    charge_electron = 4.806e-10/c
#    cte_uni = int_v/(2*np.pi*c)
    
    d_micros = d_nano*1e-3
    Rp = hBn_Rp(E,epsi1,epsi3)
    lambda_p_v = hBn_lambda_p(E,epsi1,epsi3)*d_micros
    kp = 2*np.pi/lambda_p_v
    alfa_p = kp/omegac

      

    
    kx = omegac*int_v
#    expo = np.exp(-np.sqrt(kx**2 + kp**2)*(np.abs(b) + 2*zp))
    expo = np.exp(-kp*(np.abs(b) + 2*zp))
    
#    den = np.sqrt(kx**2 + kp**2) - kp
    ky = np.sqrt(kp**2 - kx**2)
    kp_2 = np.sqrt(kp**2)
    term_kp = 1 + kp/kp_2
    
#    term_extra = 2*np.pi*1j*Rp*kp*np.abs(kp)*expo/ky

    px = - alffa_eff_x*2*1j*omegac*int_v*(np.pi*1j*Rp*kp*expo/ky)
    
    return px


def INT2_ana(omegac,epsi1,epsi3,d_nano,int_v,b,zp,omega0,kappa_factor_omega0,kappa_r_factor):     
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

    alffa = alpha_function(epsi1,omegac,omega0,kappa_factor_omega0,kappa_r_factor)
    rtaself_x, rtaself_y, rtaself_z  =  green_self_num_mix(omegac,epsi1,epsi3,d_nano,zp)
    alffa_eff_x = (1/alffa -  rtaself_x)**(-1)
    alffa_eff_y = (1/alffa -  rtaself_y)**(-1)
    alffa_eff_z = (1/alffa -  rtaself_z)**(-1)

#    charge_electron = 4.806e-10/c
#    cte_uni = int_v/(2*np.pi*c)
    
    d_micros = d_nano*1e-3
    Rp = hBn_Rp(E,epsi1,epsi3)
    lambda_p_v = hBn_lambda_p(E,epsi1,epsi3)*d_micros
    kp = 2*np.pi/lambda_p_v
    alfa_p = kp/omegac


    
    kx = omegac*int_v
#    expo = np.exp(-np.sqrt(kx**2 + kp**2)*(np.abs(b) + 2*zp))
    expo = np.exp(-kp*(np.abs(b) + 2*zp))
    
#    den = np.sqrt(kx**2 + kp**2) - kp
    ky = np.sqrt(kp**2 - kx**2)
    kp_2 = np.sqrt(kp**2)
    term_kp = 1 + kp/kp_2
    
#    term_extra = 2*np.pi*1j*Rp*kp*np.abs(kp)*expo/ky

    py = - alffa_eff_y*2*1j*(np.pi*1j*Rp*kp*expo)
    
    return py


def INT3_ana(omegac,epsi1,epsi3,d_nano,int_v,b,zp,omega0,kappa_factor_omega0,kappa_r_factor):     
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

    alffa = alpha_function(epsi1,omegac,omega0,kappa_factor_omega0,kappa_r_factor)
    rtaself_x, rtaself_y, rtaself_z  =  green_self_num_mix(omegac,epsi1,epsi3,d_nano,zp)
    alffa_eff_x = (1/alffa -  rtaself_x)**(-1)
    alffa_eff_y = (1/alffa -  rtaself_y)**(-1)
    alffa_eff_z = (1/alffa -  rtaself_z)**(-1)

#    charge_electron = 4.806e-10/c
#    cte_uni = int_v/(2*np.pi*c)
    
    d_micros = d_nano*1e-3
    Rp = hBn_Rp(E,epsi1,epsi3)
    lambda_p_v = hBn_lambda_p(E,epsi1,epsi3)*d_micros
    kp = 2*np.pi/lambda_p_v
    alfa_p = kp/omegac

          
    kx = omegac*int_v
#    expo = np.exp(-np.sqrt(kx**2 + kp**2)*(np.abs(b) + 2*zp))
    expo = np.exp(-kp*(np.abs(b) + 2*zp))
    
#    den = np.sqrt(kx**2 + kp**2) - kp
    ky = np.sqrt(kp**2 - kx**2)
    kp_2 = np.sqrt(kp**2)
    term_kp = 1 + kp/kp_2
    
#    term_extra = 2*np.pi*1j*Rp*kp*np.abs(kp)*expo/ky
    
    pz = alffa_eff_z*(np.pi*2*1j*Rp*(kp**2)*expo/ky)
    
    return pz

#%%

def INT1_num(omegac,epsi1,epsi3,d_nano,int_v,b,zp,omega0,kappa_factor_omega0,kappa_r_factor):     
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
    

    alffa = alpha_function(epsi1,omegac,omega0,kappa_factor_omega0,kappa_r_factor)
    rtaself_x, rtaself_y, rtaself_z  =  green_self_num_mix(omegac,epsi1,epsi3,d_nano,zp)
    alffa_eff_x = (1/alffa -  rtaself_x)**(-1)
    
#    charge_electron = alffa_eff*4.806e-10*int_v/(2*np.pi)

    
    cota_inf = 0.01/k0
    cota_sup = 600/k0   
    
    alpha_x = int_v    
    
#    term6 = np.sign(z)*pz*Rp*kp_2*J0*exp_electron
    alpha_parallel = lambda u: np.sqrt(alpha_x**2 + u**2)
    
    
#    rp_num = lambda u: epsi2*1j*kparallel(u) - epsi1*1j*kparallel(u) - cond*kparallel(u)**2/k0
#    rp_den = lambda u: epsi2*1j*kparallel(u) + epsi1*1j*kparallel(u) - cond*kparallel(u)**2/k0
#    rp = lambda u: rp_num(u)/rp_den(u)


    epsi_x = epsilon_x(E)
    epsi_z = epsilon_z(E)
    
#    if np.imag(np.sqrt(epsi_x*epsi_z))>0:
#        epsi_m = np.sqrt(epsi_x*epsi_z)
#    else:
#        epsi_m = -np.sqrt(epsi_x*epsi_z)
#    

    epsi_HBN_par = epsi_x
    epsi_HBN_perp = epsi_z

    d_micros = d_nano*1e-3

    kz1 = lambda u : np.sqrt(epsi1 - u**2) if (u**2 <= epsi1) else 1j*np.sqrt(u**2 - epsi1)
    kz2 = lambda u : np.sqrt(epsi_HBN_par - (epsi_HBN_par/epsi_HBN_perp)*u**2) if (u**2 <= epsi_HBN_perp) else 1j*np.sqrt((epsi_HBN_par/epsi_HBN_perp)*u**2 - epsi_HBN_par)
    kz3 = lambda u : np.sqrt(epsi3 - u**2) if (u**2 <= epsi3) else 1j*np.sqrt(u**2 - epsi3)

    r12 = lambda u : (kz1(u)*epsi_x - kz2(u)*epsi1)/(kz1(u)*epsi_x + kz2(u)*epsi1)
    r21 = lambda u : (kz2(u)*epsi1 - kz1(u)*epsi_x)/(kz2(u)*epsi1 + kz1(u)*epsi_x)
    r23 = lambda u : (kz2(u)*epsi3 - kz3(u)*epsi_x)/(kz2(u)*epsi3 + kz3(u)*epsi_x)
    

    exp_fresnel = lambda u: np.exp(1j*2*kz2(u)*omegac*d_micros)
    
    cte_t = np.sqrt(epsi1*epsi_x)
    t12 = lambda u : 2*kz1(u)*cte_t/(kz1(u)*epsi_x + kz2(u)*epsi1)
    t21 = lambda u : 2*kz2(u)*cte_t/(kz2(u)*epsi1 + kz1(u)*epsi_x)

    rp_num = lambda u: t12(u)*t21(u)*r23(u)*exp_fresnel(u)
    rp_den = lambda u: 1/(omegac**2) - r21(u)*r23(u)*exp_fresnel(u)
    rp = lambda u: r12(u) +  rp_num(u)/rp_den(u)   
 
      
    expo = lambda u: np.exp(-np.sqrt(alpha_x**2 + u**2)*k0*(2*zp + np.abs(b)))
    
    int_f_re_x = lambda u: np.real(rp(u)*expo(u)/alpha_parallel(u))
    int_f_im_x = lambda u: np.imag(rp(u)*expo(u)/alpha_parallel(u))
    
    INT_re_x,err = integrate.quad(int_f_re_x, cota_inf, cota_sup) 
    INT_im_x,err = integrate.quad(int_f_im_x, cota_inf, cota_sup) 
    
    INT_x = INT_re_x + 1j*INT_im_x
    
    
    px = -alffa_eff_x*1j*omegac*int_v*INT_x
    
    return px

#%%

def INT2_num(omegac,epsi1,epsi3,d_nano,int_v,b,zp,omega0,kappa_factor_omega0,kappa_r_factor):     
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
    

    alffa = alpha_function(epsi1,omegac,omega0,kappa_factor_omega0,kappa_r_factor)
    rtaself_x, rtaself_y, rtaself_z  =  green_self_num_mix(omegac,epsi1,epsi3,d_nano,zp)

    alffa_eff_y = (1/alffa -  rtaself_y)**(-1)

    
#    charge_electron = alffa_eff*4.806e-10*int_v/(2*np.pi)

    
    cota_inf = 0.01/k0
    cota_sup = 600/k0   
    
    alpha_x = int_v    
    
#    term6 = np.sign(z)*pz*Rp*kp_2*J0*exp_electron
    alpha_parallel = lambda u: np.sqrt(alpha_x**2 + u**2)
    
    
#    rp_num = lambda u: epsi2*1j*kparallel(u) - epsi1*1j*kparallel(u) - cond*kparallel(u)**2/k0
#    rp_den = lambda u: epsi2*1j*kparallel(u) + epsi1*1j*kparallel(u) - cond*kparallel(u)**2/k0
#    rp = lambda u: rp_num(u)/rp_den(u)


    epsi_x = epsilon_x(E)
    epsi_z = epsilon_z(E)
    
#    if np.imag(np.sqrt(epsi_x*epsi_z))>0:
#        epsi_m = np.sqrt(epsi_x*epsi_z)
#    else:
#        epsi_m = -np.sqrt(epsi_x*epsi_z)
#    

    epsi_HBN_par = epsi_x
    epsi_HBN_perp = epsi_z

    d_micros = d_nano*1e-3

    kz1 = lambda u : np.sqrt(epsi1 - u**2) if (u**2 <= epsi1) else 1j*np.sqrt(u**2 - epsi1)
    kz2 = lambda u : np.sqrt(epsi_HBN_par - (epsi_HBN_par/epsi_HBN_perp)*u**2) if (u**2 <= epsi_HBN_perp) else 1j*np.sqrt((epsi_HBN_par/epsi_HBN_perp)*u**2 - epsi_HBN_par)
    kz3 = lambda u : np.sqrt(epsi3 - u**2) if (u**2 <= epsi3) else 1j*np.sqrt(u**2 - epsi3)

    r12 = lambda u : (kz1(u)*epsi_x - kz2(u)*epsi1)/(kz1(u)*epsi_x + kz2(u)*epsi1)
    r21 = lambda u : (kz2(u)*epsi1 - kz1(u)*epsi_x)/(kz2(u)*epsi1 + kz1(u)*epsi_x)
    r23 = lambda u : (kz2(u)*epsi3 - kz3(u)*epsi_x)/(kz2(u)*epsi3 + kz3(u)*epsi_x)
    

    exp_fresnel = lambda u: np.exp(1j*2*kz2(u)*omegac*d_micros)
    
    cte_t = np.sqrt(epsi1*epsi_x)
    t12 = lambda u : 2*kz1(u)*cte_t/(kz1(u)*epsi_x + kz2(u)*epsi1)
    t21 = lambda u : 2*kz2(u)*cte_t/(kz2(u)*epsi1 + kz1(u)*epsi_x)

    rp_num = lambda u: t12(u)*t21(u)*r23(u)*exp_fresnel(u)
    rp_den = lambda u: 1/(omegac**2) - r21(u)*r23(u)*exp_fresnel(u)
    rp = lambda u: r12(u) +  rp_num(u)/rp_den(u)   
 
      
    expo = lambda u: np.exp(-np.sqrt(alpha_x**2 + u**2)*k0*(2*zp + np.abs(b)))
    
    
    
    int_f_re_y = lambda u: np.real(rp(alpha_parallel(u))*expo(u)*u/alpha_parallel(u))
    int_f_im_y = lambda u: np.imag(rp(alpha_parallel(u))*expo(u)*u/alpha_parallel(u))
    
    INT_re_y,err = integrate.quad(int_f_re_y, cota_inf, cota_sup) 
    INT_im_y,err = integrate.quad(int_f_im_y, cota_inf, cota_sup) 
    
    INT_y = INT_re_y + 1j*INT_im_y
    
    py = -alffa_eff_y*1j*k0*INT_y
    
    return py


#%%
    
def INT3_num(omegac,epsi1,epsi3,d_nano,int_v,b,zp,omega0,kappa_factor_omega0,kappa_r_factor):     
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
    

    alffa = alpha_function(epsi1,omegac,omega0,kappa_factor_omega0,kappa_r_factor)
    rtaself_x, rtaself_y, rtaself_z  =  green_self_num_mix(omegac,epsi1,epsi3,d_nano,zp)

    alffa_eff_z = (1/alffa -  rtaself_z)**(-1)
    
#    charge_electron = alffa_eff*4.806e-10*int_v/(2*np.pi)

    
    cota_inf = 0.01/k0
    cota_sup = 600/k0   
    
    alpha_x = int_v    
    
#    term6 = np.sign(z)*pz*Rp*kp_2*J0*exp_electron
    alpha_parallel = lambda u: np.sqrt(alpha_x**2 + u**2)
    
    
#    rp_num = lambda u: epsi2*1j*kparallel(u) - epsi1*1j*kparallel(u) - cond*kparallel(u)**2/k0
#    rp_den = lambda u: epsi2*1j*kparallel(u) + epsi1*1j*kparallel(u) - cond*kparallel(u)**2/k0
#    rp = lambda u: rp_num(u)/rp_den(u)


    epsi_x = epsilon_x(E)
    epsi_z = epsilon_z(E)
    
#    if np.imag(np.sqrt(epsi_x*epsi_z))>0:
#        epsi_m = np.sqrt(epsi_x*epsi_z)
#    else:
#        epsi_m = -np.sqrt(epsi_x*epsi_z)
#    

    epsi_HBN_par = epsi_x
    epsi_HBN_perp = epsi_z

    d_micros = d_nano*1e-3

    kz1 = lambda u : np.sqrt(epsi1 - u**2) if (u**2 <= epsi1) else 1j*np.sqrt(u**2 - epsi1)
    kz2 = lambda u : np.sqrt(epsi_HBN_par - (epsi_HBN_par/epsi_HBN_perp)*u**2) if (u**2 <= epsi_HBN_perp) else 1j*np.sqrt((epsi_HBN_par/epsi_HBN_perp)*u**2 - epsi_HBN_par)
    kz3 = lambda u : np.sqrt(epsi3 - u**2) if (u**2 <= epsi3) else 1j*np.sqrt(u**2 - epsi3)

    r12 = lambda u : (kz1(u)*epsi_x - kz2(u)*epsi1)/(kz1(u)*epsi_x + kz2(u)*epsi1)
    r21 = lambda u : (kz2(u)*epsi1 - kz1(u)*epsi_x)/(kz2(u)*epsi1 + kz1(u)*epsi_x)
    r23 = lambda u : (kz2(u)*epsi3 - kz3(u)*epsi_x)/(kz2(u)*epsi3 + kz3(u)*epsi_x)
    

    exp_fresnel = lambda u: np.exp(1j*2*kz2(u)*omegac*d_micros)
    
    cte_t = np.sqrt(epsi1*epsi_x)
    t12 = lambda u : 2*kz1(u)*cte_t/(kz1(u)*epsi_x + kz2(u)*epsi1)
    t21 = lambda u : 2*kz2(u)*cte_t/(kz2(u)*epsi1 + kz1(u)*epsi_x)

    rp_num = lambda u: t12(u)*t21(u)*r23(u)*exp_fresnel(u)
    rp_den = lambda u: 1/(omegac**2) - r21(u)*r23(u)*exp_fresnel(u)
    rp = lambda u: r12(u) +  rp_num(u)/rp_den(u)   
 
      
    expo = lambda u: np.exp(-np.sqrt(alpha_x**2 + u**2)*k0*(2*zp + np.abs(b)))
    
    int_f_re_z = lambda u: np.real(rp(alpha_parallel(u))*expo(u))
    int_f_im_z = lambda u: np.imag(rp(alpha_parallel(u))*expo(u))
    
    INT_re_z,err = integrate.quad(int_f_re_z, cota_inf, cota_sup) 
    INT_im_z,err = integrate.quad(int_f_im_z, cota_inf, cota_sup) 
    
    INT_z = INT_re_z + 1j*INT_im_z    
    
    
    pz = alffa_eff_z*k0*INT_z 
    
    return pz

#%%


def INT1_pole_aprox(omegac,epsi1,epsi3,d_nano,int_v,b,zp,omega0,kappa_factor_omega0,kappa_r_factor):     
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
    Rp = hBn_Rp(E,epsi1,epsi3)
    lambda_p_v = hBn_lambda_p(E,epsi1,epsi3)*d_micros
    kp = 2*np.pi/lambda_p_v
    alfa_p = kp/omegac

    alffa = alpha_function(epsi1,omegac,omega0,kappa_factor_omega0,kappa_r_factor)
    rtaself_x, rtaself_y, rtaself_z  =  green_self_num_mix(omegac,epsi1,epsi3,d_nano,zp) 
    alffa_eff_x = (1/alffa -  rtaself_x)**(-1)  

 
    cota_inf = 0.01/k0
    cota_sup = 600/k0   
    
    alpha_x = int_v    
    
#    term6 = np.sign(z)*pz*Rp*kp_2*J0*exp_electron
    alpha_parallel = lambda u: np.sqrt(alpha_x**2 + u**2)

    rp = lambda u: Rp*alfa_p/(alpha_parallel(u) - alfa_p)
      
    expo = lambda u: np.exp(-np.sqrt(alpha_x**2 + u**2)*k0*(2*zp + np.abs(b)))
    
    int_f_re_x = lambda u: np.real(rp(u)*expo(u)/alpha_parallel(u))
    int_f_im_x = lambda u: np.imag(rp(u)*expo(u)/alpha_parallel(u))
    
    INT_re_x,err = integrate.quad(int_f_re_x, cota_inf, cota_sup) 
    INT_im_x,err = integrate.quad(int_f_im_x, cota_inf, cota_sup) 
    
    INT_x = INT_re_x + 1j*INT_im_x
    
    px = - alffa_eff_x*1j*omegac*int_v*(INT_x)
    
    return px





def INT2_pole_aprox(omegac,epsi1,epsi3,d_nano,int_v,b,zp,omega0,kappa_factor_omega0,kappa_r_factor):     
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
    Rp = hBn_Rp(E,epsi1,epsi3)
    lambda_p_v = hBn_lambda_p(E,epsi1,epsi3)*d_micros
    kp = 2*np.pi/lambda_p_v
    alfa_p = kp/omegac

    alffa = alpha_function(epsi1,omegac,omega0,kappa_factor_omega0,kappa_r_factor)
    rtaself_x, rtaself_y, rtaself_z  =  green_self_num_mix(omegac,epsi1,epsi3,d_nano,zp) 

    alffa_eff_y = (1/alffa -  rtaself_y)**(-1)
 
    cota_inf = 0.01/k0
    cota_sup = 600/k0   
    
    alpha_x = int_v    
    
#    term6 = np.sign(z)*pz*Rp*kp_2*J0*exp_electron
    alpha_parallel = lambda u: np.sqrt(alpha_x**2 + u**2)

    rp = lambda u: Rp*alfa_p/(alpha_parallel(u) - alfa_p)
      
    expo = lambda u: np.exp(-np.sqrt(alpha_x**2 + u**2)*k0*(2*zp + np.abs(b)))
    
    
    
    int_f_re_y = lambda u: np.real(rp(u)*expo(u)*u/alpha_parallel(u))
    int_f_im_y = lambda u: np.imag(rp(u)*expo(u)*u/alpha_parallel(u))
    
    INT_re_y,err = integrate.quad(int_f_re_y, cota_inf, cota_sup) 
    INT_im_y,err = integrate.quad(int_f_im_y, cota_inf, cota_sup) 
    
    INT_y = INT_re_y + 1j*INT_im_y
        

    
    py = - alffa_eff_y*1j*(k0*INT_y)
    
    
    return py






def INT3_pole_aprox(omegac,epsi1,epsi3,d_nano,int_v,b,zp,omega0,kappa_factor_omega0,kappa_r_factor):     
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
    Rp = hBn_Rp(E,epsi1,epsi3)
    lambda_p_v = hBn_lambda_p(E,epsi1,epsi3)*d_micros
    kp = 2*np.pi/lambda_p_v
    alfa_p = kp/omegac

    alffa = alpha_function(epsi1,omegac,omega0,kappa_factor_omega0,kappa_r_factor)
    rtaself_x, rtaself_y, rtaself_z  =  green_self_num_mix(omegac,epsi1,epsi3,d_nano,zp) 
    alffa_eff_z = (1/alffa -  rtaself_z)**(-1)    

 
    cota_inf = 0.01/k0
    cota_sup = 600/k0   
    
    alpha_x = int_v    
    
#    term6 = np.sign(z)*pz*Rp*kp_2*J0*exp_electron
    alpha_parallel = lambda u: np.sqrt(alpha_x**2 + u**2)

    rp = lambda u: Rp*alfa_p/(alpha_parallel(u) - alfa_p)
      
    expo = lambda u: np.exp(-np.sqrt(alpha_x**2 + u**2)*k0*(2*zp + np.abs(b)))
    



    int_f_re_z = lambda u: np.real(rp(u)*expo(u))
    int_f_im_z = lambda u: np.imag(rp(u)*expo(u))
    
    INT_re_z,err = integrate.quad(int_f_re_z, cota_inf, cota_sup) 
    INT_im_z,err = integrate.quad(int_f_im_z, cota_inf, cota_sup) 
    
    INT_z = INT_re_z + 1j*INT_im_z


    
    
    pz = alffa_eff_z*(k0*INT_z)
    
    return pz















#%%