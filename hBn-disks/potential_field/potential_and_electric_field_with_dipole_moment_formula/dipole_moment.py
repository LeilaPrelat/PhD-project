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
path_constants =  path_basic.replace('/potential_field/potential_and_electric_field_with_dipole_moment_formula','')
#print('Importar modulos necesarios para este codigo')

try:
    sys.path.insert(1, path_constants)
    from hBn_PP import hBn_lambda_p, hBn_Rp, epsilon_x, epsilon_z, hBn_lambda_p_Gself_image, hBn_Rp_Gself_image
except ModuleNotFoundError:
    print('graphene_sigma.py no se encuentra en ' + path_constants)

try:
    sys.path.insert(1, path_constants)
    from green_self_image import green_self_pole_aprox, green_self_ana2, green_self_num
except ModuleNotFoundError:
    print('graphene_sigma.py no se encuentra en ' + path_constants)


try:
    sys.path.insert(1, path_constants)
    from constants import constantes
except ModuleNotFoundError:
    print('constants.py no se encuentra en ' + path_constants)

pi,hb,c,alfac,mu1,mu2 = constantes()
aux = c*hb


#%%
#
#def alpha_function(epsi1,omegac,omega0,kappa_factor_omega0,kappa_r_factor):
#    """
#    Parameters
#    ----------
#    epsilon1 : permeabilidad electrica del medio 1
#    omegac : frequencia in units of 1/micrometers 
#    omega0 : resonant frequency 
#    kappa_factor_omega0 : kappa (collision frequency) = kappa_factor*omega0
#    kappa_r_factor : kappa_r = kappa_r_factor*kappa with kappa_r_factor<1 
#    Returns
#    -------
#    lorentzian model for polarizabilty 
#    """
#    omega = omegac*c
#    n1 = epsi1*mu1
#    cte1 = np.sqrt(n1)
#    k1 = cte1
#    k1_3 = k1**3
#    
#    kappa = kappa_factor_omega0*omega0
#    kappa_r = kappa_r_factor*kappa
#    A = 3*kappa_r*0.25/k1_3
#    
#    den = (omega0 - omega)**2 + (kappa/2)**2
#    num = omega0 - omega + 1j*kappa/2
#
#    rta = A*num/den
#    return rta

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
    

def green_self_num_mix(omegac,epsi1,epsi3,d_nano,zp) :
    
    real_part = np.real(green_self_num(omegac,epsi1,epsi3,d_nano,zp))
    imaginary_part = np.imag(green_self_num(omegac,epsi1,epsi3,d_nano,zp))
    
    return real_part + 1j*imaginary_part


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

#    charge_electron = 4.806e-10/c
#    cte_uni = int_v/(2*np.pi*c)
    
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
    
    
    px = alffa_eff_x*1j*omegac*int_v*(K0 - np.pi*1j*Rp*kp*expo/ky)
    
    py = alffa_eff_y*1j*(2*1j*omegac*int_v*K1 - np.pi*1j*Rp*kp*expo)
    
    pz = alffa_eff_z*(-omegac*int_v*K1 + np.pi*1j*Rp*(kp**2)*expo/ky )
    
    return px, py, pz









def dipole_moment_anav2_resonance(omegac,epsi1,epsi3,d_nano,int_v,b,zp):     
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

    rtaself_x, rtaself_y, rtaself_z  =  green_self_ana2(omegac,epsi1,epsi3,d_nano,zp)
    alffa_eff_x = (-2*omegac**3/(3*epsi1) -  np.imag(rtaself_x))**(-1)
    alffa_eff_y = (-2*omegac**3/(3*epsi1) -  np.imag(rtaself_y))**(-1)
    alffa_eff_z = (-2*omegac**3/(3*epsi1) -  np.imag(rtaself_z))**(-1)

#    charge_electron = 4.806e-10/c
#    cte_uni = int_v/(2*np.pi*c)
    
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
    
    
    px = alffa_eff_x*1j*omegac*int_v*(K0 - np.pi*1j*Rp*kp*expo/ky)
    
    py = alffa_eff_y*1j*(2*1j*omegac*int_v*K1 - np.pi*1j*Rp*kp*expo)
    
    pz = alffa_eff_z*(-omegac*int_v*K1 + np.pi*1j*Rp*(kp**2)*expo/ky )
    
    return px, py, pz


#%%

def dipole_moment_num(omegac,epsi1,epsi3,d_nano,int_v,b,zp,omega0,kappa_factor_omega0,kappa_r_factor):     
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
    

    epsilon = 1
    D_disk_nano = 120
    alffa_x =  polarizability_parallel(E,D_disk_nano,epsilon)
    alffa_z =  polarizability_perp(E,D_disk_nano,epsilon)
    
    rtaself_x, rtaself_y, rtaself_z  =  green_self_num(omegac,epsi1,epsi3,d_nano,zp)
    alffa_eff_x = (1/alffa_x -  rtaself_x)**(-1)
    alffa_eff_y = (1/alffa_x -  rtaself_y)**(-1)
    alffa_eff_z = (1/alffa_z -  rtaself_z)**(-1)
    
#    charge_electron = alffa_eff*4.806e-10*int_v/(2*np.pi)
    charge_electron = 4.806e-10/c
    cte_uni = int_v/(2*np.pi*c)
    
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
    kz2 = lambda u : np.sqrt(epsi_HBN_par - (epsi_HBN_par/epsi_HBN_perp)*u**2) 
    kz3 = lambda u : np.sqrt(epsi3 - u**2) if (u**2 <= epsi3) else 1j*np.sqrt(u**2 - epsi3)

#    if np.imag(kz1) > 0 :
#        kz1_v2 = lambda u: kz1(u)
#    else:
#        kz1_v2 = lambda u: -kz1(u)
#        
#    if np.imag(kz2) > 0 :
#        kz2_v2 = lambda u: kz2(u)
#    else:
#        kz2_v2 = lambda u: -kz2(u)
#    
#
#    if np.imag(kz3) > 0 :
#        kz3_v2 = lambda u: kz3(u)
#    else:
#        kz3_v2 = lambda u: -kz3(u)

    r12 = lambda u : (kz1(u)*epsi_x - kz2(u)*epsi1)/(kz1(u)*epsi_x + kz2(u)*epsi1)
    r21 = lambda u : (kz2(u)*epsi1 - kz1(u)*epsi_x)/(kz2(u)*epsi1 + kz1(u)*epsi_x)
    r23 = lambda u : (kz2(u)*epsi3 - kz3(u)*epsi_x)/(kz2(u)*epsi3 + kz3(u)*epsi_x)
    

    exp_fresnel = lambda u: np.exp(1j*2*kz2(u)*omegac*d_micros)
    
    cte_t = np.sqrt(epsi1*epsi_x)
    t12 = lambda u : 2*kz1(u)*cte_t/(kz1(u)*epsi_x + kz2(u)*epsi1)
    t21 = lambda u : 2*kz2(u)*cte_t/(kz2(u)*epsi1 + kz1(u)*epsi_x)

    rp_num = lambda u: t12(u)*t21(u)*r23(u)*exp_fresnel(u)
    rp_den = lambda u: 1 - r21(u)*r23(u)*exp_fresnel(u)
    rp = lambda u: r12(u) +  rp_num(u)/rp_den(u)   
 
      
    expo = lambda u: np.exp(-np.sqrt(alpha_x**2 + u**2)*k0*(2*zp + np.abs(b)))
    
    int_f_re_x = lambda u: np.real(rp(alpha_parallel(u))*expo(u)/alpha_parallel(u))
    int_f_im_x = lambda u: np.imag(rp(alpha_parallel(u))*expo(u)/alpha_parallel(u))
    
    INT_re_x,err = integrate.quad(int_f_re_x, cota_inf, cota_sup) 
    INT_im_x,err = integrate.quad(int_f_im_x, cota_inf, cota_sup) 
    
    INT_x = INT_re_x + 1j*INT_im_x
    
    
    int_f_re_y = lambda u: np.real(rp(alpha_parallel(u))*expo(u)*u/alpha_parallel(u))
    int_f_im_y = lambda u: np.imag(rp(alpha_parallel(u))*expo(u)*u/alpha_parallel(u))
    
    INT_re_y,err = integrate.quad(int_f_re_y, cota_inf, cota_sup) 
    INT_im_y,err = integrate.quad(int_f_im_y, cota_inf, cota_sup) 
    
    INT_y = INT_re_y + 1j*INT_im_y
        


    int_f_re_z = lambda u: np.real(rp(alpha_parallel(u))*expo(u))
    int_f_im_z = lambda u: np.imag(rp(alpha_parallel(u))*expo(u))
    
    INT_re_z,err = integrate.quad(int_f_re_z, cota_inf, cota_sup) 
    INT_im_z,err = integrate.quad(int_f_im_z, cota_inf, cota_sup) 
    
    INT_z = INT_re_z + 1j*INT_im_z


    arg = np.abs(b)*omegac*int_v
    K1 = special.kn(1,arg)
    K0 = special.kn(0,arg)    
    
    
    px = alffa_eff_x*1j*omegac*int_v*(K0 - INT_x)
    
    py = alffa_eff_y*1j*(2*1j*omegac*int_v*K1 - k0*INT_y)
    
    pz = alffa_eff_z*(-omegac*int_v*K1 + k0*INT_z )
    
    return px, py, pz









def dipole_moment_num_resonance(omegac,epsi1,epsi3,d_nano,int_v,b,zp):     
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
    

    rtaself_x, rtaself_y, rtaself_z  =  green_self_num(omegac,epsi1,epsi3,d_nano,zp)
    alffa_eff_x = (-2*omegac**3/(3*epsi1) -  np.imag(rtaself_x))**(-1)
    alffa_eff_y = (-2*omegac**3/(3*epsi1) -  np.imag(rtaself_y))**(-1)
    alffa_eff_z = (-2*omegac**3/(3*epsi1) -  np.imag(rtaself_z))**(-1)
    
#    charge_electron = alffa_eff*4.806e-10*int_v/(2*np.pi)
    charge_electron = 4.806e-10/c
    cte_uni = int_v/(2*np.pi*c)
    
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
    kz2 = lambda u : np.sqrt(epsi_HBN_par - (epsi_HBN_par/epsi_HBN_perp)*u**2) 
    kz3 = lambda u : np.sqrt(epsi3 - u**2) if (u**2 <= epsi3) else 1j*np.sqrt(u**2 - epsi3)

#    if np.imag(kz1) > 0 :
#        kz1_v2 = lambda u: kz1(u)
#    else:
#        kz1_v2 = lambda u: -kz1(u)
#        
#    if np.imag(kz2) > 0 :
#        kz2_v2 = lambda u: kz2(u)
#    else:
#        kz2_v2 = lambda u: -kz2(u)
#    
#
#    if np.imag(kz3) > 0 :
#        kz3_v2 = lambda u: kz3(u)
#    else:
#        kz3_v2 = lambda u: -kz3(u)

    r12 = lambda u : (kz1(u)*epsi_x - kz2(u)*epsi1)/(kz1(u)*epsi_x + kz2(u)*epsi1)
    r21 = lambda u : (kz2(u)*epsi1 - kz1(u)*epsi_x)/(kz2(u)*epsi1 + kz1(u)*epsi_x)
    r23 = lambda u : (kz2(u)*epsi3 - kz3(u)*epsi_x)/(kz2(u)*epsi3 + kz3(u)*epsi_x)
    

    exp_fresnel = lambda u: np.exp(1j*2*kz2(u)*omegac*d_micros)
    
    cte_t = np.sqrt(epsi1*epsi_x)
    t12 = lambda u : 2*kz1(u)*cte_t/(kz1(u)*epsi_x + kz2(u)*epsi1)
    t21 = lambda u : 2*kz2(u)*cte_t/(kz2(u)*epsi1 + kz1(u)*epsi_x)

    rp_num = lambda u: t12(u)*t21(u)*r23(u)*exp_fresnel(u)
    rp_den = lambda u: 1 - r21(u)*r23(u)*exp_fresnel(u)
    rp = lambda u: r12(u) +  rp_num(u)/rp_den(u)   
      
    expo = lambda u: np.exp(-np.sqrt(alpha_x**2 + u**2)*k0*(2*zp + np.abs(b)))
    
    int_f_re_x = lambda u: np.real(rp(alpha_parallel(u))*expo(u)/alpha_parallel(u))
    int_f_im_x = lambda u: np.imag(rp(alpha_parallel(u))*expo(u)/alpha_parallel(u))
    
    INT_re_x,err = integrate.quad(int_f_re_x, cota_inf, cota_sup) 
    INT_im_x,err = integrate.quad(int_f_im_x, cota_inf, cota_sup) 
    
    INT_x = INT_re_x + 1j*INT_im_x
    
    
    int_f_re_y = lambda u: np.real(rp(alpha_parallel(u))*expo(u)*u/alpha_parallel(u))
    int_f_im_y = lambda u: np.imag(rp(alpha_parallel(u))*expo(u)*u/alpha_parallel(u))
    
    INT_re_y,err = integrate.quad(int_f_re_y, cota_inf, cota_sup) 
    INT_im_y,err = integrate.quad(int_f_im_y, cota_inf, cota_sup) 
    
    INT_y = INT_re_y + 1j*INT_im_y
        


    int_f_re_z = lambda u: np.real(rp(alpha_parallel(u))*expo(u))
    int_f_im_z = lambda u: np.imag(rp(alpha_parallel(u))*expo(u))
    
    INT_re_z,err = integrate.quad(int_f_re_z, cota_inf, cota_sup) 
    INT_im_z,err = integrate.quad(int_f_im_z, cota_inf, cota_sup) 
    
    INT_z = INT_re_z + 1j*INT_im_z


    arg = np.abs(b)*omegac*int_v
    K1 = special.kn(1,arg)
    K0 = special.kn(0,arg)    
    
    
    px = alffa_eff_x*1j*omegac*int_v*(K0 - INT_x)
    
    py = alffa_eff_y*1j*(2*1j*omegac*int_v*K1 - k0*INT_y)
    
    pz = alffa_eff_z*(-omegac*int_v*K1 + k0*INT_z )
    
    return px, py, pz



#%%


def dipole_moment_pole_aprox(omegac,epsi1,epsi3,d_nano,int_v,b,zp,omega0,kappa_factor_omega0,kappa_r_factor):     
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
    rtaself_x, rtaself_y, rtaself_z  =  green_self_pole_aprox(omegac,epsi1,epsi3,d_nano,zp) 
    alffa_eff_x = (1/alffa_x -  rtaself_x)**(-1)
    alffa_eff_y = (1/alffa_x -  rtaself_y)**(-1)
    alffa_eff_z = (1/alffa_z -  rtaself_z)**(-1)    

 
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
    
    
    int_f_re_y = lambda u: np.real(rp(u)*expo(u)*u/alpha_parallel(u))
    int_f_im_y = lambda u: np.imag(rp(u)*expo(u)*u/alpha_parallel(u))
    
    INT_re_y,err = integrate.quad(int_f_re_y, cota_inf, cota_sup) 
    INT_im_y,err = integrate.quad(int_f_im_y, cota_inf, cota_sup) 
    
    INT_y = INT_re_y + 1j*INT_im_y
        


    int_f_re_z = lambda u: np.real(rp(u)*expo(u))
    int_f_im_z = lambda u: np.imag(rp(u)*expo(u))
    
    INT_re_z,err = integrate.quad(int_f_re_z, cota_inf, cota_sup) 
    INT_im_z,err = integrate.quad(int_f_im_z, cota_inf, cota_sup) 
    
    INT_z = INT_re_z + 1j*INT_im_z


    arg = np.abs(b)*omegac*int_v
    K1 = special.kn(1,arg)
    K0 = special.kn(0,arg)    
    
    
    px = alffa_eff_x*1j*omegac*int_v*(K0 - INT_x)
    
    py = alffa_eff_y*1j*(2*1j*omegac*int_v*K1 - k0*INT_y)
    
    pz = alffa_eff_z*(-omegac*int_v*K1 + k0*INT_z )
    
    return px,py,pz

#%%
    



def dipole_moment_pole_aprox_resonance(omegac,epsi1,epsi3,d_nano,int_v,b,zp):     
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

    rtaself_x, rtaself_y, rtaself_z  =  green_self_pole_aprox(omegac,epsi1,epsi3,d_nano,zp)
    alffa_eff_x = (-2*omegac**3/(3*epsi1) -  np.imag(rtaself_x))**(-1)
    alffa_eff_y = (-2*omegac**3/(3*epsi1) -  np.imag(rtaself_y))**(-1)
    alffa_eff_z = (-2*omegac**3/(3*epsi1) -  np.imag(rtaself_z))**(-1)
 
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
    
    
    int_f_re_y = lambda u: np.real(rp(u)*expo(u)*u/alpha_parallel(u))
    int_f_im_y = lambda u: np.imag(rp(u)*expo(u)*u/alpha_parallel(u))
    
    INT_re_y,err = integrate.quad(int_f_re_y, cota_inf, cota_sup) 
    INT_im_y,err = integrate.quad(int_f_im_y, cota_inf, cota_sup) 
    
    INT_y = INT_re_y + 1j*INT_im_y
        


    int_f_re_z = lambda u: np.real(rp(u)*expo(u))
    int_f_im_z = lambda u: np.imag(rp(u)*expo(u))
    
    INT_re_z,err = integrate.quad(int_f_re_z, cota_inf, cota_sup) 
    INT_im_z,err = integrate.quad(int_f_im_z, cota_inf, cota_sup) 
    
    INT_z = INT_re_z + 1j*INT_im_z


    arg = np.abs(b)*omegac*int_v
    K1 = special.kn(1,arg)
    K0 = special.kn(0,arg)    
    
    
    px = alffa_eff_x*1j*omegac*int_v*(K0 - INT_x)
    
    py = alffa_eff_y*1j*(2*1j*omegac*int_v*K1 - k0*INT_y)
    
    pz = alffa_eff_z*(-omegac*int_v*K1 + k0*INT_z )
    
    return px,py,pz






