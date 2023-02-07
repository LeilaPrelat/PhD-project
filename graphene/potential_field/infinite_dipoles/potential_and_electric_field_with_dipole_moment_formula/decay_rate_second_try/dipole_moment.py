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
path_constants =  path_basic.replace('/potential_field/infinite_dipoles/potential_and_electric_field_with_dipole_moment_formula/decay_rate_second_try','')
#print('Importar modulos necesarios para este codigo')

try:
    sys.path.insert(1, path_constants)
    from graphene_sigma import sigma_DL
except ModuleNotFoundError:
    print('graphene_sigma.py no se encuentra en ' + path_basic)


try:
    sys.path.insert(1, path_basic)
    from green_self_image import green_self_num,green_self_ana2,green_self_num_integral_inside_light_cone
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

def dipole_moment_anav2(omegac,epsi1,epsi2,hbmu,hbgama,int_v,b,zp,omega0,kappa_factor_omega0,kappa_r_factor):     
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
    rtaself_x, rtaself_y, rtaself_z  =  green_self_num(omegac,epsi1,epsi2,hbmu,hbgama,zp)
    alffa_eff_x = (1/alffa -  rtaself_x)**(-1)
    alffa_eff_y = (1/alffa -  rtaself_y)**(-1)
    alffa_eff_z = (1/alffa -  rtaself_z)**(-1)

    charge_electron = 4.806e-10/c
    cte_uni = int_v/(2*np.pi*c)
    

    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    Rp = 2*epsi1/(epsi1 + epsi2)
    alfa_p = 1j*(epsi1 + epsi2)/(cond*cte1)
    kp = alfa_p*k1

      
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


def dipole_moment_anav2_res(omegac,epsi1,epsi2,hbmu,hbgama,int_v,b,zp):     
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

    rtaself_x, rtaself_y, rtaself_z  = green_self_ana2(omegac,epsi1,epsi2,hbmu,hbgama,zp)
    alffa_eff_x = 1j*(2*omegac**3/(3*epsi1) +  np.imag(rtaself_x))**(-1)
    alffa_eff_y = 1j*(2*omegac**3/(3*epsi1) +  np.imag(rtaself_y))**(-1)
    alffa_eff_z = 1j*(2*omegac**3/(3*epsi1) +  np.imag(rtaself_z))**(-1)
    

#    alffa_eff_x = 1
#    alffa_eff_y = 1
#    alffa_eff_z = 1
##    
    

    charge_electron = 4.806e-10/c
    cte_uni = int_v/(2*np.pi*c)
    

    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    Rp = 2*epsi1/(epsi1 + epsi2)
    alfa_p = 1j*(epsi1 + epsi2)/(cond*cte1)
    kp = alfa_p*k1

      
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



def dipole_moment_num_res(omegac,epsi1,epsi2,hbmu,hbgama,int_v,b,zp):     
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
    
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama) #no es necesario dividir por c 



    rtaself_x, rtaself_y, rtaself_z  =  green_self_num(omegac,epsi1,epsi2,hbmu,hbgama,zp)
    alffa_eff_x = 1j*(2*omegac**3/(3*epsi1) +  np.imag(rtaself_x))**(-1)
    alffa_eff_y = 1j*(2*omegac**3/(3*epsi1) +  np.imag(rtaself_y))**(-1)
    alffa_eff_z = 1j*(2*omegac**3/(3*epsi1) +  np.imag(rtaself_z))**(-1)
    
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


    rp_num = lambda u: epsi2*1j*alpha_parallel(u) - epsi1*1j*alpha_parallel(u) - cond*(alpha_parallel(u)**2)
    rp_den = lambda u: epsi2*1j*alpha_parallel(u) + epsi1*1j*alpha_parallel(u) - cond*(alpha_parallel(u)**2)
    rp = lambda u: rp_num(u)/rp_den(u)   
 
      
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
    
    return px, py, pz



def dipole_moment_anav2_for_decay_rate_res(omegac,epsi1,epsi2,hbmu,hbgama,int_v,b,zp):     
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

    rtaself_x, rtaself_y, rtaself_z  = green_self_ana2(omegac,epsi1,epsi2,hbmu,hbgama,zp)
    alffa_eff_x = 1j*(2*omegac**3/(3*epsi1) +  np.imag(rtaself_x))**(-1)
    alffa_eff_y = 1j*(2*omegac**3/(3*epsi1) +  np.imag(rtaself_y))**(-1)
    alffa_eff_z = 1j*(2*omegac**3/(3*epsi1) +  np.imag(rtaself_z))**(-1)
    

#    alffa_eff_x = 1
#    alffa_eff_y = 1
#    alffa_eff_z = 1
##    
    

    charge_electron = 4.806e-10/c
    cte_uni = int_v/(2*np.pi*c)
    

    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    Rp = 2*epsi1/(epsi1 + epsi2)
    alfa_p = 1j*(epsi1 + epsi2)/(cond*cte1)
    kp = alfa_p*k1

      
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





def dipole_moment_anav2_for_decay_rate_resonance_dir(omegac,int_v,b,zp):     
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

    arg = np.abs(b)*omegac*int_v
    K1 = special.kn(1,arg)
    K0 = special.kn(0,arg)

    px = 1j*omegac*int_v*K0
    
    py = -2*omegac*int_v*K1
    
    pz = -omegac*int_v*K1
    
    return px, py, pz







#%%
    



def dipole_moment_sin_integrar_en_y(omegac,epsi1,epsi2,hbmu,hbgama,int_v,b,zp,theta,omega0,kappa_factor_omega0,kappa_r_factor):     
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
    
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    Rp = 2*epsi1/(epsi1 + epsi2)
    alfa_p = 1j*(epsi1 + epsi2)/(cond)
    kp = alfa_p*k0

    alffa = alpha_function(epsi1,omegac,omega0,kappa_factor_omega0,kappa_r_factor)
    rtaself_x, rtaself_y, rtaself_z  =  green_self_ana2(omegac,epsi1,epsi2,hbmu,hbgama,zp) 
    alffa_eff_x = (1/alffa -  rtaself_x)**(-1)
    alffa_eff_y = (1/alffa -  rtaself_y)**(-1)
    alffa_eff_z = (1/alffa -  rtaself_z)**(-1)    

 
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


def dipole_moment_sin_integrar_en_y_resonance(omegac,epsi1,epsi2,hbmu,hbgama,int_v,b,zp,theta):     
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
    
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    Rp = 2*epsi1/(epsi1 + epsi2)
    alfa_p = 1j*(epsi1 + epsi2)/(cond)
#    kp = alfa_p*k0

    rtaself_x, rtaself_y, rtaself_z  = green_self_ana2(omegac,epsi1,epsi2,hbmu,hbgama,zp)
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
    
