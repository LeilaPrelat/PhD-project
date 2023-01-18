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


def alpha_function_num(omegac,epsi1,epsi2,hbmu,hbgama,px,py,pz,zp,omega0,kappa_factor_omega0,kappa_r_factor):
    """
    Parameters
    ----------
    omegac : omega/c = k0 en 1/micrometros    
    epsi1 : epsilon del medio de arriba del plano
    epsi2 : epsilon del medio de abajo del plano
    hbmu : chemical potential in eV  
    hbgama : collision frequency in eV
    int_v : fraccion de la vel de la c que vale la velocidad uniforme del electron en x, si v = c/int
    b: electron position en z, b < 0 (dipolo en 0 y eje de z apunta hacia abajo) 
    zp : coordenada zp del plano, zp > 0  (dipolo en 0 y eje de z apunta hacia abajo) 
    x : punto x donde miramos el campo en micrones
    omega0 : resonant frequency 
    kappa_factor_omega0 : kappa (collision frequency) = kappa_factor*omega0
    kappa_r_factor : kappa_r = kappa_r_factor*kappa with kappa_r_factor<1    
    Returns
    -------
    alfa effectivo en QE approx
    """

    E = omegac*aux  
    n1 = epsi1*mu1
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1
    k1_3 = k1**3

    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)

    rp_num = lambda u: epsi2*1j*u - epsi1*1j*u - cte1*cond*(u**2)
    rp_den = lambda u: epsi2*1j*u + epsi1*1j*u - cte1*cond*(u**2)
    rp = lambda u: rp_num(u)/rp_den(u)
    
    rs_num = lambda u: 1j*u - 1j*u - cond/cte1
    rs_den = lambda u: 1j*u + 1j*u + cond/cte1
    rs = lambda u: rs_num(u)/rs_den(u)

    cota_sup1 = 1000*k1
####       
    cte_x = k1_3*0.5 #signo menos

    alfa = alpha_function(epsi1,omegac,omega0,kappa_factor_omega0,kappa_r_factor)
    
    z_dip_barra_self = k1*2*zp  #NO ESTOY SEGURA DE ESTA EXPONENCIAL EH 
    expB_self = lambda u: np.exp(-u*z_dip_barra_self) 

    IntselfB_function_re_x = lambda u: np.real(((u**2)*rp(u) + rs(u))*expB_self(u))
    IntselfB_function_im_x = lambda u: np.imag(((u**2)*rp(u) + rs(u))*expB_self(u))

    intselfB_re_x,err = integrate.quad(IntselfB_function_re_x, 0.1, cota_sup1)
    intselfB_im_x,err = integrate.quad(IntselfB_function_im_x, 0.1, cota_sup1)
    
    
    rtaself_x = (py + px)*(intselfB_re_x + 1j*intselfB_im_x)*cte_x
    

    IntselfB_function_re_z = lambda u: np.real((u**2)*rp(u)*expB_self(u))
    IntselfB_function_im_z = lambda u: np.imag((u**2)*rp(u)*expB_self(u))

    intselfB_re_z,err = integrate.quad(IntselfB_function_re_z, 0.1, cota_sup1)
    intselfB_im_z,err = integrate.quad(IntselfB_function_im_z, 0.1, cota_sup1)    
    
    rtaself_z = pz*(intselfB_re_z + 1j*intselfB_im_z)*cte_x*2
    
    
    rtaself = rtaself_x + rtaself_z
    
    rta_final = (1/alfa - rtaself)**(-1)

    return rta_final

#%%

def alpha_function_ana(omegac,epsi1,epsi2,hbmu,hbgama,px,py,pz,zp,omega0,kappa_factor_omega0,kappa_r_factor):
    """
    Parameters
    ----------
    omegac : omega/c = k0 en 1/micrometros    
    epsi1 : epsilon del medio de arriba del plano
    epsi2 : epsilon del medio de abajo del plano
    hbmu : chemical potential in eV  
    hbgama : collision frequency in eV
    int_v : fraccion de la vel de la c que vale la velocidad uniforme del electron en x, si v = c/int
    b: electron position en z, b < 0 (dipolo en 0 y eje de z apunta hacia abajo) 
    zp : coordenada zp del plano, zp > 0  (dipolo en 0 y eje de z apunta hacia abajo) 
    x : punto x donde miramos el campo en micrones
    omega0 : resonant frequency 
    kappa_factor_omega0 : kappa (collision frequency) = kappa_factor*omega0
    kappa_r_factor : kappa_r = kappa_r_factor*kappa with kappa_r_factor<1    
    Returns
    -------
    alfa effectivo en QE approx
    """

    E = omegac*aux  
    n1 = epsi1*mu1
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1
    k1_3 = k1**3


    alfa = alpha_function(epsi1,omegac,omega0,kappa_factor_omega0,kappa_r_factor)
    
    z_dip_barra_self = k1*2*zp  #NO ESTOY SEGURA DE ESTA EXPONENCIAL EH 

    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)

    alfa_p = (omegac/k1)*1j*(epsi1 + epsi2)/cond
 
    Rp = 2*epsi1/(epsi1 + epsi2)
    
    expB = np.exp(-alfa_p*z_dip_barra_self) 
 
    
    final_x = np.pi*1j*k1_3*Rp*(alfa_p**3)
    final_y = final_x
    final_z = 2*np.pi*1j*k1_3*Rp*(alfa_p**3)
    
    final = (px*final_x + py*final_y + pz*final_z)*expB
    
    rta_final = (1/alfa - final)**(-1)

    return rta_final

#%%

def dipole_moment_ana(omegac,epsi1,epsi2,hbmu,hbgama,int_v,px,py,pz,b,zp,omega0,kappa_factor_omega0,kappa_r_factor):     
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
#    k1_2 = (k0*cte1)**2
 #   n_v1 = int_v/cte1

    alffa_eff = alpha_function_ana(omegac,epsi1,epsi2,hbmu,hbgama,px,py,pz,zp,omega0,kappa_factor_omega0,kappa_r_factor) 

    charge_electron = alffa_eff*4.806e-10*int_v/(2*np.pi)
#    charge_electron = alffa_eff*int_v/(2*np.pi)
    

    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    Rp = 2*epsi1/(epsi1 + epsi2)
    alfa_p = 1j*(epsi1 + epsi2)/(cond*cte1)
    kp = alfa_p*k1

      
    arg = np.abs(b)*omegac*int_v
    K1 = special.kn(1,arg)
    
    kx = omegac*int_v
    expo = np.exp(-np.sqrt(kx**2 + kp**2)*(np.abs(b) + 2*zp))
    expo = np.exp(-kp*(np.abs(b) + 2*zp))
    
    den = np.sqrt(kx**2 + kp**2) - kp
    
    px = charge_electron*1j*(omegac*int_v*K1 - 2*np.pi*1j*Rp*kp*expo/den)
    
    py = px
    
    pz = charge_electron*(-omegac*int_v*K1 - 2*np.pi*1j*Rp*kp*expo/den)
    
    return px,py,pz


#%%

def dipole_moment_num(omegac,epsi1,epsi2,hbmu,hbgama,int_v,px,py,pz,b,zp,omega0,kappa_factor_omega0,kappa_r_factor):     
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
    # k1 = omegac*cte1
#    k1_2 = (k0*cte1)**2
 #   n_v1 = int_v/cte1
    
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)

    alffa_eff = alpha_function_num(omegac,epsi1,epsi2,hbmu,hbgama,px,py,pz,zp,omega0,kappa_factor_omega0,kappa_r_factor)
    
    charge_electron = alffa_eff*4.806e-10*int_v/(2*np.pi)
#    charge_electron = alffa_eff*int_v/(2*np.pi)
    
    cota_inf = 1
    cota_sup = 1001   
    
    
#    term6 = np.sign(z)*pz*Rp*kp_2*J0*exp_electron

    rp_num = lambda u: epsi2*1j*u - epsi1*1j*u - cte1*cond*u**2
    rp_den = lambda u: epsi2*1j*u + epsi1*1j*u - cte1*cond*u**2
    rp = lambda u: rp_num(u)/rp_den(u)
 
      
    arg = np.abs(b)*omegac*int_v
    K1 = special.kn(1,arg)
    
    kx = omegac*int_v

    expo = lambda u: np.exp(-np.sqrt(kx**2 + u**2)*(np.abs(b) + 2*zp))
    
    int_f_re = lambda u: np.real(rp(u)*expo(u))
    int_f_im = lambda u: np.imag(rp(u)*expo(u))
    
    INT_re,err = integrate.quad(int_f_re, cota_inf, cota_sup) 
    INT_im,err = integrate.quad(int_f_im, cota_inf, cota_sup) 
    
    INT = INT_re + 1j*INT_im
    
    px = charge_electron*1j*(omegac*int_v*K1 - INT)
    
    py = px
    
    pz = charge_electron*(-omegac*int_v*K1 - INT)    
    
    return px,py,pz

#%%