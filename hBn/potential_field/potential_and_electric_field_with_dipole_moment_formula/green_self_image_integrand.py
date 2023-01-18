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
path_constants =  path_basic.replace('/potential_field/potential_and_electric_field_with_dipole_moment_formula','')
#print('Importar modulos necesarios para este codigo')

try:
    sys.path.insert(1, path_constants)
    from hBn_PP import hBn_lambda_p, hBn_Rp, epsilon_x, epsilon_z
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


def green_self_num(omegac,epsi1,epsi3,d_nano,zp,u):
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

    kz1 = np.sqrt(epsi1 - u**2) if (u**2 <= epsi1) else 1j*np.sqrt(u**2 - epsi1)
    kz2 = np.sqrt(epsi_HBN_par - (epsi_HBN_par/epsi_HBN_perp)*u**2) 
    kz3 =  np.sqrt(epsi3 - u**2) if (u**2 <= epsi3) else 1j*np.sqrt(u**2 - epsi3)

    r12 =  (kz1*epsi_x - kz2*epsi1)/(kz1*epsi_x + kz2*epsi1)
    r21 = (kz2*epsi1 - kz1*epsi_x)/(kz2*epsi1 + kz1*epsi_x)
    r23 = (kz2*epsi3 - kz3*epsi_x)/(kz2*epsi3 + kz3*epsi_x)
    

    exp_fresnel = np.exp(1j*2*kz2*omegac*d_micros)
    
    cte_t = np.sqrt(epsi1*epsi_x)
    t12 = 2*kz1*cte_t/(kz1*epsi_x + kz2*epsi1)
    t21 = 2*kz2*cte_t/(kz2*epsi1 + kz1*epsi_x)

    rp_num =  t12*t21*r23*exp_fresnel
    rp_den =  1 - r21*r23*exp_fresnel
    rp =  r12 +  rp_num/rp_den   
    
#    print(np.real(rp(10)))
    
#    rs_num = lambda u: 1j*u - 1j*u - cond/cte1
#    rs_den = lambda u: 1j*u + 1j*u + cond/cte1
#    rs = lambda u: rs_num(u)/rs_den(u)

    
####       
    cte_x = k1_3*0.5 #signo menos

    
    z_dip_barra_self = k1*2*zp  #NO ESTOY SEGURA DE ESTA EXPONENCIAL EH 
    expB_self =  np.exp(-u*z_dip_barra_self) 

    IntselfB_function_re_xx =  u**2*np.real(rp)*expB_self
    IntselfB_function_im_xx =  u**2*np.imag(rp)*expB_self

    
    rtaself_x = (IntselfB_function_re_xx + 1j*IntselfB_function_im_xx)*cte_x
    rtaself_y = rtaself_x
    rtaself_z = 2*rtaself_x
    

    return rtaself_x, rtaself_y, rtaself_z

#%%
    

def green_self_pole_aprox(omegac,epsi1,epsi3,d_nano,zp,u):
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

    d_micros = d_nano*1e-3
    Rp = hBn_Rp(E,epsi1,epsi3)
    lambda_p_v = hBn_lambda_p(E,epsi1,epsi3)*d_micros
    kp = 2*np.pi/lambda_p_v
    alfa_p = kp/omegac
    
    rp = Rp*alfa_p/(u - alfa_p)     
#    print(np.real(rp(10)))
    
#    rs_num = lambda u: 1j*u - 1j*u - cond/cte1
#    rs_den = lambda u: 1j*u + 1j*u + cond/cte1
#    rs = lambda u: rs_num(u)/rs_den(u)

    
####       
    cte_x = k1_3*0.5 #signo menos

    
    z_dip_barra_self = k1*2*zp  #NO ESTOY SEGURA DE ESTA EXPONENCIAL EH 
    expB_self =  np.exp(-u*z_dip_barra_self) 

    IntselfB_function_re_xx =  u**2*np.real(rp)*expB_self
    IntselfB_function_im_xx = u**2*np.imag(rp)*expB_self
    
    
    rtaself_x = (IntselfB_function_re_xx + 1j*IntselfB_function_im_xx)*cte_x
    rtaself_y = rtaself_x
    rtaself_z = 2*rtaself_x
    

    return rtaself_x, rtaself_y, rtaself_z


#%%
