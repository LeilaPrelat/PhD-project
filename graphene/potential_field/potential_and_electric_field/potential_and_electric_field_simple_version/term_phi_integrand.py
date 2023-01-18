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
path_constants =  path_basic.replace('potential_field/potential_and_electric_field/potential_and_electric_field_simple_version/potential_1_dipolo_analytical_vs_num','')
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


def intergand_num(omegac,epsi1,epsi2,hbmu,hbgama,x,y,zp,u): # u es alfa parallel y rp es fresnel coefficient
    """
    Parameters
    ----------
    omegac : omega/c = k0 en 1/micrometros    
    epsi1 : epsilon del medio de arriba del plano
    epsi2 : epsilon del medio de abajo del plano
    hbmu : chemical potential in eV  
    hbgama : collision frequency in eV
    int_v : fraccion de la vel de la c que vale la velocidad uniforme del electron en x, si v = c/int
    zp : coordenada zp del plano, zp > 0  (dipolo en 0 y eje de z apunta hacia abajo) 
    x : punto x donde miramos el campo en micrones    
    Returns
    -------
    alfa effectivo en QE approx
    """

    E = omegac*aux 
    k0 = omegac
    n1 = epsi1*mu1
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1
    k1_2 = k1**2

    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)

    rp_num = epsi2*1j*u - epsi1*1j*u - cond*(u**2)
    rp_den = epsi2*1j*u + epsi1*1j*u - cond*(u**2)
    rp = rp_num/rp_den

    R = np.sqrt(x**2 + y**2)
    J1 = special.jv(1,u*k1*R)
    
    z_dip_barra_self = k1*2*zp  
    expB_self = np.exp(-u*z_dip_barra_self) 

    IntselfB_function_re_x = np.real(u*rp*J1*expB_self)
    IntselfB_function_im_x = np.imag(u*rp*J1*expB_self)
 
    rta_final =  k1_2*(IntselfB_function_re_x + 1j*IntselfB_function_im_x)

    return rta_final

#%%


def intergand_num2(omegac,epsi1,epsi2,hbmu,hbgama,x,y,zp,u): # u es alfa parallel y rp es pole aprox
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
    k0 = omegac
    n1 = epsi1*mu1
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1
    k1_2 = k1**2

    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    omega = omegac*c

    alfa_p = 1j*(epsi1 + epsi2)/cond
 
    Rp = 2*epsi1/(epsi1 + epsi2)
#    rp_num = epsi2*1j*u - epsi1*1j*u - cond*(u**2)
#    rp_den = epsi2*1j*u + epsi1*1j*u - cond*(u**2)
    rp = Rp*alfa_p/(u - alfa_p)
    
#    rs_num = lambda u: 1j*u - 1j*u - cond/cte1
#    rs_den = lambda u: 1j*u + 1j*u + cond/cte1
#    rs = lambda u: rs_num(u)/rs_den(u)

    # cota_sup1 = 1000*k1
    R = np.sqrt(x**2 + y**2)
    J1 = special.jv(1,u*k1*R)
    
    z_dip_barra_self = k1*2*zp  
    expB_self = np.exp(-u*z_dip_barra_self) 

    IntselfB_function_re_x = np.real(u*rp*J1*expB_self)
    IntselfB_function_im_x = np.imag(u*rp*J1*expB_self)
 
    rta_final =  k1_2*(IntselfB_function_re_x + 1j*IntselfB_function_im_x)

    return rta_final

#%%