#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila
"""
import numpy as np
import sys
import os 
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


def green_tensor_self_num(omegac,epsi1,epsi2,hbmu,hbgama,px,py,pz,zp):
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
    k1_3 = k1**3

    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)

    rp_num = lambda u: epsi2*1j*u - epsi1*1j*u - cte1*cond*(u**2)
    rp_den = lambda u: epsi2*1j*u + epsi1*1j*u - cte1*cond*(u**2)
    rp = lambda u: rp_num(u)/rp_den(u)
    
#    rs_num = lambda u: 1j*u - 1j*u - cond/cte1
#    rs_den = lambda u: 1j*u + 1j*u + cond/cte1
#    rs = lambda u: rs_num(u)/rs_den(u)

    cota_sup1 = 1000*k1
####       
    cte_x = k1_3*0.5 #signo menos
    
    z_dip_barra_self = k1*2*zp  #NO ESTOY SEGURA DE ESTA EXPONENCIAL EH 
    expB_self = lambda u: np.exp(-u*z_dip_barra_self) 

    IntselfB_function_re_x = lambda u: np.real(((u**2)*rp(u))*expB_self(u))
    IntselfB_function_im_x = lambda u: np.imag(((u**2)*rp(u))*expB_self(u))

    intselfB_re_x,err = integrate.quad(IntselfB_function_re_x, 1, cota_sup1)
    intselfB_im_x,err = integrate.quad(IntselfB_function_im_x, 1, cota_sup1)
    
    
    rtaself_x = (py + px)*(intselfB_re_x + 1j*intselfB_im_x)*cte_x
    
    
    
    

    IntselfB_function_re_z = lambda u: np.real((u**2)*rp(u)*expB_self(u))
    IntselfB_function_im_z = lambda u: np.imag((u**2)*rp(u)*expB_self(u))

    intselfB_re_z,err = integrate.quad(IntselfB_function_re_z, 1, cota_sup1)
    intselfB_im_z,err = integrate.quad(IntselfB_function_im_z, 1, cota_sup1)    
    
    rtaself_z = pz*(intselfB_re_z + 1j*intselfB_im_z)*cte_x*2
    
    
    rtaself = rtaself_x + rtaself_z
    
    rta_final =  - rtaself

    return rta_final

#%%

def green_tensor_self_ana(omegac,epsi1,epsi2,hbmu,hbgama,px,py,pz,zp):
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
    
    z_dip_barra_self = k1*2*zp  #NO ESTOY SEGURA DE ESTA EXPONENCIAL EH 

    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)

    alfa_p = 1j*(epsi1 + epsi2)/cond
 
    Rp = 2*epsi1/(epsi1 + epsi2)
    
    expB = np.exp(-alfa_p*z_dip_barra_self) 
 
    
    final_x = np.pi*1j*k1_3*Rp*(alfa_p**3)
    final_y = final_x
    final_z = 2*np.pi*1j*k1_3*Rp*(alfa_p**3)
    
    final = (px*final_x + py*final_y + pz*final_z)*expB
    
    rta_final = - final

    return rta_final

#%%
