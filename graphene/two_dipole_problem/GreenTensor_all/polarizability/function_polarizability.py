
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila

definir polarizabilidad
"""
import numpy as np
import sys
import os 
from scipy.special import spherical_jn,spherical_yn
from scipy import integrate
from scipy import special

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/polarizability','')
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

def alpha_function_eff_QE(omegac,epsi1,epsi2,hbmu,hbgama,int_v,b,zp,xD,omega0,kappa_factor_omega0,kappa_r_factor):
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
    alfa effectivo/alfa en QE approx
    """

    E = omegac*aux  
    n1 = epsi1*mu1
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1
    k1_2 = k1**2

    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)

    rp_num = lambda u: epsi2*1j*u - epsi1*1j*u - cte1*cond*(u**2)
    rp_den = lambda u: epsi2*1j*u + epsi1*1j*u - cte1*cond*(u**2)
    rp = lambda u: rp_num(u)/rp_den(u)
    
    rs_num = lambda u: 1j*u - 1j*u - cond/cte1
    rs_den = lambda u: 1j*u + 1j*u + cond/cte1
    rs = lambda u: rs_num(u)/rs_den(u)

    cota_sup1 = 1000
    cota_sup2 = 80*k1
####       
    cte = 1j*k1_2*0.5 #signo menos

    alfa = alpha_function(epsi1,omegac,omega0,kappa_factor_omega0,kappa_r_factor)
    
    z_dip_barra_self = k1*(2*zp)  #NO ESTOY SEGURA DE ESTA EXPONENCIAL EH 
    expB_self = lambda u: np.exp(-u*z_dip_barra_self) 
    exp_electron = lambda w: np.exp(1j*w*int_v/n1)
    
    IntselfB_function_re = lambda w,u: np.real(((u**2)*rp(u) + rs(u))*expB_self(u)*exp_electron(w))
    IntselfB_function_im = lambda w,u: np.imag(((u**2)*rp(u) + rs(u))*expB_self(u)*exp_electron(w))

    intselfB_re,err = integrate.dblquad(IntselfB_function_re, 1, cota_sup1, -cota_sup2, cota_sup2)
    intselfB_im,err = integrate.dblquad(IntselfB_function_im, 1, cota_sup1, -cota_sup2, cota_sup2)
    
    rtaself = -(intselfB_re + 1j*intselfB_im)*cte
    
    rta_final = (1/alfa - rtaself)**(-1)

    return rta_final/alfa

#%%

def alpha_function_eff(omegac,epsi1,epsi2,hbmu,hbgama,int_v,b,zp,xD,omega0,kappa_factor_omega0,kappa_r_factor): 
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
    alfa effectivo/alfa 
    """

    E = omegac*aux      
    n1 = epsi1*mu1
    n2 = epsi2*mu2
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1
    k1_2 = k1**2
    
    exp_electron = lambda w: np.exp(1j*w*int_v/n1)
    aux2 = n2/n1
    alpha_z1 = lambda u: np.sqrt(1-u**2) if u<1 else 1j*np.sqrt(u**2-1)
    alpha_z2 = lambda u: np.sqrt(aux2-u**2) if u<aux2 else 1j*np.sqrt(u**2-aux2)
        
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)

    rp_num = lambda u: epsi2*alpha_z1(u) - epsi1*alpha_z2(u) + cte1*cond*alpha_z1(u)*alpha_z2(u)
    rp_den = lambda u: epsi2*alpha_z1(u) + epsi1*alpha_z2(u) + cte1*cond*alpha_z1(u)*alpha_z2(u)
    rp = lambda u: rp_num(u)/rp_den(u)
    
    rs_num = lambda u: alpha_z1(u) - alpha_z2(u) - cond/cte1
    rs_den = lambda u: alpha_z1(u) + alpha_z2(u) + cond/cte1
    rs = lambda u: rs_num(u)/rs_den(u)

    cota_sup1 = 0.95
    cota_inf1 = 1.05
#    cota_sup2 = 1e5
    cota_sup1A = 1000
    cota_sup2A = 80*k1
         
    alfa = alpha_function(epsi1,omegac,omega0,kappa_factor_omega0,kappa_r_factor)
    cte = k1_2*0.5
    
    z_dip_barra_self = k1*2*zp  #NO ESTOY SEGURA DE ESTA EXPONENCIAL EH 
    exp_self = lambda u: np.exp(1j*alpha_z1(u)*z_dip_barra_self) 
    
    arg = lambda u : u*xD*k1*np.sqrt(2)
    J0_self = lambda u: special.jn(0,arg(u)) 
    J2_self = lambda u: special.jn(2,arg(u)) 
    
    Intself_A_function_re = lambda w,u: np.real(u*exp_self(u)*exp_electron(w)*(-rp(u)*alpha_z1(u) + rs(u)/alpha_z1(u))*(J0_self(u) + J2_self(u)))
    Intself_A_function_im = lambda w,u: np.imag(u*exp_self(u)*exp_electron(w)*(-rp(u)*alpha_z1(u) + rs(u)/alpha_z1(u))*(J0_self(u) + J2_self(u)))

    intselfA_re,err = integrate.dblquad(Intself_A_function_re, 0, cota_sup1, -cota_sup2A, cota_sup2A)
    intselfA_im,err = integrate.dblquad(Intself_A_function_im, 0, cota_sup1, -cota_sup2A, cota_sup2A)

    Intself_B_function_re = lambda w,u: np.real(u*exp_self(u)*exp_electron(w)*(-rp(u)*alpha_z1(u) + rs(u)/alpha_z1(u))*(J0_self(u) + J2_self(u)))
    Intself_B_function_im = lambda w,u: np.imag(u*exp_self(u)*exp_electron(w)*(-rp(u)*alpha_z1(u) + rs(u)/alpha_z1(u))*(J0_self(u) + J2_self(u)))

    intselfB_re,err = integrate.dblquad(Intself_B_function_re, cota_inf1, cota_sup1A, -cota_sup2A, cota_sup2A)
    intselfB_im,err = integrate.dblquad(Intself_B_function_im, cota_inf1, cota_sup1A, -cota_sup2A, cota_sup2A)
    
    rtaself = (intselfB_re + intselfA_re + 1j*(intselfA_im + intselfB_im))*cte

    rta_final = (1/alfa - rtaself)**(-1)

    return rta_final/alfa


#%%

def alpha_functionMIE(omegac,epsi1,R):
    """
    Parameters
    ----------
    omegac : omega/c = k0 en 1/micrometros    
    epsilon1 : permeabilidad electrica del medio 1
    R : radius of the sphere in micrometers
    Returns
    -------
    using the Mie coefficient of the order = 1 for the polarizability
    """
    # omega = omegac*c
    k0 = omegac
    n1 = epsi1*mu1
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1
    k1_3 = k1**3
    
    arg1 = k0*R*cte1
    arg0 = k0*R

    j1 = spherical_jn(1,arg1)
    der_j1 = spherical_jn(1,arg1,1)
    
    j0 = spherical_jn(1,arg0)
    der_j0 = spherical_jn(1,arg0,1)
    
    
    y0 = spherical_yn(1,arg0)
    der_y0 = spherical_yn(1,arg0,1)
    
    h0_mas = 1j*j0 - y0
    der_h0_mas = 1j*der_j0 - der_y0

    term1 = j1 + arg1*der_j1
    term0 = j0 + arg0*der_j0
    term2 = h0_mas + arg0*der_h0_mas
    
    num = -j0*term1 + n1*term0*j1
    
    den = h0_mas*term1 - n1*term2*j1 

    return 3*0.5*k1_3*num/den

#%%
