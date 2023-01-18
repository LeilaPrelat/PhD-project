
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila

dipole moment
"""
from scipy import integrate
from scipy import special
import numpy as np
import sys
import os 

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/dipole_moment','')
#print('Importar modulos necesarios para este codigo')

try:
    sys.path.insert(1, path_constants)
    from Ag_sigma import sigma_DL
except ModuleNotFoundError:
    print('Ag_sigma.py no se encuentra en ' + path_basic)
    
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
    k1 = omegac*cte1
    k1_3 = k1**3
    
    kappa = kappa_factor_omega0*omega0
    kappa_r = kappa_r_factor*kappa
    A = 3*kappa_r*0.25/k1_3
    
    den = omega0 - omega - 1j*kappa/2

    rta = A/den
    return rta

#%%

def Efield_tot_QE(omegac,epsi1,epsi2,omega_bulk,gamma_in,d,epsilon_b,int_v,b,zp,x,y,z,omega0,kappa_factor_omega0,kappa_r_factor):
    """
    Parameters
    ----------
    omegac : omega/c = k0 en 1/micrometros    
    epsi1 : epsilon del medio de arriba del plano
    epsi2 : epsilon del medio de abajo del plano
    omega_bulk : frecuencia de bulk Hz
    gamma_in : frecuencia de colision en eV
    d : espesor del plano Ag en micrometros
    epsilon_b : permeabilidad electrica del plano 
    int_v : fraccion de la vel de la c que vale la velocidad uniforme del electron en x, si v = c/int
    b: electron position en z, b < 0 (dipolo en 0 y eje de z apunta hacia abajo) 
    zp : coordenada zp del plano, zp > 0  (dipolo en 0 y eje de z apunta hacia abajo) 
    x : punto x donde miramos el campo en micrones
    y : punto y donde miramos el campo en micrones
    z : punto z donde miramos el campo en micrones 
    omega0 : resonant frequency 
    kappa_factor_omega0 : kappa (collision frequency) = kappa_factor*omega0
    kappa_r_factor : kappa_r = kappa_r_factor*kappa with kappa_r_factor<1    
    
    Returns
    -------
    dipole moment (green tensor direct + reflejado)
    habiendo aplicado QE al green tensor + alfa
    """
#    E = omegac*aux      
    omega = omegac*c #=omega/c    
    n1 = epsi1*mu1
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1
    # k1_2 = k1**2
    k1_2 = k1**2
    x_bar = x*k1
    y_bar = y*k1
# green tensor evaluated in : DIPOLO EN EL 0,0,0 
#    x,y,z = 0, 0, 0
#    yD,zD = 0, b

    # green tensor evaluated in x,y,z = 0, 0, 0 and yD,zD = 0, b
    z_dip_barra_dir = k1*np.abs(z-b)               # green tensor evaluated in x,y,z = 0, 0, 0 and yD,zD = 0, b
    atan_dir = lambda w : np.cos(2*np.arctan2(np.abs(y_bar),np.abs(x_bar - w)))  
    
    expB_dir = lambda u: np.exp(-u*z_dip_barra_dir) 
    
    sin_electron = lambda w: np.sin(w*int_v/n1)
    cos_electron = lambda w: np.cos(w*int_v/n1)
    
    
    J0_dir = lambda w,u: special.jn(0,u*np.sqrt((x_bar-w)**2 + y_bar**2)) 
    J2_dir = lambda w,u: special.jn(2,u*np.sqrt((x_bar-w)**2 + y_bar**2)) 

    cota_sup1 = 3001
    
    cota_sup1 = 1001
    cota_sup2 = 80*k1
############ I0 5 ########################################################
    Intdir_B_function_re = lambda w,u: expB_dir(u)*cos_electron(w)*(J0_dir(w,u) + atan_dir(w)*J2_dir(w,u))*(u**2 - 1)
    Intdir_B_function_im = lambda w,u: expB_dir(u)*sin_electron(w)*(J0_dir(w,u) + atan_dir(w)*J2_dir(w,u))*(u**2 - 1)

    intdirB_re,err = integrate.dblquad(Intdir_B_function_re, 1, cota_sup1, -cota_sup2, cota_sup2)
    intdirB_im,err = integrate.dblquad(Intdir_B_function_im, 1, cota_sup1, -cota_sup2, cota_sup2)
#########################################################################

    charge_e = 1.602*1e-19*c
    cte_aux = 1j*(charge_e*omega)/(2*np.pi*int_v)

    cte = k1_2*0.5*cte_aux #signo menos
    rtadir = (intdirB_re + 1j*intdirB_im)*cte 
  
    ####### reflejado ################################

    z_dip_barra_ref = k1*(np.abs(z) + np.abs(b) + 2*zp)  #NO ESTOY SEGURA DE ESTA EXPONENCIAL EH 
    atan_ref = lambda w : np.cos(2*np.arctan2(np.abs(y_bar),np.abs(x_bar + w)))      
    
    cond = 4*np.pi*sigma_DL(omega,omega_bulk,gamma_in,d,epsilon_b)/c

    rp_num = lambda u: epsi2*1j*u - epsi1*1j*u - cte1*cond*(u**2)
    rp_den = lambda u: epsi2*1j*u + epsi1*1j*u - cte1*cond*(u**2)
    rp = lambda u: rp_num(u)/rp_den(u)
    
    rs_num = lambda u: 1j*u - 1j*u - cond/cte1
    rs_den = lambda u: 1j*u + 1j*u + cond/cte1
    rs = lambda u: rs_num(u)/rs_den(u)

    expB_ref = lambda u: np.exp(-u*z_dip_barra_ref) 
    exp_electron = lambda w: np.exp(1j*w*int_v/n1)
    J0_ref = lambda w,u: special.jn(0,u*np.sqrt((x_bar+w)**2 + y_bar**2)) 
    J2_ref = lambda w,u: special.jn(2,u*np.sqrt((x_bar+w)**2 + y_bar**2)) 

    cota_sup1 = 1000
#    cota_sup2 = 80*k1
############ I0 5 + I2 5 #################################################################################################
#    Int5_B_function_re = lambda w,u: np.real((J0(w,u) + J2(w,u))*rs(u)*expB(u)*exp_electron(w))
#    Int5_B_function_im = lambda w,u: np.imag((J0(w,u) + J2(w,u))*rs(u)*expB(u)*exp_electron(w))

#    int5B_re,err = integrate.dblquad(Int5_B_function_re, 0, cota_sup1, -cota_sup2, cota_sup2)
#    int5B_im,err = integrate.dblquad(Int5_B_function_im, 0, cota_sup1, -cota_sup2, cota_sup2)
############ I0 5 + I2 5 + I0 6 + I2 6 #################################################################################################
    IntrefB_function_re = lambda w,u: np.real((J0_ref(w,u) + atan_ref(w)*J2_ref(w,u))*((u**2)*rp(u) + rs(u))*expB_ref(u)*exp_electron(w))
    IntrefB_function_im = lambda w,u: np.imag((J0_ref(w,u) + atan_ref(w)*J2_ref(w,u))*((u**2)*rp(u) + rs(u))*expB_ref(u)*exp_electron(w))

    intrefB_re,err = integrate.dblquad(IntrefB_function_re, 1, cota_sup1, -cota_sup2, cota_sup2)
    intrefB_im,err = integrate.dblquad(IntrefB_function_im, 1, cota_sup1, -cota_sup2, cota_sup2)
    
#    rta5 = (int5B_re + 1j*int5B_im)*cte 
    rtaref = (intrefB_re + 1j*intrefB_im)*cte 
    
    ########### self image dipole ################### caso particular del reflejado pero con J0 = 1 y J2 = 0 ####################
    
    z_dip_barra_self = k1*(2*zp)  #NO ESTOY SEGURA DE ESTA EXPONENCIAL EH 
    expB_self = lambda u: np.exp(-u*z_dip_barra_self) 
    
    IntselfB_function_re = lambda w,u: np.real(((u**2)*rp(u) + rs(u))*expB_self(u)*exp_electron(w))
    IntselfB_function_im = lambda w,u: np.imag(((u**2)*rp(u) + rs(u))*expB_self(u)*exp_electron(w))

    intselfB_re,err = integrate.dblquad(IntselfB_function_re, 1, cota_sup1, -cota_sup2, cota_sup2)
    intselfB_im,err = integrate.dblquad(IntselfB_function_im, 1, cota_sup1, -cota_sup2, cota_sup2)
    
    rtaself = (intselfB_re + 1j*intselfB_im)*cte
    
    alffa = alpha_function(epsi1,omegac,omega0,kappa_factor_omega0,kappa_r_factor)
    
    den  = 1/alffa + rtaself # el menos se convierte en un +
    
    return (rtadir + rtaref)/den


#%%

def Efield_tot(omegac,epsi1,epsi2,omega_bulk,gamma_in,d,epsilon_b,int_v,b,zp,x,y,z,omega0,kappa_factor_omega0,kappa_r_factor):
    """
    Parameters
    ----------
    omegac : omega/c = k0 en 1/micrometros    
    epsi1 : epsilon del medio de arriba del plano
    epsi2 : epsilon del medio de abajo del plano
    omega_bulk : frecuencia de bulk Hz
    gamma_in : frecuencia de colision en eV
    d : espesor del plano Ag en micrometros
    epsilon_b : permeabilidad electrica del plano 
    int_v : fraccion de la vel de la c que vale la velocidad uniforme del electron en x, si v = c/int
    b: electron position en z, b < 0 (dipolo en 0 y eje de z apunta hacia abajo) 
    zp : coordenada zp del plano, zp > 0  (dipolo en 0 y eje de z apunta hacia abajo) 
    x : punto x donde miramos el campo en micrones
    y : punto y donde miramos el campo en micrones
    z : punto z donde miramos el campo en micrones 
    omega0 : resonant frequency 
    kappa_factor_omega0 : kappa (collision frequency) = kappa_factor*omega0
    kappa_r_factor : kappa_r = kappa_r_factor*kappa with kappa_r_factor<1    

    Returns
    -------
    dipole moment (green tensor direct + reflejado)
    """
#    E = omegac*aux  
    omega = omegac*c #=omega/c    
    n1 = epsi1*mu1
    n2 = epsi2*mu2
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1
    k1_2 = k1**2
    x_bar = x*k1
    y_bar = y*k1

# green tensor evaluated in : DIPOLO EN EL 0,0,0 
#    x,y,z = 0, 0, 0
#    yD,zD = 0, b

     # green tensor evaluated in x,y,z = 0, 0, 0 and yD,zD = 0, b
    z_dip_barra_dir = k1*np.abs(z-b)               # green tensor evaluated in x,y,z = 0, 0, 0 and yD,zD = 0, b
    atan_dir = lambda w : np.cos(2*np.arctan2(np.abs(y_bar),np.abs(x_bar - w)))  

    alpha_z1A = lambda u: np.sqrt(1 - u**2)
    alpha_z1B = lambda u: np.sqrt(u**2 - 1)
    
    expA = lambda u: np.exp(1j*alpha_z1A(u)*z_dip_barra_dir) 
    expB = lambda u: np.exp(-alpha_z1B(u)*z_dip_barra_dir)
     
    exp_electron = lambda w: np.exp(1j*w*int_v/n1)
    J0_dir = lambda w,u: special.jn(0,u*np.sqrt((x_bar-w)**2 + y_bar**2)) 
    J2_dir = lambda w,u: special.jn(2,u*np.sqrt((x_bar-w)**2 + y_bar**2)) 
    
    cota_sup1 = 0.95
    cota_inf1 = 1.05
#    cota_sup2 = 1e5
    cota_sup1A = 3000
    
    cota_sup1A = 1000
    cota_sup2A = 80*k1

############ I0 5 ########################################################
    Intdir_A_function_re = lambda w,u: np.real((J0_dir(w,u) + atan_dir(w)*J2_dir(w,u))*u*expA(u)*exp_electron(w)*(alpha_z1A(u) + 1/alpha_z1A(u)))
    Intdir_A_function_im = lambda w,u: np.imag((J0_dir(w,u) + atan_dir(w)*J2_dir(w,u))*u*expA(u)*exp_electron(w)*(alpha_z1A(u) + 1/alpha_z1A(u)))

    intdirA_re,err = integrate.dblquad(Intdir_A_function_re, 0, cota_sup1, -cota_sup2A, cota_sup2A)
    intdirA_im,err = integrate.dblquad(Intdir_A_function_im, 0, cota_sup1, -cota_sup2A, cota_sup2A)

    Intdir_B_function_re = lambda w,u: np.real(1j*(J0_dir(w,u) + atan_dir(w)*J2_dir(w,u))*u*expB(u)*exp_electron(w)*(alpha_z1B(u) - 1/alpha_z1B(u)))
    Intdir_B_function_im = lambda w,u: np.imag(1j*(J0_dir(w,u) + atan_dir(w)*J2_dir(w,u))*u*expB(u)*exp_electron(w)*(alpha_z1B(u) - 1/alpha_z1B(u)))

    intdirB_re,err = integrate.dblquad(Intdir_B_function_re, cota_inf1, cota_sup1A, -cota_sup2A, cota_sup2A)
    intdirB_im,err = integrate.dblquad(Intdir_B_function_im, cota_inf1, cota_sup1A, -cota_sup2A, cota_sup2A)
    
##########################################################################

    charge_e = 1.602*1e-19*c
    cte_aux = (charge_e*omega)/(2*np.pi*int_v) ## se cancela el i con el de la cte y se cancela el menos

    cte = k1_2*0.5*cte_aux
    rtadir = (intdirB_re + intdirA_re + 1j*(intdirB_im + intdirA_im))*cte 

    ####### reflejado ################################

    z_dip_barra_ref =  k1*(np.abs(z) + np.abs(b) + 2*zp)  #NO ESTOY SEGURA DE ESTA EXPONENCIAL EH 
    atan_ref = lambda w : np.cos(2*np.arctan2(np.abs(y_bar),np.abs(x_bar + w)))  
    aux2 = n2/n1
    alpha_z1 = lambda u: np.sqrt(1-u**2) if u<1 else 1j*np.sqrt(u**2-1)
    alpha_z2 = lambda u: np.sqrt(aux2-u**2) if u**2<aux2 else 1j*np.sqrt(u**2-aux2)
        
    cond = 4*np.pi*sigma_DL(omega,omega_bulk,gamma_in,d,epsilon_b)/c

    rp_num = lambda u: epsi2*alpha_z1(u) - epsi1*alpha_z2(u) + cte1*cond*alpha_z1(u)*alpha_z2(u)
    rp_den = lambda u: epsi2*alpha_z1(u) + epsi1*alpha_z2(u) + cte1*cond*alpha_z1(u)*alpha_z2(u)
    rp = lambda u: rp_num(u)/rp_den(u)
    
    rs_num = lambda u: alpha_z1(u) - alpha_z2(u) - cond/cte1
    rs_den = lambda u: alpha_z1(u) + alpha_z2(u) + cond/cte1
    rs = lambda u: rs_num(u)/rs_den(u)

    exp = lambda u: np.exp(1j*alpha_z1(u)*z_dip_barra_ref) 
    J0_ref = lambda w,u: special.jn(0,u*np.sqrt((x_bar+w)**2 + y_bar**2)) 
    J2_ref = lambda w,u: special.jn(2,u*np.sqrt((x_bar+w)**2 + y_bar**2)) 
    
#    cota_sup1 = 0.95
#    cota_inf1 = 1.05
    
    cota_sup1A = 1000
    cota_sup2A = 80*k1

############ I0 5 + I2 5  + I0 6 + I2 6 #################################################################################################
    Intref_A_function_re = lambda w,u: np.real((J0_ref(w,u) + atan_ref(w)*J2_ref(w,u))*u*exp(u)*exp_electron(w)*(-rp(u)*alpha_z1(u) + rs(u)/alpha_z1(u)))
    Intref_A_function_im = lambda w,u: np.imag((J0_ref(w,u) + atan_ref(w)*J2_ref(w,u))*u*exp(u)*exp_electron(w)*(-rp(u)*alpha_z1(u) + rs(u)/alpha_z1(u)))

    intrefA_re,err = integrate.dblquad(Intref_A_function_re, 0, cota_sup1, -cota_sup2A, cota_sup2A)
    intrefA_im,err = integrate.dblquad(Intref_A_function_im, 0, cota_sup1, -cota_sup2A, cota_sup2A)

    Intref_B_function_re = lambda w,u: np.real((J0_ref(w,u) + atan_ref(w)*J2_ref(w,u))*u*exp(u)*exp_electron(w)*(-rp(u)*alpha_z1(u) + rs(u)/alpha_z1(u)))
    Intref_B_function_im = lambda w,u: np.imag((J0_ref(w,u) + atan_ref(w)*J2_ref(w,u))*u*exp(u)*exp_electron(w)*(-rp(u)*alpha_z1(u) + rs(u)/alpha_z1(u)))

    intrefB_re,err = integrate.dblquad(Intref_B_function_re, cota_inf1, cota_sup1A, -cota_sup2A, cota_sup2A)
    intrefB_im,err = integrate.dblquad(Intref_B_function_im, cota_inf1, cota_sup1A, -cota_sup2A, cota_sup2A)
        
    # cte = k1_2*0.5*cte_aux #(*)
    rtaref = (intrefB_re + intrefA_re + 1j*(intrefA_im + intrefB_im))*cte


    ########### self image dipole ################### caso particular del reflejado pero con J0 = 1 y J2 = 0 ####################
    

    z_dip_barra_self = k1*(2*zp)  #NO ESTOY SEGURA DE ESTA EXPONENCIAL EH 
    exp_self = lambda u: np.exp(1j*alpha_z1(u)*z_dip_barra_self) 
    
    Intself_A_function_re = lambda w,u: np.real(u*exp_self(u)*exp_electron(w)*(-rp(u)*alpha_z1(u) + rs(u)/alpha_z1(u)))
    Intself_A_function_im = lambda w,u: np.imag(u*exp_self(u)*exp_electron(w)*(-rp(u)*alpha_z1(u) + rs(u)/alpha_z1(u)))

    intselfA_re,err = integrate.dblquad(Intself_A_function_re, 0, cota_sup1, -cota_sup2A, cota_sup2A)
    intselfA_im,err = integrate.dblquad(Intself_A_function_im, 0, cota_sup1, -cota_sup2A, cota_sup2A)

    Intself_B_function_re = lambda w,u: np.real(u*exp_self(u)*exp_electron(w)*(-rp(u)*alpha_z1(u) + rs(u)/alpha_z1(u)))
    Intself_B_function_im = lambda w,u: np.imag(u*exp_self(u)*exp_electron(w)*(-rp(u)*alpha_z1(u) + rs(u)/alpha_z1(u)))

    intselfB_re,err = integrate.dblquad(Intself_B_function_re, cota_inf1, cota_sup1A, -cota_sup2A, cota_sup2A)
    intselfB_im,err = integrate.dblquad(Intself_B_function_im, cota_inf1, cota_sup1A, -cota_sup2A, cota_sup2A)
    
    rtaself = (intselfB_re + intselfA_re + 1j*(intselfA_im + intselfB_im))*cte


    alffa = alpha_function(epsi1,omegac,omega0,kappa_factor_omega0,kappa_r_factor)
    
    den  = 1/alffa - rtaself 
    
    return (rtadir + rtaref)/den
    

#%%

