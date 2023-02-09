
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila

campo externo directo + reflejado numerico con la convencion de z hacia abajo
en z = 0
No hay solucion analitica porque es para cualquier x,y (solo hay sol analitica en x=y=0)
Dipolo en el 0,0
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

"""
def alpha_function(omegac,epsi1,epsi2,hbmu,hbargama,z_plane):
    
    E = omegac*aux
    n1 = epsi1*mu1
    n2 = epsi2*mu2
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1
    k1_3 = k1**3
    
    aux2 = n2/n1
    
    def alpha_z1(u): 
        return np.sqrt(1-u**2) if u<1 else 1j*np.sqrt(u**2-1)
    
    def alpha_z2(u):
        return np.sqrt(aux2-u**2) if u<aux2 else 1j*np.sqrt(u**2-aux2)
    
    sigmatot = sigma_DL(E,hbmu,hbargama)
    cond3 = sigmatot*alfac*4*pi  #sigma devuelve la conductividad teniendo que multiplicar por alfac*c ---> no hay que dividir por c

    
    exp = lambda u: np.exp(1j*k1*(alpha_z1(u)-alpha_z2(u))*z_plane) 

    num = lambda u: 2*mu2*alpha_z1(u)*2
    den = lambda u: mu2*alpha_z1(u) + mu1*alpha_z2(u) - cond3*mu2*mu1/cte1 
    func = lambda u: num(u)/den(u) 

    F_re = lambda u: np.real(func(u)*exp(u))
    F_im = lambda u: np.imag(func(u)*exp(u))

    cota_supB = 600
    cota_supA = 0.95
    cota_infB = 1.05
    
#    cota_sup2A = 80*k1

    intA_re,err = integrate.quad(F_re, 0, cota_supA)
    intB_re,err = integrate.quad(F_re, cota_infB, cota_supB)
    
    intA_im,err = integrate.quad(F_im, 0, cota_supA)
    intB_im,err = integrate.quad(F_im, cota_infB, cota_supB)

    cte = 3*k1_3/2

    final = intA_re + intB_re + 1j*(intA_im + intB_im)
    
    return final*cte

#%%

def alpha_functionQE(omegac,epsi1,epsi2,hbmu,hbargama,z_plane):
    
    E = omegac*aux
    n1 = epsi1*mu1
    # n2 = epsi2*mu2
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1
    k1_3 = k1**3
    
    def alpha_z1(u): 
        return 1j*u
    
    def alpha_z2(u):
        return 1j*u
    
    sigmatot = sigma_DL(E,hbmu,hbargama)
    cond3 = sigmatot*alfac*4*pi  #sigma devuelve la conductividad teniendo que multiplicar por alfac*c ---> no hay que dividir por c

    
    exp = lambda u: np.exp(1j*k1*(alpha_z1(u)-alpha_z2(u))*z_plane) 

    num = lambda u: 2*mu2*alpha_z1(u)*2
    den = lambda u: mu2*alpha_z1(u) + mu1*alpha_z2(u) - cond3*mu2*mu1/cte1 
    func = lambda u: num(u)/den(u) 

    F_re = lambda u: np.real(func(u)*exp(u))
    F_im = lambda u: np.imag(func(u)*exp(u))

    cota_supB = 600
    # cota_supA = 0.95
    # cota_infB = 1.05
    
#    cota_sup2A = 80*k1

    int_re,err = integrate.quad(F_re, 1, cota_supB)
    int_im,err = integrate.quad(F_im, 1, cota_supB)

    cte = 3*k1_3/2

    final = int_re + int_re + 1j*(int_im + int_im)
    
    return final*cte

"""
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

def Efield_tot_QE(omegac,epsi1,epsi2,hbmu,hbgama,int_v,b,zp,xD,omega0,kappa_factor_omega0,kappa_r_factor):
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
    dipole moment (green tensor direct + reflejado)
    habiendo aplicado QE al green tensor + alfa
    """
    E = omegac*aux      
    # omega = omegac*c #=omega/c    
    n1 = epsi1*mu1
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1
    # k1_2 = k1**2
    k1_2 = k1**2
    xD_bar = xD*k1
    # y_bar = y*k1
# green tensor evaluated in : DIPOLO EN EL 0,0,0 
#    x,y,z = 0, 0, 0
#    yD,zD = 0, b

    # green tensor evaluated in x,y,z = 0, 0, 0 and yD,zD = 0, b
    z_dip_barra_dir = k1*np.abs(-b)  ### evaluado en -b CREO
    atan_dir = lambda w : np.cos(2*np.arctan2(0,np.abs(xD_bar - w)))  
    
    expB_dir = lambda u: np.exp(-u*z_dip_barra_dir) 
    
    sin_electron = lambda w: np.sin(w*int_v/n1)
    cos_electron = lambda w: np.cos(w*int_v/n1)
    
    
    J0_dir = lambda w,u: special.jn(0,u*(xD_bar-w)) 
    J2_dir = lambda w,u: special.jn(2,u*(xD_bar-w)) 

    cota_sup1 = 3001
    cota_sup2 = 80*k1
############ I0 5 ########################################################
    Intdir_B_function_re = lambda w,u: expB_dir(u)*cos_electron(w)*(J0_dir(w,u) + atan_dir(w)*J2_dir(w,u))*(u**2 - 1)
    Intdir_B_function_im = lambda w,u: expB_dir(u)*sin_electron(w)*(J0_dir(w,u) + atan_dir(w)*J2_dir(w,u))*(u**2 - 1)

    intdirB_re,err = integrate.dblquad(Intdir_B_function_re, 1, cota_sup1, -cota_sup2, cota_sup2)
    intdirB_im,err = integrate.dblquad(Intdir_B_function_im, 1, cota_sup1, -cota_sup2, cota_sup2)
#########################################################################

    pi2 = np.pi**2
    cte_aux = alfac*(c**2)/(pi2*int_v)

    cte = 1j*k1_2*0.5 #signo menos
    rtadir = (intdirB_re + 1j*intdirB_im)*cte 
  
    ####### reflejado ################################

    z_dip_barra_ref = k1*(b + 2*zp)  #NO ESTOY SEGURA DE ESTA EXPONENCIAL EH 
    atan_ref = lambda w : np.cos(2*np.arctan2(0,np.abs(xD_bar + w)))      
    
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)

    rp_num = lambda u: epsi2*1j*u - epsi1*1j*u - cte1*cond*(u**2)
    rp_den = lambda u: epsi2*1j*u + epsi1*1j*u - cte1*cond*(u**2)
    rp = lambda u: rp_num(u)/rp_den(u)
    
    rs_num = lambda u: 1j*u - 1j*u - cond/cte1
    rs_den = lambda u: 1j*u + 1j*u + cond/cte1
    rs = lambda u: rs_num(u)/rs_den(u)

    expB_ref = lambda u: np.exp(-u*z_dip_barra_ref) 
    exp_electron = lambda w: np.exp(1j*w*int_v/n1)
    J0_ref = lambda w,u: special.jn(0,u*(xD_bar+w)) 
    J2_ref = lambda w,u: special.jn(2,u*(xD_bar+w)) 

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
    rtaref = -(intrefB_re + 1j*intrefB_im)*cte 
    
    ########### self image dipole ################### caso particular del reflejado pero con J0 = 1 y J2 = 0 ####################
    
    z_dip_barra_self = k1*(2*zp)  #NO ESTOY SEGURA DE ESTA EXPONENCIAL EH 
    expB_self = lambda u: np.exp(-u*z_dip_barra_self) 
    
    IntselfB_function_re = lambda w,u: np.real(((u**2)*rp(u) + rs(u))*expB_self(u)*exp_electron(w))
    IntselfB_function_im = lambda w,u: np.imag(((u**2)*rp(u) + rs(u))*expB_self(u)*exp_electron(w))

    intselfB_re,err = integrate.dblquad(IntselfB_function_re, 1, cota_sup1, -cota_sup2, cota_sup2)
    intselfB_im,err = integrate.dblquad(IntselfB_function_im, 1, cota_sup1, -cota_sup2, cota_sup2)
    
    rtaself = -(intselfB_re + 1j*intselfB_im)*cte
    
    alffa = alpha_function(epsi1,omegac,omega0,kappa_factor_omega0,kappa_r_factor)
    
    den  = 1/alffa + rtaself # el menos se convierte en un +
    
    return (rtadir + rtaref)**2*cte_aux/den


#%%

def Efield_tot(omegac,epsi1,epsi2,hbmu,hbgama,int_v,b,zp,xD,omega0,kappa_factor_omega0,kappa_r_factor):
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
    dipole moment (green tensor direct + reflejado)
    """
    E = omegac*aux  
    # omega = omegac*c #=omega/c    
    n1 = epsi1*mu1
    n2 = epsi2*mu2
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1
    k1_2 = k1**2
    xD_bar = xD*k1
    # y_bar = y*k1

# green tensor evaluated in : DIPOLO EN EL 0,0,0 
#    x,y,z = 0, 0, 0
#    yD,zD = 0, b

     # green tensor evaluated in x,y,z = 0, 0, 0 and yD,zD = 0, b
    z_dip_barra_dir = k1*np.abs(-b)               # green tensor evaluated in x,y,z = 0, 0, 0 and yD,zD = 0, b
    atan_dir = lambda w : np.cos(2*np.arctan2(0,np.abs(xD_bar - w)))  

    alpha_z1A = lambda u: np.sqrt(1 - u**2)
    alpha_z1B = lambda u: np.sqrt(u**2 - 1)
    
    expA = lambda u: np.exp(1j*alpha_z1A(u)*z_dip_barra_dir) 
    expB = lambda u: np.exp(-alpha_z1B(u)*z_dip_barra_dir)
     
    exp_electron = lambda w: np.exp(1j*w*int_v/n1)
    J0_dir = lambda w,u: special.jn(0,u*(xD_bar-w))
    J2_dir = lambda w,u: special.jn(2,u*(xD_bar-w))
    
    cota_sup1 = 0.95
    cota_inf1 = 1.05
#    cota_sup2 = 1e5
    cota_sup1A = 3000
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
    pi2 = np.pi**2
    cte_aux = alfac*(c**2)/(pi2*int_v)

    cte = k1_2*0.5
    rtadir = (intdirB_re + intdirA_re + 1j*(intdirB_im + intdirA_im))*cte 

    ####### reflejado ################################

    z_dip_barra_ref =  k1*(b + 2*zp)  #NO ESTOY SEGURA DE ESTA EXPONENCIAL EH 
    atan_ref = lambda w : np.cos(2*np.arctan2(0,np.abs(xD_bar + w)))  
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

    exp = lambda u: np.exp(1j*alpha_z1(u)*z_dip_barra_ref) 
    J0_ref = lambda w,u: special.jn(0,u*(xD_bar+w))
    J2_ref = lambda w,u: special.jn(2,u*(xD_bar+w)) 
    
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
    
    return (rtadir + rtaref)**2*cte_aux/den
    

#%%

