
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
from scipy import integrate
from scipy import special

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/two_dipoles' ,'')
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

def alpha_function_eff_num(omegac,epsi1,epsi2,hbmu,hbgama,int_v,zp,xD,yD,zD,omega0,kappa_factor_omega0,kappa_r_factor):
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
    k1_3 = k1**3

    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)

    rp_num = lambda u: epsi2*1j*u - epsi1*1j*u - cte1*cond*(u**2)
    rp_den = lambda u: epsi2*1j*u + epsi1*1j*u - cte1*cond*(u**2)
    rp = lambda u: rp_num(u)/rp_den(u)

    cota_sup1 = 1000*k1
####       
    cte = k1_3*0.5 #signo menos

    alfa = alpha_function(epsi1,omegac,omega0,kappa_factor_omega0,kappa_r_factor)
    
    z_dip_barra_self = k1*(2*np.abs(zp) + np.abs(zD))
    expB_self = lambda u: np.exp(-u*z_dip_barra_self)
    
    term = np.sqrt((2*np.abs(xD))**2 + (2*np.abs(yD))**2)
    arg = lambda u: u*k1*term
    phi_self = np.arctan2(np.abs(yD), np.abs(xD))
    cos_phi = np.cos(2*phi_self)
        
    J0_self = lambda u: special.jn(0,arg(u)) 
    J2_self = lambda u: special.jn(2,arg(u)) 
    bessel_self = lambda u: J0_self(u) + J2_self(u)*cos_phi
    
    IntselfB_function_re = lambda u: np.real((u**2)*rp(u)*expB_self(u)*bessel_self(u))
    IntselfB_function_im = lambda u: np.imag((u**2)*rp(u)*expB_self(u)*bessel_self(u))

    intselfB_re,err = integrate.quad(IntselfB_function_re, 1, cota_sup1)
    intselfB_im,err = integrate.quad(IntselfB_function_im, 1, cota_sup1)
    
    rtaself = -(intselfB_re + 1j*intselfB_im)*cte
    
    rta_final = (1/alfa - rtaself)**(-1)

    return rta_final/alfa

#%%


#def alpha_function_eff_ana(omegac,epsi1,epsi2,hbmu,hbgama,int_v,zp,xD,yD,zD,omega0,kappa_factor_omega0,kappa_r_factor):
#    """
#    Parameters
#    ----------
#    omegac : omega/c = k0 en 1/micrometros    
#    epsi1 : epsilon del medio de arriba del plano
#    epsi2 : epsilon del medio de abajo del plano
#    hbmu : chemical potential in eV  
#    hbgama : collision frequency in eV
#    int_v : fraccion de la vel de la c que vale la velocidad uniforme del electron en x, si v = c/int
#    b: electron position en z, b < 0 (dipolo en 0 y eje de z apunta hacia abajo) 
#    zp : coordenada zp del plano, zp > 0  (dipolo en 0 y eje de z apunta hacia abajo) 
#    x : punto x donde miramos el campo en micrones
#    omega0 : resonant frequency 
#    kappa_factor_omega0 : kappa (collision frequency) = kappa_factor*omega0
#    kappa_r_factor : kappa_r = kappa_r_factor*kappa with kappa_r_factor<1    
#    Returns
#    -------
#    alfa effectivo/alfa en QE approx
#    """
#
#    E = omegac*aux  
#    n1 = epsi1*mu1
#    cte1 = np.sqrt(n1)
#    k1 = omegac*cte1
#    k1_3 = k1**3
#
#    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
#
#    Rp = 2*epsi1/(epsi1 + epsi2)
#    alfa_p = 1j*(epsi1 + epsi2)/(cond*cte1)
#    rp = Rp*alfa_p
##
#####       
#    cte = k1_3*0.5 #signo menos
#
#    alfa = alpha_function(epsi1,omegac,omega0,kappa_factor_omega0,kappa_r_factor)
#    
#    z_dip_barra_self = k1*(2*np.abs(zp) + 2*np.abs(zD))
#    expB_self = np.exp(-alfa_p*z_dip_barra_self)
#    
#    term = np.sqrt((2*np.abs(xD))**2 + (2*np.abs(yD))**2)
#    arg = alfa_p*k1*term
#    phi_self = np.arctan2(np.abs(yD), np.abs(xD))
#    cos_phi = np.cos(2*phi_self)
#        
#    J0_self = special.jn(0,arg) 
#    J2_self = special.jn(2,arg) 
#    bessel_self = J0_self + J2_self*cos_phi
#    
#    IntselfB_function_re = np.real((alfa_p**2)*rp*expB_self*bessel_self)
#    IntselfB_function_im = np.imag((alfa_p**2)*rp*expB_self*bessel_self)
#
##    intselfB_re,err = integrate.quad(IntselfB_function_re, 1, cota_sup1)
##    intselfB_im,err = integrate.quad(IntselfB_function_im, 1, cota_sup1)
#    
#    rtaself = -(IntselfB_function_re + 1j*IntselfB_function_im)*cte
#    
#    rta_final = (1/alfa - rtaself)**(-1)
#
#    return rta_final/alfa

#%%

def green_dir_ana(omegac,epsi1,xD1,yD1,zD1,xD2,yD2,zD2):
    
    """
    Parameters
    ----------
    omegac : omega/c = k0 en 1/micrometros    
    epsi1 : epsilon del medio de arriba del plano
    x : coordenada x 
    y : coordenada y 
    z : coordenada z
    xD : coordenada x del dipolo 
    yD : coordenada y del dipolo 
    zD : coordenada z del dipolo 
    Returns
    -------
    Gxx direct (self interaction of the dipole)
    con z hacia abajo (convencion del paper)
    analitico luego de aplicar QE
    """
    
    n1 = epsi1*mu1
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1
    k1_3 = k1**3

    Rbarra = k1*np.sqrt((xD1-xD2)**2 + (yD1-yD2)**2)  #adimensional
    Rbarra_2 = Rbarra**2

    z_dip_barra = k1*np.abs(zD1-zD2) #adimensional
    z_dip_barra_2 = z_dip_barra**2
    
    phi = np.arctan2(np.abs(yD1-yD2),np.abs(xD1-xD2))

    aux0 = z_dip_barra**2 + Rbarra**2
    aux1 = aux0**(-1/2)
    aux2 = aux0**(1/2) - z_dip_barra
    
    I0_5 = 0.5*k1_3*aux1
    I2_5 = 0.5*(k1_3/Rbarra_2)*np.cos(2*phi)*aux1*(aux2**2) 
    I0_6 = 0.5*k1_3*(3*z_dip_barra_2/(aux0**(5/2)) - aux0**(-3/2))

    # derivada del I6 orden 2 ### overleaf
    term1 = 3*z_dip_barra_2*(aux2**2)/(aux0**(5/2))
    term2 = 2*z_dip_barra*aux2*(2*z_dip_barra*aux1-2)/(aux0**(3/2))
    term3 = aux2**2/(aux0**(3/2))
    term4 = aux2*(-2*z_dip_barra_2*(aux1**3) + 2*aux1)*aux1
    term5 = 2*((z_dip_barra*aux1-1)**2)*aux1
    
    I2_6 = 0.5*(k1_3/Rbarra_2)*np.cos(2*phi)*(term1 - term2 - term3 + term4 + term5)
    termf_final = I0_5 + I2_5 - I0_6 - I2_6
    
    return termf_final  
    

#%%

#def green_ref_ana(omegac,epsi1,epsi2,hbmu,hbgama,xD1,yD1,zD1,xD2,yD2,zD2,zp):
#    
#    """
#    Parameters
#    ----------
#    omegac : omega/c = k0 en 1/micrometros    
#    epsi1 : epsilon del medio de arriba del plano
#    x : coordenada x 
#    y : coordenada y 
#    z : coordenada z
#    xD : coordenada x del dipolo 
#    yD : coordenada y del dipolo 
#    zD : coordenada z del dipolo 
#    Returns
#    -------
#    Gxx direct (self interaction of the dipole)
#    con z hacia abajo (convencion del paper)
#    analitico luego de aplicar QE
#    """
#    E = omegac*aux
#    k0 = omegac #=omega/c
#    
#    n1 = epsi1*mu1
#    cte1 = np.sqrt(n1)
#    k1 = omegac*cte1
#    k1_3 = k1**3
#
#    Rbarra = k1*np.sqrt((xD1 + xD2)**2 + (yD1 + yD2)**2)  #adimensional
##    Rbarra_2 = Rbarra**2
#
#    z_dip_barra = k1*np.abs(zD1 + zD2) #adimensional
##    z_dip_barra_2 = z_dip_barra**2
#    
#    phi = np.arctan2(np.abs(yD1 + yD2), np.abs(xD1 + xD2))
#
#    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
#    
#    alpha_p = k0*1j*(epsi1 + epsi2)/(cond*k1)
#    Rp = 2*epsi1/(epsi1 + epsi2)
#  
#    J0 = special.jn(0,alpha_p*Rbarra) 
#    J2 = special.jn(2,alpha_p*Rbarra) 
#
#    alpha_p3 = alpha_p**3
#    
#    expo = np.exp(-alpha_p*(z_dip_barra + 2*zp))
#    
#    rta = k1_3*0.5*Rp*alpha_p3*(J0 + J2*np.cos(2*phi))*expo
#    
#    return rta  

#%%
    
def green_tensor_ref_rp_fresnel_QE(omegac,epsi1,epsi2,hbmu,hbgama,x,y,z,xe,ye,ze,zp):
    """
    Parameters
    ----------
    omegac : omega/c = k0 en 1/micrometros    
    epsi1 : epsilon del medio de arriba del plano
    epsi2 : epsilon del medio de abajo del plano
    hbmu : chemical potential in eV  
    hbgama : collision frequency in eV
    x : coordenada x 
    y : coordenada y 
    z : coordenada z
    xD : coordenada x del dipolo 
    yD : coordenada y del dipolo 
    zD : coordenada z del dipolo 
    zp : coordenada zp del plano
    Returns
    -------
    Gxx reflejado (superficie)
    con z hacia abajo (convencion del paper)
    numerico con la aprox QE (no hay alpha_z)
    """
    E = omegac*aux
    k0 = omegac #=omega/c
    # x_y = ky/k0
    n1 = epsi1*mu1
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1
    k1_3 = (k0*cte1)**3
    Rbarra =  k1*np.sqrt((x + xe)**2 + (y + ye)**2)
    z_dip_barra = k1*np.abs(np.abs(z) + 2*zp + np.abs(ze))
    phi = np.arctan2(np.abs(y + ye),np.abs(x + xe))
    
    expB = lambda u: np.exp(-u*z_dip_barra) 
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    
    rp_num = lambda u: epsi2*1j*u - epsi1*1j*u - cte1*cond*u**2
    rp_den = lambda u: epsi2*1j*u + epsi1*1j*u - cte1*cond*u**2
    rp = lambda u: rp_num(u)/rp_den(u)

    J0 = lambda u: special.jn(0,u*Rbarra) 
    J2 = lambda u: special.jn(2,u*Rbarra) 
    
    
    cota_inf = 1
    cota_sup = 1001*k1    
    cota_sup = 20
    
############ I0 5 ##################################
################### I0 6 ################ ###################
    Int06_B_function_re = lambda u: np.real(rp(u)*J0(u)*expB(u)*u**2)
    Int06_B_function_im = lambda u: np.imag(rp(u)*J0(u)*expB(u)*u**2)

    int06B_re,err = integrate.quad(Int06_B_function_re, cota_inf, cota_sup)
    int06B_im,err = integrate.quad(Int06_B_function_im, cota_inf, cota_sup)
################### I2 6 ################ ###################
    Int26_B_function_re = lambda u: np.real(rp(u)*J2(u)*expB(u)*u**2)
    Int26_B_function_im = lambda u: np.imag(rp(u)*J2(u)*expB(u)*u**2)

    int26B_re,err = integrate.quad(Int26_B_function_re, cota_inf, cota_sup)
    int26B_im,err = integrate.quad(Int26_B_function_im, cota_inf, cota_sup)

    
    int0_6 = int06B_re + 1j*int06B_im
    int2_6 = int26B_re + 1j*int26B_im
    
    cte = 0.5*k1_3
    cte2 = np.cos(2*phi)  
    
    return -int0_6*cte -int2_6*cte*cte2

#%%

def alpha_function_eff_ana_2dipo(omegac,epsi1,epsi2,hbmu,hbgama,int_v,zp,xD1,yD1,zD1,xD2,yD2,zD2,omega0,kappa_factor_omega0,kappa_r_factor):
    
    """
    Parameters
    ----------
    omegac : omega/c = k0 en 1/micrometros    
    epsi1 : epsilon del medio de arriba del plano
    x : coordenada x 
    y : coordenada y 
    z : coordenada z
    xD : coordenada x del dipolo 
    yD : coordenada y del dipolo 
    zD : coordenada z del dipolo 
    Returns
    -------
    Gxx direct (self interaction of the dipole)
    con z hacia abajo (convencion del paper)
    analitico luego de aplicar QE
    """
    
    G_ref = green_tensor_ref_rp_fresnel_QE(omegac,epsi1,epsi2,hbmu,hbgama,xD1,yD1,zD1,xD2,yD2,zD2,zp)
    G_dir = green_dir_ana(omegac,epsi1,xD1,yD1,zD1,xD2,yD2,zD2)
    
    alpha_eff_1 = alpha_function_eff_num(omegac,epsi1,epsi2,hbmu,hbgama,int_v,zp,xD1,yD1,zD1,omega0,kappa_factor_omega0,kappa_r_factor)
    alpha_eff_2 = alpha_function_eff_num(omegac,epsi1,epsi2,hbmu,hbgama,int_v,zp,xD2,yD2,zD2,omega0,kappa_factor_omega0,kappa_r_factor)
    
    Gtot = (G_dir + G_ref)**2
    
    rta = 1/alpha_eff_1 - alpha_eff_2*Gtot
    
    return rta  

#%%
    
def Efield_dir_ana_paper149(omegac,epsi1,int_v,b,x):
    """
    Parameters
    ----------
    omegac : omega/c = k0 en 1/micrometros    
    epsi1 : epsilon del medio de arriba del plano
    int_v : fraccion de la vel de la c que vale la velocidad uniforme del electron en x, si v = c/int
    b : electron position en z, b < 0 (dipolo en 0 y eje de z apunta hacia abajo)  
    Returns
    -------
    External field direct del paper 149
    """
    
    # x_y = ky/k0
#    n1 = epsi1*mu1
#    cte1 = np.sqrt(n1)
#    k1 = omegac*cte1
#    k1_2 = k1**2
    
    arg = np.abs(b)*omegac*int_v
    K0 = special.kn(0,arg)
    # cte = omega/(v*np.abs(b))
    
#    omega = omegac*c
    v = c/int_v

#    charge_e = 3.33564*1e-10
#    charge_e = 1
#    cte_aux = 1j*omegac*charge_e*int_v/np.abs(v)

    charge_e = 1
    cte_aux = 1j*2*omegac*charge_e*(int_v**2)/c

    expo = np.exp(1j*omegac*int_v*x)
    cte = 1 - int_v**(-2)

    arg = np.abs(b)*omegac*int_v
    K0 = special.kn(0,arg)
    
    return K0*cte_aux*cte*expo


#%%
    
def Efield_ref_semi_ana(omegac,epsi1,epsi2,hbmu,hbgama,x,y,z,int_v,b,zp):
    """
    Parameters
    ----------
    omegac : omega/c = k0 en 1/micrometros    
    epsi1 : epsilon del medio de arriba del plano
    int_v : fraccion de la vel de la c que vale la velocidad uniforme del electron en x, si v = c/int
    b : electron position en z, b < 0 (dipolo en 0 y eje de z apunta hacia abajo)  
    Returns
    -------
    External field direct del paper 149
    """
    E = omegac*aux    
    # x_y = ky/k0
    k0 = omegac
    n1 = epsi1*mu1
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1
    k1_2 = k1**2
    
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    
    alpha_p = k0*1j*(epsi1 + epsi2)/(cond*k1)
    Rp = 2*epsi1/(epsi1 + epsi2)
    
    z_dip_barra = k1*(np.abs(b) + 2*np.abs(zp) + np.abs(z))  #NO ESTOY SEGURA DE ESTA EXPONENCIAL EH , b
    expo = np.exp(-alpha_p*z_dip_barra)
    
    x_bar, y_bar = k1*x, k1*y
    
    Rref = lambda w: np.sqrt((x_bar + w)**2 + y_bar**2)
    
    phi = lambda w : np.arctan2(np.abs(y),np.abs(np.abs(w) + np.abs(x)))
    atan = lambda w : np.cos(2*phi(w))
    # cte = omega/(v*np.abs(b))
    
    J0 = lambda w: special.jn(0,alpha_p*Rref(w)) 
    J2 = lambda w: special.jn(2,alpha_p*Rref(w)) 
    
    exp_electron = lambda w: np.exp(1j*w*int_v) 
        
    function_int = lambda w: (J0(w) + J2(w)*atan(w))*exp_electron(w)
    
    omega = omegac*c
    v = c/int_v

#    charge_e = 3.33564*1e-10
    charge_e = 1
    
    alpha_p3 = alpha_p**3
    sng_v = np.sign(v)
    cte_final = alpha_p3*k1_2*sng_v*0.5*charge_e*expo*Rp/omega
    
    final_int,err = integrate.quad(function_int, -80*k1, 80*k1) 
    
    return final_int*cte_final

#%%
    
def p1_2dip(omegac,epsi1,epsi2,hbmu,hbgama,int_v,zp,b,xD1,yD1,zD1,xD2,yD2,zD2,omega0,kappa_factor_omega0,kappa_r_factor):
    
    """
    Parameters
    ----------
    omegac : omega/c = k0 en 1/micrometros    
    epsi1 : epsilon del medio de arriba del plano
    x : coordenada x 
    y : coordenada y 
    z : coordenada z
    xD : coordenada x del dipolo 
    yD : coordenada y del dipolo 
    zD : coordenada z del dipolo 
    Returns
    -------
    Gxx direct (self interaction of the dipole)
    con z hacia abajo (convencion del paper)
    analitico luego de aplicar QE
    """
    
    alpha_eff_eff_1 = alpha_function_eff_ana_2dipo(omegac,epsi1,epsi2,hbmu,hbgama,int_v,zp,xD1,yD1,zD1,xD2,yD2,zD2,omega0,kappa_factor_omega0,kappa_r_factor)
    
    Edir1 = Efield_dir_ana_paper149(omegac,epsi1,int_v,b,xD1)
    Eref1 = Efield_ref_semi_ana(omegac,epsi1,epsi2,hbmu,hbgama,xD1,yD1,zD1,int_v,b,zp)
    Etot1 = Edir1 + Eref1
    
    alpha_eff2 = alpha_function_eff_num(omegac,epsi1,epsi2,hbmu,hbgama,int_v,zp,xD2,yD2,zD2,omega0,kappa_factor_omega0,kappa_r_factor)
    
    Edir2 = Efield_dir_ana_paper149(omegac,epsi1,int_v,b,xD2)
    Eref2 = Efield_ref_semi_ana(omegac,epsi1,epsi2,hbmu,hbgama,xD2,yD2,zD2,int_v,b,zp)
    Etot2 = Edir2 + Eref2    
     
#    G_ref = green_ref_ana(omegac,epsi1,epsi2,hbmu,hbgama,xD1,yD1,zD1,xD2,yD2,zD2,zp)
    
    G_ref = green_tensor_ref_rp_fresnel_QE(omegac,epsi1,epsi2,hbmu,hbgama,xD1,yD1,zD1,xD2,yD2,zD2,zp)
    G_dir = green_dir_ana(omegac,epsi1,xD1,yD1,zD1,xD2,yD2,zD2)
    
    Gtot = G_dir + G_ref
    
    rta =  Etot1 + alpha_eff2*Etot2*Gtot
    
    return rta*alpha_eff_eff_1 

#%%
