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


def green_self_num(omegac,epsi1,epsi2,hbmu,hbgama,zp):
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

    rp_num = lambda u: epsi2*1j*u - epsi1*1j*u - cte1*cond*u**2
    rp_den = lambda u: epsi2*1j*u + epsi1*1j*u - cte1*cond*u**2
    rp = lambda u: rp_num(u)/rp_den(u)   
#    print(np.real(rp(10)))
    
#    rs_num = lambda u: 1j*u - 1j*u - cond/cte1
#    rs_den = lambda u: 1j*u + 1j*u + cond/cte1
#    rs = lambda u: rs_num(u)/rs_den(u)


    cota_sup1 = 400
    cota_inf1 = 0.01
    
####       
    cte_x = k1_3*0.5 #signo menos

    
    z_dip_barra_self = k1*2*zp  #NO ESTOY SEGURA DE ESTA EXPONENCIAL EH 
    expB_self = lambda u: np.exp(-u*z_dip_barra_self) 

    IntselfB_function_re_xx = lambda u: np.real((u**2)*rp(u)*expB_self(u))
    IntselfB_function_im_xx = lambda u: np.imag((u**2)*rp(u)*expB_self(u))

    intselfB_re_x,err = integrate.quad(IntselfB_function_re_xx, cota_inf1, cota_sup1)
    intselfB_im_x,err = integrate.quad(IntselfB_function_im_xx, cota_inf1, cota_sup1)
#    
#    print(intselfB_re_x)
    
    
    rtaself_x = (intselfB_re_x + 1j*intselfB_im_x)*cte_x
    rtaself_y = rtaself_x
    rtaself_z = 2*rtaself_x
    

    return rtaself_x, rtaself_y, rtaself_z




def green_self_num_integral_inside_light_cone(omegac,epsi1,epsi2,hbmu,hbgama,zp):
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

    rp_num = lambda u: epsi2*1j*u - epsi1*1j*u - cte1*cond*u**2
    rp_den = lambda u: epsi2*1j*u + epsi1*1j*u - cte1*cond*u**2
    rp = lambda u: rp_num(u)/rp_den(u)   
#    print(np.real(rp(10)))
    
#    rs_num = lambda u: 1j*u - 1j*u - cond/cte1
#    rs_den = lambda u: 1j*u + 1j*u + cond/cte1
#    rs = lambda u: rs_num(u)/rs_den(u)


    cota_sup1 = 1
    cota_inf1 = 0.01
    
####       
    cte_x = k1_3*0.5 #signo menos

    
    z_dip_barra_self = k1*2*zp  #NO ESTOY SEGURA DE ESTA EXPONENCIAL EH 
    expB_self = lambda u: np.exp(-u*z_dip_barra_self) 

    IntselfB_function_re_xx = lambda u: np.real((u**2)*rp(u)*expB_self(u))
    IntselfB_function_im_xx = lambda u: np.imag((u**2)*rp(u)*expB_self(u))

    intselfB_re_x,err = integrate.quad(IntselfB_function_re_xx, cota_inf1, cota_sup1)
    intselfB_im_x,err = integrate.quad(IntselfB_function_im_xx, cota_inf1, cota_sup1)
#    
#    print(intselfB_re_x)
    
    
    rtaself_x = (intselfB_re_x + 1j*intselfB_im_x)*cte_x
    rtaself_y = rtaself_x
    rtaself_z = 2*rtaself_x
    

    return rtaself_x, rtaself_y, rtaself_z
  






#%%
    

def green_self_pole_aprox(omegac,epsi1,epsi2,hbmu,hbgama,zp):
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
    Rp = 2*epsi1/(epsi1 + epsi2)
    alfa_p = 1j*(epsi1 + epsi2)/(cond*cte1)

    rp = lambda u: Rp*alfa_p/(u - alfa_p)   
#    print(np.real(rp(10)))
    
#    rs_num = lambda u: 1j*u - 1j*u - cond/cte1
#    rs_den = lambda u: 1j*u + 1j*u + cond/cte1
#    rs = lambda u: rs_num(u)/rs_den(u)

    cota_sup1 = 400
    cota_inf1 = 0.01
    
####       
    cte_x = k1_3*0.5 #signo menos

    
    z_dip_barra_self = k1*2*zp  #NO ESTOY SEGURA DE ESTA EXPONENCIAL EH 
    expB_self = lambda u: np.exp(-u*z_dip_barra_self) 

    IntselfB_function_re_xx = lambda u: u**2*np.real(rp(u))*expB_self(u)
    IntselfB_function_im_xx = lambda u: u**2*np.imag(rp(u))*expB_self(u)

    intselfB_re_x,err = integrate.quad(IntselfB_function_re_xx, cota_inf1, cota_sup1)
    intselfB_im_x,err = integrate.quad(IntselfB_function_im_xx, cota_inf1, cota_sup1)
    
#    print(intselfB_re_x)
    
    
    rtaself_x = (intselfB_re_x + 1j*intselfB_im_x)*cte_x
    rtaself_y = rtaself_x
    rtaself_z = 2*rtaself_x
    

    return rtaself_x, rtaself_y, rtaself_z

#%%

def green_self_ana(omegac,epsi1,epsi2,hbmu,hbgama,zp):     
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
#
    

    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    Rp = 2*epsi1/(epsi1 + epsi2)
    alfa_p = 1j*(epsi1 + epsi2)/(cond*cte1)
    kp = alfa_p*k1

    kp_3 = kp**3
    
#    kp = np.real(kp)   ### neglecting the imaginary part of kp para la formula de EELS (la del campo)
#    alfa_p = np.real(alfa_p)
      

    cte_x = 1j*np.pi*Rp*kp_3*0.5 #signo menos
    

    z_dip_barra_self = k1*2*zp  #NO ESTOY SEGURA DE ESTA EXPONENCIAL EH 
    expB_self = np.exp(-alfa_p*z_dip_barra_self) 

    
    rtaself_x = expB_self*cte_x
    rtaself_y = rtaself_x
    rtaself_z = 2*rtaself_x
    

    return rtaself_x, rtaself_y, rtaself_z


#%%



def green_self_ana2(omegac,epsi1,epsi2,hbmu,hbgama,zp):     
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
#
    

    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    Rp = 2*epsi1/(epsi1 + epsi2)
    alfa_p = 1j*(epsi1 + epsi2)/(cond*cte1)
    kp = alfa_p*k1

    kp_3 = kp**3
    
#    kp = np.real(kp)   ### neglecting the imaginary part of kp para la formula de EELS (la del campo)
#    alfa_p = np.real(alfa_p)
      


    

    arg = -2*kp*zp  #NO ESTOY SEGURA DE ESTA EXPONENCIAL EH 
    expo_self = np.exp(arg) 

    
    
    
    rtaself_x = Rp*kp*0.5*( (1 + 2*kp*zp)/(4*zp**2) + kp**2*special.exp1(arg)*expo_self  )
    
    

    rtaself_y = rtaself_x
    rtaself_z = 2*rtaself_x
    
    return rtaself_x, rtaself_y, rtaself_z







def green_self_ana3(omegac,epsi1,epsi2,hbmu,hbgama,zp):     
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
#
    

    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    Rp = 2*epsi1/(epsi1 + epsi2)
    alfa_p = 1j*(epsi1 + epsi2)/(cond*cte1)
    kp = alfa_p*k1

    kp_3 = kp**3
    
#    kp = np.real(kp)   ### neglecting the imaginary part of kp para la formula de EELS (la del campo)
#    alfa_p = np.real(alfa_p)
      


    

    arg = -2*kp*zp  #NO ESTOY SEGURA DE ESTA EXPONENCIAL EH 
    expo_self = np.exp(arg) 

    
    
    
    rtaself_x = Rp*kp*0.5*( (1 + 2*kp*zp)/(4*zp**2) + kp**2*np.pi*1j*expo_self  )
    
    

    rtaself_y = rtaself_x
    rtaself_z = 2*rtaself_x
    
    return rtaself_x, rtaself_y, rtaself_z





