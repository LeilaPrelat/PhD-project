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
path_constants =  path_basic.replace('/potential_field/potential_and_electric_field_with_dipole_moment_formula','')
#print('Importar modulos necesarios para este codigo')



try:
    sys.path.insert(1, path_constants)
    from dipole_moment import dipole_moment_anav2
except ModuleNotFoundError:
    print('dipole_moment.py no se encuentra en ' + path_basic) 
    
    
  # # # # # falta multiplicar el momento dipolar por e/(2*pi*v`)
    
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

limitt = 100

#%%

def EELS_ana_INT1(omegac,epsi1,epsi3,d_nano,int_v,b,zp,omega0,kappa_factor_omega0,kappa_r_factor):     
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


    
    px_v,py_v,pz_v = dipole_moment_anav2(omegac,epsi1,epsi3,d_nano,int_v,b,zp,omega0,kappa_factor_omega0,kappa_r_factor)
    
    v = c/int_v
    cte_dip_moment = alfac*c/(2*np.pi*v) ### 
    cte_dip_moment2 = alfac*int_v/(2*np.pi)

    v_2 = (c/int_v)**2
    cte_dip_moment2 = alfac*c/(pi*v_2)
    
    px_v,py_v,pz_v = px_v*cte_dip_moment2, py_v*cte_dip_moment2, pz_v*cte_dip_moment2

    d_micros = d_nano*1e-3
    Rp = hBn_Rp(E,epsi1,epsi3)
    lambda_p_v = hBn_lambda_p(E,epsi1,epsi3)*d_micros
    kp = 2*np.pi/lambda_p_v


    expo = np.exp(-kp*(-b + 2*zp))


    omega_v = omegac*int_v

    rta = px_v*2*pi*1j*Rp*kp*expo


    return rta/pi
#    else:
 #       return 0



def EELS_ana2_INT1(omegac,epsi1,epsi3,d_nano,int_v,b,zp,omega0,kappa_factor_omega0,kappa_r_factor):     
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


    
    px_v,py_v,pz_v = dipole_moment_anav2(omegac,epsi1,epsi2,hbmu,hbgama,int_v,b,zp,omega0,kappa_factor_omega0,kappa_r_factor)

    v = c/int_v
    cte_dip_moment = alfac*c/(2*np.pi*v) ### 
    cte_dip_moment2 = alfac*int_v/(2*np.pi)

    v_2 = (c/int_v)**2
    cte_dip_moment2 = alfac*c/(pi*v_2)
    
    px_v,py_v,pz_v = px_v*cte_dip_moment2, py_v*cte_dip_moment2, pz_v*cte_dip_moment2


    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    Rp = 2*epsi1/(epsi1 + epsi2)
    alfa_p = 1j*(epsi1 + epsi2)/(cond*cte1)
    kp = alfa_p*k1

    expo = np.exp(-kp*(-b + 2*zp))

    arg = (-b + 2*zp)*(omegac*int_v - kp) 
#    arg = -kp*(-b + 2*zp) 

    expo = np.exp(arg)

    omega_v = omegac*int_v
    
    rta = px_v*Rp*kp*special.exp1(arg)*expo
    
    return rta/pi


def EELS_num_INT1(omegac,epsi1,epsi3,d_nano,int_v,b,zp,omega0,kappa_factor_omega0,kappa_r_factor):     
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

    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    Rp = 2*epsi1/(epsi1 + epsi2)
    alfa_p = 1j*(epsi1 + epsi2)/(cond*cte1)
    
#    alfa_p = np.real(alfa_p)
    px_v,py_v,pz_v = dipole_moment_anav2(omegac,epsi1,epsi2,hbmu,hbgama,int_v,b,zp,omega0,kappa_factor_omega0,kappa_r_factor) 
    
    v = c/int_v
    cte_dip_moment = alfac*c/(2*np.pi*v) ### 
    cte_dip_moment2 = alfac*int_v/(2*np.pi)

    v_2 = (c/int_v)**2
    cte_dip_moment2 = alfac*c/(pi*v_2)
    
    px_v,py_v,pz_v = px_v*cte_dip_moment2, py_v*cte_dip_moment2, pz_v*cte_dip_moment2

    rp = lambda alpha_parallel : Rp*alfa_p/(alpha_parallel-alfa_p)
    
    rp_num = lambda u: epsi2*1j*u - epsi1*1j*u - cte1*cond*u**2
    rp_den = lambda u: epsi2*1j*u + epsi1*1j*u - cte1*cond*u**2
    rp = lambda u: rp_num(u)/rp_den(u)       
    
    
    exp_electron = lambda alpha_parallel: np.exp(-alpha_parallel*omegac*(np.abs(b) + 2*zp))

    cota_inf, cota_sup = int_v + 0.001, 1000
    

    term0_re_f = lambda u: np.real(rp(u)*exp_electron(u))
    term0_im_f = lambda u: np.imag(rp(u)*exp_electron(u))
    
    term0_int_re,err = integrate.quad(term0_re_f, cota_inf , cota_sup, limit = limitt) 
    term0_int_im,err = integrate.quad(term0_im_f, cota_inf , cota_sup, limit = limitt) 
   
    aux0 = term0_int_re + 1j*term0_int_im

    
    omega_v = omegac*int_v
        
    
    rta = px_v*aux0*k1

    return rta/pi


#%%
    

def EELS_ana_INT2(omegac,epsi1,epsi3,d_nano,int_v,b,zp,omega0,kappa_factor_omega0,kappa_r_factor):     
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

    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    Rp = 2*epsi1/(epsi1 + epsi2)
    alfa_p = 1j*(epsi1 + epsi2)/(cond)
    kp = alfa_p*omegac
    kp_3 = kp**3


    px_v,py_v,pz_v = dipole_moment_anav2(omegac,epsi1,epsi2,hbmu,hbgama,int_v,b,zp,omega0,kappa_factor_omega0,kappa_r_factor)


    v = c/int_v
    cte_dip_moment = alfac*c/(2*np.pi*v) ### 
    cte_dip_moment2 = alfac*int_v/(2*np.pi)

    v_2 = (c/int_v)**2
    cte_dip_moment2 = alfac*c/(pi*v_2)
   
    px_v,py_v,pz_v = px_v*cte_dip_moment2, py_v*cte_dip_moment2, pz_v*cte_dip_moment2

    expo = np.exp(-kp*(-b + 2*zp))    
    
    den1 = np.sqrt((omegac*int_v)**2 - kp**2)
    den2 = omegac*int_v + den1



    omega_v = omegac*int_v
    
    rta = -px_v*np.pi*1j*Rp*kp_3*expo/(den1*den2)

    return rta/pi


def EELS_num_INT2(omegac,epsi1,epsi3,d_nano,int_v,b,zp,omega0,kappa_factor_omega0,kappa_r_factor):     
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

    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    Rp = 2*epsi1/(epsi1 + epsi2)
    alfa_p = 1j*(epsi1 + epsi2)/(cond)


    k1 = omegac

    px_v,py_v,pz_v = dipole_moment_anav2(omegac,epsi1,epsi2,hbmu,hbgama,int_v,b,zp,omega0,kappa_factor_omega0,kappa_r_factor)


    v = c/int_v
    cte_dip_moment = alfac*c/(2*np.pi*v) ### 
    cte_dip_moment2 = alfac*int_v/(2*np.pi)

    v_2 = (c/int_v)**2
    cte_dip_moment2 = alfac*c/(pi*v_2)
   
    px_v,py_v,pz_v = px_v*cte_dip_moment2, py_v*cte_dip_moment2, pz_v*cte_dip_moment2

    rp = lambda alpha_parallel : Rp*alfa_p/(alpha_parallel-alfa_p)

    rp_num = lambda u: epsi2*1j*u - epsi1*1j*u - cond*u**2
    rp_den = lambda u: epsi2*1j*u + epsi1*1j*u - cond*u**2
    rp = lambda u: rp_num(u)/rp_den(u)   

    exp_electron = lambda alpha_parallel: np.exp(-alpha_parallel*omegac*(np.abs(b) + 2*zp))

    cota_inf, cota_sup = 0.001, int_v - 0.001
    

    term_sqr1 = lambda u : 1/np.sqrt(int_v**2 - u**2 )
    term_sqr2 = lambda u : 1/(int_v + term_sqr1(u))

    term0_re_f = lambda u: np.real(u**2*rp(u)*exp_electron(u)*term_sqr1(u)*term_sqr2(u))
    term0_im_f = lambda u: np.imag(u**2*rp(u)*exp_electron(u)*term_sqr1(u)*term_sqr2(u))
    
    term0_int_re,err = integrate.quad(term0_re_f, cota_inf , cota_sup, limit = limitt) 
    term0_int_im,err = integrate.quad(term0_im_f, cota_inf , cota_sup, limit = limitt) 
   
    aux0 = term0_int_re + 1j*term0_int_im
 

    omega_v = omegac*int_v
    
    rta =  -px_v*aux0*k1
    
    return rta/pi

#%%
    


def EELS_ana_INT3(omegac,epsi1,epsi3,d_nano,int_v,b,zp,omega0,kappa_factor_omega0,kappa_r_factor):     
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

    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    Rp = 2*epsi1/(epsi1 + epsi2)
    alfa_p = 1j*(epsi1 + epsi2)/(cond)
    kp = alfa_p*omegac
    kp_2 = kp**2


    px_v,py_v,pz_v = dipole_moment_anav2(omegac,epsi1,epsi2,hbmu,hbgama,int_v,b,zp,omega0,kappa_factor_omega0,kappa_r_factor)


    v = c/int_v
    cte_dip_moment = alfac*c/(2*np.pi*v) ### 
    cte_dip_moment2 = alfac*int_v/(2*np.pi)
    v_2 = (c/int_v)**2
    cte_dip_moment2 = alfac*c/(pi*v_2)
    
    px_v,py_v,pz_v = px_v*cte_dip_moment2, py_v*cte_dip_moment2, pz_v*cte_dip_moment2

    expo = np.exp(-kp*(-b + 2*zp))    
    
    den = np.sqrt( kp**2 - (omegac*int_v)**2 )
    

    omega_v = omegac*int_v
 
 
    rta =  pz_v*np.sign(b)*pi*1j*Rp*kp_2*expo/den

    return rta/pi

#%%

def EELS_num_INT3(omegac,epsi1,epsi3,d_nano,int_v,b,zp,omega0,kappa_factor_omega0,kappa_r_factor):     
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
    k_3 = omegac**3


    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    Rp = 2*epsi1/(epsi1 + epsi2)
    alfa_p = 1j*(epsi1 + epsi2)/(cond*cte1)
    

    
#    alfa_p = np.real(alfa_p)
    px_v,py_v,pz_v = dipole_moment_anav2(omegac,epsi1,epsi2,hbmu,hbgama,int_v,b,zp,omega0,kappa_factor_omega0,kappa_r_factor) 
    
    kp = alfa_p*k1

    v = c/int_v
    cte_dip_moment = alfac*c/(2*np.pi*v) ### 
    
    v_2 = (c/int_v)**2
    cte_dip_moment2 = alfac*c/(pi*v_2)
    
    
    px_v,py_v,pz_v = px_v*cte_dip_moment2, py_v*cte_dip_moment2, pz_v*cte_dip_moment2
   

    rp = lambda alpha_parallel : Rp*alfa_p/(alpha_parallel-alfa_p)

    rp_num = lambda u: epsi2*1j*u - epsi1*1j*u - cte1*cond*u**2
    rp_den = lambda u: epsi2*1j*u + epsi1*1j*u - cte1*cond*u**2
    rp = lambda u: rp_num(u)/rp_den(u)   

    exp_electron = lambda alpha_parallel: np.exp(-alpha_parallel*omegac*(np.abs(b) + 2*zp))

    cota_inf, cota_sup =  int_v + 0.001, 1000
    
    term_sqr = lambda u : 1/np.sqrt(u**2 - int_v**2  )

    term0_re_f = lambda u: np.real(term_sqr(u)*u*rp(u)*exp_electron(u))
    term0_im_f = lambda u: np.imag(term_sqr(u)*u*rp(u)*exp_electron(u))
    
    term0_int_re,err = integrate.quad(term0_re_f, cota_inf , cota_sup, limit = limitt) 
    term0_int_im,err = integrate.quad(term0_im_f, cota_inf , cota_sup, limit = limitt) 
   
    aux0 = term0_int_re + 1j*term0_int_im
    
#    omega_v = omegac*int_v

    rta = pz_v*np.sign(b)*aux0*k1
    
    return rta/pi


#%%



def EELS_ana_f_dir(omegac,epsi1,epsi3,d_nano,int_v,b,zp,omega0,kappa_factor_omega0,kappa_r_factor):     
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
    arg = omegac*int_v*np.abs(b)
    K0 = special.kn(0,arg)
    K1 = special.kn(1,arg)
    px_v,py_v,pz_v = dipole_moment_anav2(omegac,epsi1,epsi2,hbmu,hbgama,int_v,b,zp,omega0,kappa_factor_omega0,kappa_r_factor)


    v = c/int_v
    cte_dip_moment = alfac*c/(2*np.pi*v) ### 
    cte_dip_moment2 = alfac*c/pi
    
    px_v,py_v,pz_v = px_v*cte_dip_moment2, py_v*cte_dip_moment2, pz_v*cte_dip_moment2
    
#    omega = omegac*c
    v_3 = v**3
    omega = omegac*c
    cte_aux = omega/v_3

    
    rta =  (-px_v*K0 - 1j*pz_v*np.sign(b)*K1)*cte_aux
    return rta/pi


#%%

def EELS_parallel_f_dir(omegac,epsi1,epsi3,d_nano,int_v,b,zp,L):     
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

    n1 = epsi1*mu1
    cte1 = np.sqrt(n1)
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    Rp = 2*epsi1/(epsi1 + epsi2)
    alfa_p = 1j*(epsi1 + epsi2)/(cond*cte1)
    kp = alfa_p*omegac
    kp_2 = kp**2
    
    expo = np.exp(-2*kp*(zp-b))
    
#    omega = omegac*c    
    rta =  Rp*kp_2*expo/np.sqrt(kp_2 - (omegac*int_v)**2)
    
    v = c/int_v
    cte_final = alfac*c*L/(v**2) ## adimensional

    
    return 4*np.real(rta)*cte_final

#%%
