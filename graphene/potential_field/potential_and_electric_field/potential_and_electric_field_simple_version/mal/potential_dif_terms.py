#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila

comparar los terminos numericos y analiticos por separados

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
path_constants =  path_basic.replace('/potential_field/potential_and_electric_field/potential_and_electric_field_simple_version/potential_1_dipolo_analytical_vs_num','')
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

def electric_potential_ana1(omegac,epsi1,epsi2,hbmu,hbgama,x,y,z,zp,px,py,pz):     
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
    
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    Rp = 2*epsi1/(epsi1 + epsi2)
    alfa_p = 1j*(epsi1 + epsi2)/(cond*cte1)
    kp = alfa_p*k1
    kp_2 = kp**2 
    
    z_dip_barra = k1*(-z + 2*zp)        
    exp_electron = np.exp(-alfa_p*z_dip_barra)
    
    phi = np.arctan2(np.abs(y),np.abs(x))    
    R = np.sqrt(x**2 + y**2)
    term_px_py = px*np.cos(phi) + py*np.sin(phi)    
    # term_aux = np.abs(z)**2 + R**2
    J1 = special.jv(1,kp*R)
    # J0 = special.jv(0,kp*R)

    term1 = -term_px_py*Rp*2*np.pi*1j*kp_2*J1*exp_electron
    # term3 = -pz*np.sign(z)*(Rp*2*np.pi*1j*kp_2*J0*exp_electron )

    return term1

#%%

def electric_potential_ana2(omegac,epsi1,epsi2,hbmu,hbgama,x,y,z,zp,px,py,pz):     
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
    # n1 = epsi1*mu1
    # cte1 = np.sqrt(n1)
    k1 = omegac
#    k1_2 = (k0*cte1)**2
 #   n_v1 = int_v/cte1
    
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    Rp = 2*epsi1/(epsi1 + epsi2)
    alfa_p = 1j*(epsi1 + epsi2)/cond
    kp = alfa_p*k1
    kp_2 = kp**2 
    
    z_dip_barra = k1*(-z + 2*zp)        
    exp_electron = np.exp(-alfa_p*z_dip_barra)
    
    # phi = np.arctan2(np.abs(y),np.abs(x))    
    R = np.sqrt(x**2 + y**2)
    # term_px_py = px*np.cos(phi) + py*np.sin(phi)    
    # term_aux = np.abs(z)**2 + R**2
    # J1 = special.jv(1,kp*R)
    J0 = special.jv(0,kp*R)

    # term2 = -term_px_py*Rp*2*np.pi*1j*kp_2*J1*exp_electron
    term3 = -pz*np.sign(z)*Rp*2*np.pi*1j*kp_2*J0*exp_electron

    return term3 

#%%

def electric_potential_num1(omegac,epsi1,epsi2,hbmu,hbgama,x,y,z,zp,px,py,pz):     
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
    # n1 = epsi1*mu1
    # cte1 = np.sqrt(n1)
    k1 = omegac
#    k1_2 = (k0*cte1)**2
 #   n_v1 = int_v/cte1
    k1_2 = k1**2
    
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)

    
    phi = np.arctan2(np.abs(y),np.abs(x))    
    R = np.sqrt(x**2 + y**2)
    term_px_py = px*np.cos(phi) + py*np.sin(phi)    
    # term_aux = np.abs(z)**2 + R**2

    # term1 = -term_px_py*(np.abs(z)**2/(term_aux**(3/2)) - 1/(term_aux**(1/2)))/R
    # term3 = -pz*np.sign(z)*np.abs(z)/(term_aux**(3/2)) 
    
    cota_inf = 1
    cota_sup = 1001*k1
    
    z_dip_barra = 2*zp - z
    
    exp_electron_f = lambda u: np.exp(-u*k1*z_dip_barra)
    # J0 = lambda u : special.jv(0,u*R)
    J1 =  lambda u : special.jv(1,u*k1*R)


    rp_num = lambda u: epsi2*1j*u - epsi1*1j*u - cond*u**2
    rp_den = lambda u: epsi2*1j*u + epsi1*1j*u - cond*u**2
    rp = lambda u: rp_num(u)/rp_den(u)   
    
    term2_re_f = lambda u: np.real(J1(u)*u*rp(u)*exp_electron_f(u))
    term2_im_f = lambda u: np.imag(J1(u)*u*rp(u)*exp_electron_f(u))
    
    term2_int_re,err = integrate.quad(term2_re_f, cota_inf, cota_sup) 
    term2_int_im,err = integrate.quad(term2_im_f, cota_inf, cota_sup) 

    term2_int = term2_int_re + 1j*term2_int_im
    term2_final = -term_px_py*term2_int
    
    
    # term4_re_f = lambda u: np.real(J0(u)*u*rp(u)*exp_electron_f(u))
    # term4_im_f = lambda u: np.imag(J0(u)*u*rp(u)*exp_electron_f(u))
    
    # term4_int_re,err = integrate.quad(term4_re_f, cota_inf, cota_sup) 
    # term4_int_im,err = integrate.quad(term4_im_f, cota_inf, cota_sup)     
    
    # term4_int = term4_int_re + 1j*term4_int_im
    # term4_final = -pz*np.sign(z)*term4_int
    
    return term2_final*k1_2

#%%

def electric_potential_num2(omegac,epsi1,epsi2,hbmu,hbgama,x,y,z,zp,px,py,pz):     
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
    
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)

    
    # phi = np.arctan2(np.abs(y),np.abs(x))    
    R = np.sqrt(x**2 + y**2)
    # term_px_py = px*np.cos(phi) + py*np.sin(phi)    
    # term_aux = np.abs(z)**2 + R**2

    # term1 = -term_px_py*(np.abs(z)**2/(term_aux**(3/2)) - 1/(term_aux**(1/2)))/R
    # term3 = -pz*np.sign(z)*np.abs(z)/(term_aux**(3/2)) 
    
    cota_inf = 1
    cota_sup = 1001*k1   
    
    z_dip_barra = 2*zp - z
    
    exp_electron_f = lambda u: np.exp(-u*k1*z_dip_barra)
    J0 = lambda u : special.jv(0,u*k1*R)
    # J1 =  lambda u : special.jv(1,u*R)


    rp_num = lambda u: epsi2*1j*u - epsi1*1j*u - cond*u**2
    rp_den = lambda u: epsi2*1j*u + epsi1*1j*u - cond*u**2
    rp = lambda u: rp_num(u)/rp_den(u)   
    
    # term2_re_f = lambda u: np.real(J1(u)*u*rp(u)*exp_electron_f(u))
    # term2_im_f = lambda u: np.imag(J1(u)*u*rp(u)*exp_electron_f(u))
    
    # term2_int_re,err = integrate.quad(term2_re_f, cota_inf, cota_sup) 
    # term2_int_im,err = integrate.quad(term2_im_f, cota_inf, cota_sup) 

    # term2_int = term2_int_re + 1j*term2_int_im
    # term2_final = -term_px_py*term2_int
    
    
    term4_re_f = lambda u: np.real(J0(u)*u*rp(u)*exp_electron_f(u))
    term4_im_f = lambda u: np.imag(J0(u)*u*rp(u)*exp_electron_f(u))
    
    term4_int_re,err = integrate.quad(term4_re_f, cota_inf, cota_sup) 
    term4_int_im,err = integrate.quad(term4_im_f, cota_inf, cota_sup)     
    
    term4_int = term4_int_re + 1j*term4_int_im
    term4_final = -pz*np.sign(z)*term4_int
    
    return term4_final 

#%%