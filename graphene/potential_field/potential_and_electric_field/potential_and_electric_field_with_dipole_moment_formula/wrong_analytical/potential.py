#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila
The potential cannot depend on the position of the electron $b$ because that is inside the dipole moment
"""
import numpy as np
import sys
import os 
from scipy import special
from scipy import integrate

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/potential_field/potential/potential_and_electric_field_with_dipole_moment_formula','')
#print('Importar modulos necesarios para este codigo')

try:
    sys.path.insert(1, path_constants)
    from graphene_sigma import sigma_DL
except ModuleNotFoundError:
    print('graphene_sigma.py no se encuentra en ' + path_basic)


try:
    sys.path.insert(1, path_constants)
    from dipole_moment import dipole_moment_anav2,dipole_moment_num,dipole_moment_pole_aprox
except ModuleNotFoundError:
    print('dipole_moment.py no se encuentra en ' + path_basic)

try:
    sys.path.insert(1, path_constants)
    from constants import constantes
except ModuleNotFoundError:
    print('constants.py no se encuentra en ' + path_constants)

pi,hb,c,alfac,mu1,mu2 = constantes()
aux = c*hb


#%%

def electric_potential_ana(omegac,epsi1,epsi2,hbmu,hbgama,int_v,px,py,pz,b,zp,x,y,z,xD,yD,zD,omega0,kappa_factor_omega0,kappa_r_factor):     
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
    
    px_v,py_v,pz_v = dipole_moment_anav2(omegac,epsi1,epsi2,hbmu,hbgama,int_v,px,py,pz,b,zp,omega0,kappa_factor_omega0,kappa_r_factor)
    
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    Rp = 2*epsi1/(epsi1 + epsi2)
    alfa_p = 1j*(epsi1 + epsi2)/(cond*cte1)
    kp = alfa_p*k1
    kp_2 = kp**2 
    
    z_dip_barra = k1*(-z + 2*zp)        
    exp_electron = np.exp(-alfa_p*z_dip_barra)
    
    phi = np.arctan2(np.abs(y),np.abs(x))    
    R = np.sqrt(x**2 + y**2)
    term_px_py = px_v*np.cos(phi) + py_v*np.sin(phi)    
    term_aux = np.abs(z)**2 + R**2
    H1 = special.hankel1(1,kp*R)
    H0 = special.hankel1(0,kp*R)

    
    term1 = -term_px_py*(np.abs(z)**2/(term_aux**(3/2)) - 1/(term_aux**(1/2)))/R
    term2 = -term_px_py*Rp*np.pi*1j*kp_2*H1*exp_electron
    term3 = -pz_v*np.sign(z)*(Rp*np.pi*1j*kp_2*H0*exp_electron + np.abs(z)/(term_aux**(3/2)) )
    
    return  term1 + term2 + term3 
    
#%%

def electric_potential_num(omegac,epsi1,epsi2,hbmu,hbgama,int_v,px,py,pz,b,zp,x,y,z,xD,yD,zD,omega0,kappa_factor_omega0,kappa_r_factor):     
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
    
#    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
#    Rp = 2*epsi1/(epsi1 + epsi2)
#    alfa_p = 1j*(epsi1 + epsi2)/(cond*cte1)
#    kp = alfa_p*k1
#    kp_2 = kp**2 
    
    px_v,py_v,pz_v = dipole_moment_num(omegac,epsi1,epsi2,hbmu,hbgama,int_v,px,py,pz,b,zp,omega0,kappa_factor_omega0,kappa_r_factor)
    
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)

    
    phi = np.arctan2(np.abs(y),np.abs(x))    
    R = np.sqrt(x**2 + y**2)
    term_px_py = px_v*np.cos(phi) + py_v*np.sin(phi)    
    term_aux = np.abs(z)**2 + R**2

    term1 = -term_px_py*(np.abs(z)**2/(term_aux**(3/2)) - 1/(term_aux**(1/2)))/R
    term3 = -pz_v*np.sign(z)*np.abs(z)/(term_aux**(3/2)) 
    
    cota_inf = 1
    cota_sup = 1001   
    
    z_dip_barra = k1*(2*zp - z)
    
    exp_electron_f = lambda u: np.exp(-u*z_dip_barra)
    J0 = lambda u : special.jv(0,u*R)
    J1 =  lambda u : special.jv(1,u*R)


    rp_num = lambda u: epsi2*1j*u - epsi1*1j*u - cte1*cond*u**2
    rp_den = lambda u: epsi2*1j*u + epsi1*1j*u - cte1*cond*u**2
    rp = lambda u: rp_num(u)/rp_den(u)   
    
    term2_re_f = lambda u: np.real(J1(u)*u*rp(u)*exp_electron_f(u))
    term2_im_f = lambda u: np.imag(J1(u)*u*rp(u)*exp_electron_f(u))
    
    term2_int_re,err = integrate.quad(term2_re_f, cota_inf, cota_sup) 
    term2_int_im,err = integrate.quad(term2_im_f, cota_inf, cota_sup) 

    term2_int = term2_int_re + 1j*term2_int_im
    term2_final = -term_px_py*term2_int
    
    
    term4_re_f = lambda u: np.real(J0(u)*u*rp(u)*exp_electron_f(u))
    term4_im_f = lambda u: np.imag(J0(u)*u*rp(u)*exp_electron_f(u))
    
    term4_int_re,err = integrate.quad(term4_re_f, cota_inf, cota_sup) 
    term4_int_im,err = integrate.quad(term4_im_f, cota_inf, cota_sup)     
    
    term4_int = term4_int_re + 1j*term4_int_im
    term4_final = -pz_v*np.sign(z)*term4_int
    
    return term1 + term3 + term2_final + term4_final 
    
#%%
    

def electric_potential_pole_approx(omegac,epsi1,epsi2,hbmu,hbgama,int_v,px,py,pz,b,zp,x,y,z,xD,yD,zD,omega0,kappa_factor_omega0,kappa_r_factor):     
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
    
#    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
#    Rp = 2*epsi1/(epsi1 + epsi2)
#    alfa_p = 1j*(epsi1 + epsi2)/(cond*cte1)
#    kp = alfa_p*k1
#    kp_2 = kp**2 
    
    px_v,py_v,pz_v = dipole_moment_pole_aprox(omegac,epsi1,epsi2,hbmu,hbgama,int_v,px,py,pz,b,zp,omega0,kappa_factor_omega0,kappa_r_factor)
    
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)

    
    phi = np.arctan2(np.abs(y),np.abs(x))    
    R = np.sqrt(x**2 + y**2)
    term_px_py = px_v*np.cos(phi) + py_v*np.sin(phi)    
    term_aux = np.abs(z)**2 + R**2

    term1 = -term_px_py*(np.abs(z)**2/(term_aux**(3/2)) - 1/(term_aux**(1/2)))/R
    term3 = -pz_v*np.sign(z)*np.abs(z)/(term_aux**(3/2)) 
    
    cota_inf = 1
    cota_sup = 1001   
    
    z_dip_barra = k1*(2*zp - z)
    
    exp_electron_f = lambda u: np.exp(-u*z_dip_barra)
    J0 = lambda u : special.jv(0,u*R)
    J1 =  lambda u : special.jv(1,u*R)


    rp_num = lambda u: epsi2*1j*u - epsi1*1j*u - cte1*cond*u**2
    rp_den = lambda u: epsi2*1j*u + epsi1*1j*u - cte1*cond*u**2
    rp = lambda u: rp_num(u)/rp_den(u)   
    
    term2_re_f = lambda u: np.real(J1(u)*u*rp(u)*exp_electron_f(u))
    term2_im_f = lambda u: np.imag(J1(u)*u*rp(u)*exp_electron_f(u))
    
    term2_int_re,err = integrate.quad(term2_re_f, cota_inf, cota_sup) 
    term2_int_im,err = integrate.quad(term2_im_f, cota_inf, cota_sup) 

    term2_int = term2_int_re + 1j*term2_int_im
    term2_final = -term_px_py*term2_int
    
    
    term4_re_f = lambda u: np.real(J0(u)*u*rp(u)*exp_electron_f(u))
    term4_im_f = lambda u: np.imag(J0(u)*u*rp(u)*exp_electron_f(u))
    
    term4_int_re,err = integrate.quad(term4_re_f, cota_inf, cota_sup) 
    term4_int_im,err = integrate.quad(term4_im_f, cota_inf, cota_sup)     
    
    term4_int = term4_int_re + 1j*term4_int_im
    term4_final = -pz_v*np.sign(z)*term4_int
    
    return term1 + term3 + term2_final + term4_final 
    
#%%