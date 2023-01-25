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
    from hBn_PP import epsilon_x
except ModuleNotFoundError:
    print('hBn_PP.py no se encuentra en ' + path_constants)


try:
    sys.path.insert(1, path_constants)
    from dipole_moment import dipole_moment_ana_resonance_v1,dipole_moment_ana_resonance_v2,dipole_moment_pole_aprox_resonance_v1, dipole_moment_num_resonance,dipole_moment_pole_aprox_resonance_v2
except ModuleNotFoundError:
    print('green_self_image.py no se encuentra en ' + path_basic)


try:
    sys.path.insert(1, path_constants)
    from constants import constantes
except ModuleNotFoundError:
    print('constants.py no se encuentra en ' + path_constants)

pi,hb,c,alfac,mu1,mu2 = constantes()
aux = c*hb


#%%
    
#     rp = lambda u: alpha_parallel(u)/(alpha_parallel(u) - alfa_p)

def potential_ana_resonance_v1(omegac,epsi_silica,d_nano,int_v,b,zp,R,phi,z):     ### usando 
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
    px,py,pz en unidades de k*alfa_eff
    """

    E = omegac*aux
#    k0 = omegac #=omega/c
    # x_y = ky/k0

#    k1_2 = (k0*cte1)**2
 #   n_v1 = int_v/cte1

    px, py, pz  =  dipole_moment_ana_resonance_v1(omegac,epsi_silica,d_nano,int_v,b,zp)


#    charge_electron = 4.806e-10/c
#    cte_uni = int_v/(2*np.pi*c)
    epsi_x = epsilon_x(E)
    epsi_HBN_par = epsi_x
    
    d_micro = d_nano*1e-3
    alfa_p = epsi_silica(E)*2/(omegac*d_micro*(epsi_HBN_par-1))
    kp = alfa_p*omegac


      
    arg = kp*R
    H1 = special.hankel1(1,arg)
    H0 = special.hankel1(0,arg)
    
    expo = np.exp(-kp*(2*zp-z))
    
    
    term_p = px*np.cos(phi) + py*np.sin(phi)
    
    term_r_z = np.abs(z)**2 + R**2
    term1 = (np.abs(z)**2/(term_r_z**(3/2))  - term_r_z**(-1/2) )/R
    
    
    
    term2 = 2*np.pi*1j*H1*expo*(kp**2)
    
    term3 = 2*np.pi*1j*H0*expo*(kp**2)
    
    
    term4 = np.abs(z)/(term_r_z**(3/2))
    
    
    rta = term_p*(term1 + term2 ) - pz*np.sign(z)*(term3 + term4)
    
    cte_final = -2/(epsi_silica(E) + 1)
    
    return rta*cte_final


def potential_ana_resonance_v2(omegac,epsi_silica,d_nano,int_v,b,zp,R,phi,z):     ### usando 
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
    px,py,pz en unidades de k*alfa_eff
    """

    E = omegac*aux
#    k0 = omegac #=omega/c
    # x_y = ky/k0

#    k1_2 = (k0*cte1)**2
 #   n_v1 = int_v/cte1

    px, py, pz  =  dipole_moment_ana_resonance_v2(omegac,epsi_silica,d_nano,int_v,b,zp)


#    charge_electron = 4.806e-10/c
#    cte_uni = int_v/(2*np.pi*c)
    epsi_x = epsilon_x(E)
    epsi_HBN_par = epsi_x
    
    d_micro = d_nano*1e-3
    alfa_p = epsi_silica(E)*2/(omegac*d_micro*(epsi_HBN_par-1))
    kp = alfa_p*omegac


      
    arg = kp*R
    H1 = special.hankel1(1,arg)
    H0 = special.hankel1(0,arg)
    
    expo = np.exp(-kp*(2*zp-z))
    
    
    term_p = px*np.cos(phi) + py*np.sin(phi)
    
    term_r_z = np.abs(z)**2 + R**2
    term1 = (np.abs(z)**2/(term_r_z**(3/2))  - term_r_z**(-1/2) )/R
    
    
    
    term2 = 2*np.pi*1j*H1*expo*(kp**2)
    
    term3 = 2*np.pi*1j*H0*expo*(kp**2)
    
    
    term4 = np.abs(z)/(term_r_z**(3/2))
    
    
    rta = term_p*(term1 + term2 ) - pz*np.sign(z)*(term3 + term4)
    
    cte_final = -2/(epsi_silica(E) + 1)
    
    return rta*cte_final

#%%


def potential_num_resonance(omegac,epsi_silica,d_nano,int_v,b,zp,R,phi,z):     
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
    k0 = omegac #=omega/c
    # x_y = ky/k0
#    n1 = epsi1*mu1

    px, py, pz  =  dipole_moment_num_resonance(omegac,epsi_silica,d_nano,int_v,b,zp)


#    charge_electron = 4.806e-10/c
#    cte_uni = int_v/(2*np.pi*c)
    epsi_x = epsilon_x(E)
    epsi_HBN_par = epsi_x
    
    d_micro = d_nano*1e-3
#    alfa_p = epsi_silica(E)*2/(omegac*d_micro*(epsi_HBN_par-1))
#    kp = alfa_p*omegac

    sigma_2D = d_micro*(epsi_HBN_par - 1 )/(4*pi)  ## se cancela el i*omega del sigma 
    r_prima = lambda u: (1 - epsi_silica(E)/(2*pi*u*k0*sigma_2D) )**(-1)
    r = (1 - epsi_silica(E))/(1 + epsi_silica(E))
    exp_fresnel = lambda u:  np.exp(-2*u*k0*zp)  ## zp in micros
    r_fresnel = lambda u: r_prima(u)/(1 - r*r_prima(u)*exp_fresnel(u) )
      

    J1 = lambda u: special.jv(1,u*k0*R)
    J0 = lambda u: special.jv(0,u*k0*R)
    
    expo = lambda u :  np.exp(-u*k0*(2*zp-z))
    
    
    
    term_p = px*np.cos(phi) + py*np.sin(phi)
    
    term_r_z = np.abs(z)**2 + R**2
    term1 = (np.abs(z)**2/(term_r_z**(3/2))  - term_r_z**(-1/2) )/R
    
    
    k0_2 = k0**2
    
    term2_re = lambda u: (np.real(r_fresnel(u))*u*J1(u)*expo(u))*k0_2
    term2_im = lambda u: (np.imag(r_fresnel(u))*u*J1(u)*expo(u))*k0_2
    
    
    term3_re = lambda u: (np.real(r_fresnel(u))*u*J0(u)*expo(u))*k0_2  
    term3_im = lambda u: (np.imag(r_fresnel(u))*u*J0(u)*expo(u))*k0_2      
    
    
    cota_inf = 0.01/k0
    cota_sup = 600/k0
    
    INT2_re_y,err = integrate.quad(term2_re, cota_inf, cota_sup) 
    INT2_im_y,err = integrate.quad(term2_im, cota_inf, cota_sup) 
    
    INT3_re_y,err = integrate.quad(term3_re, cota_inf, cota_sup) 
    INT3_im_y,err = integrate.quad(term3_im, cota_inf, cota_sup) 
    
    
    
    INT2_y = INT2_re_y + 1j*INT2_im_y
    INT3_y = INT3_re_y + 1j*INT3_im_y        

    term4 = np.abs(z)/(term_r_z**(3/2))
    
    
    
#    term_extra = 2*np.pi*1j*Rp*kp*np.abs(kp)*expo/ky
    
    rta = term_p*(term1 + INT2_y ) - pz*np.sign(z)*(INT3_y + term4)
    
    cte_final = -2/(epsi_silica(E) + 1)
    
    return rta*cte_final


#%%
    
#     rp = lambda u: alpha_parallel(u)/(alpha_parallel(u) - alfa_p)
def potential_pole_aprox_resonance_v1(omegac,epsi_silica,d_nano,int_v,b,zp,R,phi,z):     
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
    k0 = omegac #=omega/c
    # x_y = ky/k0

    px, py, pz  =  dipole_moment_pole_aprox_resonance_v1(omegac,epsi_silica,d_nano,int_v,b,zp)


#    charge_electron = 4.806e-10/c
#    cte_uni = int_v/(2*np.pi*c)

    J1 = lambda u: special.jv(1,u*k0*R)
    J0 = lambda u: special.jv(0,u*k0*R)
    
    expo = lambda u :  np.exp(-u*k0*(2*zp-z))
    
    epsi_x = epsilon_x(E)
    epsi_HBN_par = epsi_x
    
    d_micro = d_nano*1e-3
    alfa_p = epsi_silica(E)*2/(omegac*d_micro*(epsi_HBN_par-1))    
    
    r_fresnel = lambda u: u/(u - alfa_p)
    
    term_p = px*np.cos(phi) + py*np.sin(phi)
    
    term_r_z = np.abs(z)**2 + R**2
    term1 = (np.abs(z)**2/(term_r_z**(3/2))  - term_r_z**(-1/2) )/R
    
    
    k0_2 = k0**2
    
    term2_re = lambda u: (np.real(r_fresnel(u))*u*J1(u)*expo(u))*k0_2
    term2_im = lambda u: (np.imag(r_fresnel(u))*u*J1(u)*expo(u))*k0_2
    
    
    term3_re = lambda u: (np.real(r_fresnel(u))*u*J0(u)*expo(u))*k0_2  
    term3_im = lambda u: (np.imag(r_fresnel(u))*u*J0(u)*expo(u))*k0_2      
    
    
    cota_inf = 0.01/k0
    cota_sup = 600/k0
    
    INT2_re_y,err = integrate.quad(term2_re, cota_inf, cota_sup) 
    INT2_im_y,err = integrate.quad(term2_im, cota_inf, cota_sup) 
    
    INT3_re_y,err = integrate.quad(term3_re, cota_inf, cota_sup) 
    INT3_im_y,err = integrate.quad(term3_im, cota_inf, cota_sup) 
    
    
    
    INT2_y = INT2_re_y + 1j*INT2_im_y
    INT3_y = INT3_re_y + 1j*INT3_im_y        

    term4 = np.abs(z)/(term_r_z**(3/2))
    
    
    
#    term_extra = 2*np.pi*1j*Rp*kp*np.abs(kp)*expo/ky
    
    rta = term_p*(term1 + INT2_y ) - pz*np.sign(z)*(INT3_y + term4)
    
    cte_final = -2/(epsi_silica(E) + 1)
    
    return rta*cte_final



#     rp = lambda u: alfa_p/(alpha_parallel(u) - alfa_p)

def potential_pole_aprox_resonance_v2(omegac,epsi_silica,d_nano,int_v,b,zp,R,phi,z):     
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
    k0 = omegac #=omega/c
    # x_y = ky/k0
#    n1 = epsi1*mu1
#    cte1 = np.sqrt(n1)
    # k1 = omegac*cte1
#    k1_2 = (k0*cte1)**2
 #   n_v1 = int_v/cte1
    

    px, py, pz  =  dipole_moment_pole_aprox_resonance_v2(omegac,epsi_silica,d_nano,int_v,b,zp)


#    charge_electron = 4.806e-10/c
#    cte_uni = int_v/(2*np.pi*c)

    J1 = lambda u: special.jv(1,u*k0*R)
    J0 = lambda u: special.jv(0,u*k0*R)
    
    expo = lambda u :  np.exp(-u*k0*(2*zp-z))
    
    epsi_x = epsilon_x(E)
    epsi_HBN_par = epsi_x
    
    d_micro = d_nano*1e-3
    alfa_p = epsi_silica(E)*2/(omegac*d_micro*(epsi_HBN_par-1))    
    
    r_fresnel = lambda u: alfa_p/(u - alfa_p)
    
    term_p = px*np.cos(phi) + py*np.sin(phi)
    
    term_r_z = np.abs(z)**2 + R**2
    term1 = (np.abs(z)**2/(term_r_z**(3/2))  - term_r_z**(-1/2) )/R
    
    
    k0_2 = k0**2
    
    term2_re = lambda u: (np.real(r_fresnel(u))*u*J1(u)*expo(u))*k0_2
    term2_im = lambda u: (np.imag(r_fresnel(u))*u*J1(u)*expo(u))*k0_2
    
    
    term3_re = lambda u: (np.real(r_fresnel(u))*u*J0(u)*expo(u))*k0_2  
    term3_im = lambda u: (np.imag(r_fresnel(u))*u*J0(u)*expo(u))*k0_2      
    
    
    cota_inf = 0.01/k0
    cota_sup = 600/k0
    
    INT2_re_y,err = integrate.quad(term2_re, cota_inf, cota_sup) 
    INT2_im_y,err = integrate.quad(term2_im, cota_inf, cota_sup) 
    
    INT3_re_y,err = integrate.quad(term3_re, cota_inf, cota_sup) 
    INT3_im_y,err = integrate.quad(term3_im, cota_inf, cota_sup) 
    
    
    
    INT2_y = INT2_re_y + 1j*INT2_im_y
    INT3_y = INT3_re_y + 1j*INT3_im_y        

    term4 = np.abs(z)/(term_r_z**(3/2))
    
    
    
#    term_extra = 2*np.pi*1j*Rp*kp*np.abs(kp)*expo/ky
    
    rta = term_p*(term1 + INT2_y ) - pz*np.sign(z)*(INT3_y + term4)
    
    cte_final = -2/(epsi_silica(E) + 1)
    
    return rta*cte_final


#%%
    
