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
    from green_self_image_dif_sign import green_self_pole_aprox_v1,green_self_pole_aprox_v2, green_self_ana_v2, green_self_ana_v1, green_self_num,green_self_num_integral_inside_light_cone, green_self_ana2_perfect_dipole, green_self_num_integral_inside_light_cone_perfect_dipole
except ModuleNotFoundError:
    print('green_self_image_dif_sig.py no se encuentra en ' + path_constants)


try:
    sys.path.insert(1, path_constants)
    from hBn_PP import epsilon_x, epsilon_z,hBn_lambda_p,hBn_Rp
except ModuleNotFoundError:
    print('hBn_PP.py no se encuentra en ' + path_constants)


try:
    sys.path.insert(1, path_constants)
    from eff_polarizability_disk import polarizability_perp, polarizability_parallel
except ModuleNotFoundError:
    print('hBn_PP.py no se encuentra en ' + path_constants)



try:
    sys.path.insert(1, path_constants)
    from constants import constantes
except ModuleNotFoundError:
    print('constants.py no se encuentra en ' + path_constants)

pi,hb,c,alfac,mu1,mu2 = constantes()
aux = c*hb


#%%

    
#     rp = lambda u: alpha_parallel(u)/(alpha_parallel(u) - alfa_p)

def dipole_moment_ana_resonance_v1(omegac,epsi_silica,d_nano_film,d_thickness_disk_nano,D_disk_nano,int_v,b,zp):     ### usando 
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

    rtaself_x, rtaself_y, rtaself_z  =  green_self_ana_v1(omegac,epsi_silica,d_nano_film,zp)
    

    alfa_parallel = polarizability_parallel(E,d_thickness_disk_nano,D_disk_nano,epsi_silica)
    alfa_perp = polarizability_perp(E,d_thickness_disk_nano,D_disk_nano,epsi_silica)

    
    alffa_eff_x = (1/alfa_parallel -  rtaself_x)**(-1)
    alffa_eff_y = alffa_eff_x
    alffa_eff_z = (1/alfa_perp -  rtaself_z)**(-1)

#    charge_electron = 4.806e-10/c
#    cte_uni = int_v/(2*np.pi*c)
    epsi_x = epsilon_x(E)
    epsi_HBN_par = epsi_x
    
    d_micro = d_nano_film*1e-3
    alfa_p = epsi_silica(E)*2/(omegac*d_micro*(epsi_HBN_par-1))
    kp = alfa_p*omegac

#    d_micros = d_nano*1e-3
#    lambda_p_v = hBn_lambda_p(E,epsi_silica(E),epsi_silica(E))*d_micros
#    kp = 2*np.pi/lambda_p_v
#    alfa_p = kp/omegac 

      
    arg = np.abs(b)*omegac*int_v
    K1 = special.kn(1,arg)
    K0 = special.kn(0,arg)
    
    kx = omegac*int_v
#    expo = np.exp(-np.sqrt(kx**2 + kp**2)*(np.abs(b) + 2*zp))
    expo = np.exp(-kp*(np.abs(b) + 2*zp))
    
    ky = np.sqrt(kp**2 - kx**2)
    kp_2 = np.sqrt(kp**2)
    term_kp = kp_2 + kp
    
#    term_kp = 2*kp
#    kp_2 = kp
#    term_extra = 2*np.pi*1j*Rp*kp*np.abs(kp)*expo/ky
    
    
    px = alffa_eff_x*1j*omegac*int_v*(K0 + np.pi*1j*term_kp*expo/ky)
    
    py = alffa_eff_y*1j*(2*1j*omegac*int_v*K1 + np.pi*1j*term_kp*expo)
    
    pz = alffa_eff_z*(-omegac*int_v*K1 - np.pi*1j*kp_2*term_kp*expo/ky )
    
    
    return px, py, pz


def dipole_moment_anav1_for_decay_rate_resonance(omegac,epsi_silica,d_nano_film,d_thickness_disk_nano,D_disk_nano,int_v,b,zp):     ### usando 
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

    rtaself_x1, rtaself_y1, rtaself_z1  =  green_self_ana_v2(omegac,epsi_silica,d_nano_film,zp)
    rtaself_x2, rtaself_y2, rtaself_z2  =  green_self_num_integral_inside_light_cone(omegac,epsi_silica,d_nano_film,zp)

    rtaself_x, rtaself_y, rtaself_z  =  rtaself_x1 - rtaself_x2, rtaself_y1 - rtaself_y2, rtaself_z1 - rtaself_z2
    
#
#    alfa_parallel = polarizability_parallel(E,d_thickness_disk_nano,D_disk_nano,epsi_silica)
#    alfa_perp = polarizability_perp(E,d_thickness_disk_nano,D_disk_nano,epsi_silica)
#
#    
#    alffa_eff_x = (1/alfa_parallel -  rtaself_x)**(-1)
#    alffa_eff_y = alffa_eff_x
#    alffa_eff_z = (1/alfa_perp -  rtaself_z)**(-1)

    alffa_eff_x = 1j*(2*omegac**3/(3*epsi_silica(E)) +  np.imag(rtaself_x))**(-1)
    alffa_eff_y = 1j*(2*omegac**3/(3*epsi_silica(E)) +  np.imag(rtaself_y))**(-1)
    alffa_eff_z = 1j*(2*omegac**3/(3*epsi_silica(E)) +  np.imag(rtaself_z))**(-1)

#    charge_electron = 4.806e-10/c
#    cte_uni = int_v/(2*np.pi*c)
    epsi_x = epsilon_x(E)
    epsi_HBN_par = epsi_x
    
    d_micro = d_nano_film*1e-3
    alfa_p = epsi_silica(E)*2/(omegac*d_micro*(epsi_HBN_par-1))
    kp = alfa_p*omegac

    d_micros = d_nano_film*1e-3
    lambda_p_v = hBn_lambda_p(E,epsi_silica(E),epsi_silica(E))*d_micros
    kp = 2*np.pi/lambda_p_v
    alfa_p = kp/omegac 

      
    arg = np.abs(b)*omegac*int_v
    K1 = special.kn(1,arg)
    K0 = special.kn(0,arg)
    
    kx = omegac*int_v
#    expo = np.exp(-np.sqrt(kx**2 + kp**2)*(np.abs(b) + 2*zp))
    expo = np.exp(-kp*(np.abs(b) + 2*zp))
    
    ky = np.sqrt(kp**2 - kx**2)
    kp_2 = np.sqrt(kp**2)
    term_kp = kp_2 + kp
    
    term_kp = 2*kp
    kp_2 = kp
    
#    term_extra = 2*np.pi*1j*Rp*kp*np.abs(kp)*expo/ky
    
    
    px = alffa_eff_x*1j*omegac*int_v*(K0 + np.pi*1j*term_kp*expo/ky)
    
    py = alffa_eff_y*1j*(2*1j*omegac*int_v*K1 + np.pi*1j*term_kp*expo)
    
    pz = alffa_eff_z*(-omegac*int_v*K1 - np.pi*1j*kp_2*term_kp*expo/ky )
    
    
    return px, py, pz

    

def dipole_moment_ana_perfect_dipole_resonance(omegac,epsi_silica,d_nano_film,d_thickness_disk_nano,D_disk_nano,int_v,b,zp):     
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

    rtaself_x1, rtaself_y1, rtaself_z1  =  green_self_ana2_perfect_dipole(omegac,epsi_silica,d_nano_film,zp)
    rtaself_x2, rtaself_y2, rtaself_z2  =  green_self_num_integral_inside_light_cone_perfect_dipole(omegac,epsi_silica,d_nano_film,zp)

    rtaself_x, rtaself_y, rtaself_z  =  rtaself_x1 - rtaself_x2, rtaself_y1 - rtaself_y2, rtaself_z1 - rtaself_z2
    

    alfa_parallel = polarizability_parallel(E,d_thickness_disk_nano,D_disk_nano,epsi_silica)
    alfa_perp = polarizability_perp(E,d_thickness_disk_nano,D_disk_nano,epsi_silica)

    
    alffa_eff_x = (1/alfa_parallel -  rtaself_x)**(-1)
    alffa_eff_y = alffa_eff_x
    alffa_eff_z = (1/alfa_perp -  rtaself_z)**(-1)

    alffa_eff_x = 1j*(2*omegac**3/(3*epsi_silica(E)) +  np.imag(rtaself_x))**(-1)
    alffa_eff_y = 1j*(2*omegac**3/(3*epsi_silica(E)) +  np.imag(rtaself_y))**(-1)
    alffa_eff_z = 1j*(2*omegac**3/(3*epsi_silica(E)) +  np.imag(rtaself_z))**(-1)
    
    
#    epsi1 = 1
    d_micros = d_nano_film*1e-3
    Rp = hBn_Rp(E,epsi_silica(E),epsi_silica(E))
    lambda_p_v = hBn_lambda_p(E,epsi_silica(E),epsi_silica(E))*d_micros
    kp = 2*np.pi/lambda_p_v
    alfa_p = kp/omegac 
            
    arg = np.abs(b)*omegac*int_v
    K1 = special.kn(1,arg)
    K0 = special.kn(0,arg)
    
    kx = omegac*int_v
#    expo = np.exp(-np.sqrt(kx**2 + kp**2)*(np.abs(b) + 2*zp))
    expo = np.exp(-kp*(np.abs(b) + 2*zp))
    
#    den = np.sqrt(kx**2 + kp**2) - kp
    ky = np.sqrt(kp**2 - kx**2)
    kp_2 = np.sqrt(kp**2)
    term_kp = 1 + kp/kp_2
    
#    term_extra = 2*np.pi*1j*Rp*kp*np.abs(kp)*expo/ky
    
    
    px = alffa_eff_x*1j*omegac*int_v*(K0 - 2*np.pi*1j*Rp*kp*expo/ky)
    
    py = alffa_eff_y*1j*(2*1j*omegac*int_v*K1 - 2*np.pi*1j*Rp*kp*expo)
    
    pz = alffa_eff_z*(-omegac*int_v*K1 + 2*np.pi*1j*Rp*(kp**2)*expo/ky )
    
    return px, py, pz


#     rp = lambda u: alfa_p/(alpha_parallel(u) - alfa_p)

def dipole_moment_ana_resonance_v2(omegac,epsi_silica,d_nano_film,d_thickness_disk_nano,D_disk_nano,int_v,b,zp):     
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

    rtaself_x, rtaself_y, rtaself_z  =  green_self_ana_v2(omegac,epsi_silica,d_nano_film,zp)

    alfa_parallel = polarizability_parallel(E,d_thickness_disk_nano,D_disk_nano,epsi_silica)
    alfa_perp = polarizability_perp(E,d_thickness_disk_nano,D_disk_nano,epsi_silica)

    
    alffa_eff_x = (1/alfa_parallel -  rtaself_x)**(-1)
    alffa_eff_y = alffa_eff_x
    alffa_eff_z = (1/alfa_perp -  rtaself_z)**(-1)


#    charge_electron = 4.806e-10/c
#    cte_uni = int_v/(2*np.pi*c)
    epsi_x = epsilon_x(E)
    epsi_HBN_par = epsi_x
    
    d_micro = d_nano_film*1e-3
    alfa_p = epsi_silica(E)*2/(omegac*d_micro*(epsi_HBN_par-1))
    kp = alfa_p*omegac


      
    arg = np.abs(b)*omegac*int_v
    K1 = special.kn(1,arg)
    K0 = special.kn(0,arg)
    
    kx = omegac*int_v
#    expo = np.exp(-np.sqrt(kx**2 + kp**2)*(np.abs(b) + 2*zp))
    expo = np.exp(-kp*(np.abs(b) + 2*zp))
    
    ky = np.sqrt(kp**2 - kx**2)
    kp_2 = np.sqrt(kp**2)
    term_kp = 1 + kp/kp_2
    
#    term_extra = 2*np.pi*1j*Rp*kp*np.abs(kp)*expo/ky
    
    Rp = 1
    px = alffa_eff_x*1j*omegac*int_v*(K0 + 2*np.pi*1j*Rp*kp*expo/ky)
    
    py = alffa_eff_y*1j*(2*1j*omegac*int_v*K1 + 2*np.pi*1j*Rp*kp*expo)
    
    pz = alffa_eff_z*(-omegac*int_v*K1 - 2*np.pi*1j*Rp*(kp**2)*expo/ky )
    
    return px, py, pz




def dipole_moment_anav2_for_decay_rate_resonance(omegac,epsi_silica,d_nano_film,d_thickness_disk_nano,D_disk_nano,int_v,b,zp):     
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



    rtaself_x1, rtaself_y1, rtaself_z1  =  green_self_ana_v2(omegac,epsi_silica,d_nano_film,zp)
    rtaself_x2, rtaself_y2, rtaself_z2  =  green_self_num_integral_inside_light_cone(omegac,epsi_silica,d_nano_film,zp)

    rtaself_x, rtaself_y, rtaself_z  =  rtaself_x1 - rtaself_x2, rtaself_y1 - rtaself_y2, rtaself_z1 - rtaself_z2
#    rtaself_x, rtaself_y, rtaself_z  =  rtaself_x1 , rtaself_y1 , rtaself_z1 
    
    alfa_parallel = polarizability_parallel(E,d_thickness_disk_nano,D_disk_nano,epsi_silica)
    alfa_perp = polarizability_perp(E,d_thickness_disk_nano,D_disk_nano,epsi_silica)

    
    alffa_eff_x = (1/alfa_parallel -  rtaself_x)**(-1)
    alffa_eff_y = alffa_eff_x
    alffa_eff_z = (1/alfa_perp -  rtaself_z)**(-1)

#    charge_electron = 4.806e-10/c
#    cte_uni = int_v/(2*np.pi*c)
    epsi_x = epsilon_x(E)
    epsi_HBN_par = epsi_x
    
    d_micro = d_nano_film*1e-3
    alfa_p = epsi_silica(E)*2/(omegac*d_micro*(epsi_HBN_par-1))
    kp = alfa_p*omegac


      
    arg = np.abs(b)*omegac*int_v
    K1 = special.kn(1,arg)
    K0 = special.kn(0,arg)
    
    kx = omegac*int_v
#    expo = np.exp(-np.sqrt(kx**2 + kp**2)*(np.abs(b) + 2*zp))
    expo = np.exp(-kp*(np.abs(b) + 2*zp))
    
    ky = np.sqrt(kp**2 - kx**2)
    kp_2 = np.sqrt(kp**2)
    term_kp = kp_2 + kp
    
#    term_extra = 2*np.pi*1j*Rp*kp*np.abs(kp)*expo/ky
    
    
    px = alffa_eff_x*1j*omegac*int_v*(K0 + np.pi*1j*term_kp*expo/ky)
    
    py = alffa_eff_y*1j*(2*1j*omegac*int_v*K1 + np.pi*1j*term_kp*expo)
    
    pz = alffa_eff_z*(-omegac*int_v*K1 - np.pi*1j*kp_2*term_kp*expo/ky )
    
    
    return px, py, pz



def dipole_moment_ana_for_decay_rate_resonance_dir(omegac,int_v,b,zp):     
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


#    k1_2 = (k0*cte1)**2
      
    arg = np.abs(b)*omegac*int_v
    K1 = special.kn(1,arg)
    K0 = special.kn(0,arg)

    
#    term_extra = 2*np.pi*1j*Rp*kp*np.abs(kp)*expo/ky
    
    
    px = 1j*omegac*int_v*K0 
    
    py = -2*omegac*int_v*K1 
                         
    pz = -omegac*int_v*K1 
    
    return px, py, pz




#%%


def dipole_moment_num_resonance(omegac,epsi_silica,d_nano_film,d_thickness_disk_nano,D_disk_nano,int_v,b,zp):     
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
    

    rtaself_x, rtaself_y, rtaself_z  =  green_self_num(omegac,epsi_silica,d_nano_film,zp)
    
    alfa_parallel = polarizability_parallel(E,d_thickness_disk_nano,D_disk_nano,epsi_silica)
    alfa_perp = polarizability_perp(E,d_thickness_disk_nano,D_disk_nano,epsi_silica)

    
    alffa_eff_x = (1/alfa_parallel -  rtaself_x)**(-1)
    alffa_eff_y = alffa_eff_x
    alffa_eff_z = (1/alfa_perp -  rtaself_z)**(-1)

    
#    charge_electron = alffa_eff*4.806e-10*int_v/(2*np.pi)
    
    cota_inf = 0.01*k0
    cota_sup = 600*k0   
    
    alpha_x = int_v    
    
#    term6 = np.sign(z)*pz*Rp*kp_2*J0*exp_electron
    alpha_parallel = lambda u: np.sqrt(alpha_x**2 + u**2)
    
    

    epsi_x = epsilon_x(E)
    epsi_HBN_par = epsi_x

    r = (1 - epsi_silica(E))/(1 + epsi_silica(E))

    expB_self = lambda u: np.exp(-alpha_parallel(u)*omegac*2*zp) 
    
    d_micro = d_nano_film*1e-3
    alfa_p = epsi_silica(E)*2/(omegac*d_micro*(epsi_HBN_par-1))


    rp = lambda u: -alpha_parallel(u)/(alpha_parallel(u)*(1-r*expB_self(u)) - alfa_p)
 
      
    expo = lambda u: np.exp(-np.sqrt(alpha_x**2 + u**2)*k0*(2*zp + np.abs(b)))
    
    int_f_re_x = lambda u: np.real(rp(alpha_parallel(u))*expo(u)/alpha_parallel(u))
    int_f_im_x = lambda u: np.imag(rp(alpha_parallel(u))*expo(u)/alpha_parallel(u))
    
    INT_re_x,err = integrate.quad(int_f_re_x, cota_inf, cota_sup) 
    INT_im_x,err = integrate.quad(int_f_im_x, cota_inf, cota_sup) 
    
    INT_x = INT_re_x + 1j*INT_im_x
    
    
    int_f_re_y = lambda u: np.real(rp(alpha_parallel(u))*expo(u)*u/alpha_parallel(u))
    int_f_im_y = lambda u: np.imag(rp(alpha_parallel(u))*expo(u)*u/alpha_parallel(u))
    
    INT_re_y,err = integrate.quad(int_f_re_y, cota_inf, cota_sup) 
    INT_im_y,err = integrate.quad(int_f_im_y, cota_inf, cota_sup) 
    
    INT_y = INT_re_y + 1j*INT_im_y
        


    int_f_re_z = lambda u: np.real(rp(alpha_parallel(u))*expo(u))
    int_f_im_z = lambda u: np.imag(rp(alpha_parallel(u))*expo(u))
    
    INT_re_z,err = integrate.quad(int_f_re_z, cota_inf, cota_sup) 
    INT_im_z,err = integrate.quad(int_f_im_z, cota_inf, cota_sup) 
    
    INT_z = INT_re_z + 1j*INT_im_z


    arg = np.abs(b)*omegac*int_v
    K1 = special.kn(1,arg)
    K0 = special.kn(0,arg)    
    
    
    px = alffa_eff_x*1j*omegac*int_v*(K0 - INT_x)
    
    py = alffa_eff_y*1j*(2*1j*omegac*int_v*K1 - k0*INT_y)
    
    pz = alffa_eff_z*(-omegac*int_v*K1 + k0*INT_z )
    
    return px, py, pz


#%%
    
#     rp = lambda u: alpha_parallel(u)/(alpha_parallel(u) - alfa_p)
def dipole_moment_pole_aprox_resonance_v1(omegac,epsi_silica,d_nano_film,d_thickness_disk_nano,D_disk_nano,int_v,b,zp):     
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
    


    rtaself_x, rtaself_y, rtaself_z  =  green_self_pole_aprox_v1(omegac,epsi_silica,d_nano_film,zp)

    alfa_parallel = polarizability_parallel(E,d_thickness_disk_nano,D_disk_nano,epsi_silica)
    alfa_perp = polarizability_perp(E,d_thickness_disk_nano,D_disk_nano,epsi_silica)

    
    alffa_eff_x = (1/alfa_parallel -  rtaself_x)**(-1)
    alffa_eff_y = alffa_eff_x
    alffa_eff_z = (1/alfa_perp -  rtaself_z)**(-1)


    epsi_x = epsilon_x(E)
    epsi_HBN_par = epsi_x
    
    d_micro = d_nano_film*1e-3
    alfa_p = epsi_silica(E)*2/(omegac*d_micro*(epsi_HBN_par-1))

    cota_inf = 0.01/k0
    cota_sup = 600/k0   
    
    alpha_x = int_v    
    
#    term6 = np.sign(z)*pz*Rp*kp_2*J0*exp_electron
    alpha_parallel = lambda u: np.sqrt(alpha_x**2 + u**2)

    rp = lambda u: -alpha_parallel(u)/(alpha_parallel(u) - alfa_p)
      
    expo = lambda u: np.exp(-np.sqrt(alpha_x**2 + u**2)*k0*(2*zp + np.abs(b)))
    
    int_f_re_x = lambda u: np.real(rp(u)*expo(u)/alpha_parallel(u))
    int_f_im_x = lambda u: np.imag(rp(u)*expo(u)/alpha_parallel(u))
    
    INT_re_x,err = integrate.quad(int_f_re_x, cota_inf, cota_sup) 
    INT_im_x,err = integrate.quad(int_f_im_x, cota_inf, cota_sup) 
    
    INT_x = INT_re_x + 1j*INT_im_x
    
    
    int_f_re_y = lambda u: np.real(rp(u)*expo(u)*u/alpha_parallel(u))
    int_f_im_y = lambda u: np.imag(rp(u)*expo(u)*u/alpha_parallel(u))
    
    INT_re_y,err = integrate.quad(int_f_re_y, cota_inf, cota_sup) 
    INT_im_y,err = integrate.quad(int_f_im_y, cota_inf, cota_sup) 
    
    INT_y = INT_re_y + 1j*INT_im_y
        


    int_f_re_z = lambda u: np.real(rp(u)*expo(u))
    int_f_im_z = lambda u: np.imag(rp(u)*expo(u))
    
    INT_re_z,err = integrate.quad(int_f_re_z, cota_inf, cota_sup) 
    INT_im_z,err = integrate.quad(int_f_im_z, cota_inf, cota_sup) 
    
    INT_z = INT_re_z + 1j*INT_im_z


    arg = np.abs(b)*omegac*int_v
    K1 = special.kn(1,arg)
    K0 = special.kn(0,arg)    
    
    
    px = alffa_eff_x*1j*omegac*int_v*(K0 - INT_x)
    
    py = alffa_eff_y*1j*(2*1j*omegac*int_v*K1 - k0*INT_y)
    
    pz = alffa_eff_z*(-omegac*int_v*K1 + k0*INT_z )
    
    return px,py,pz












def dipole_moment_pole_aprox_resonance_v1_for_decay_rate(omegac,epsi_silica,d_nano_film,d_thickness_disk_nano,D_disk_nano,int_v,b,zp):     
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
    


    rtaself_x1, rtaself_y1, rtaself_z1  =  green_self_ana_v2(omegac,epsi_silica,d_nano_film,zp)
    rtaself_x2, rtaself_y2, rtaself_z2  =  green_self_num_integral_inside_light_cone(omegac,epsi_silica,d_nano_film,zp)

    rtaself_x, rtaself_y, rtaself_z  =  rtaself_x1 - rtaself_x2, rtaself_y1 - rtaself_y2, rtaself_z1 - rtaself_z2

    alfa_parallel = polarizability_parallel(E,d_thickness_disk_nano,D_disk_nano,epsi_silica)
    alfa_perp = polarizability_perp(E,d_thickness_disk_nano,D_disk_nano,epsi_silica)

    
    alffa_eff_x = (1/alfa_parallel -  rtaself_x)**(-1)
    alffa_eff_y = alffa_eff_x
    alffa_eff_z = (1/alfa_perp -  rtaself_z)**(-1)


    epsi_x = epsilon_x(E)
    epsi_HBN_par = epsi_x
    
    d_micro = d_nano_film*1e-3
    alfa_p = epsi_silica(E)*2/(omegac*d_micro*(epsi_HBN_par-1))

    cota_inf = 0.01/k0
    cota_sup = 600/k0   
    
    alpha_x = int_v    
    
#    term6 = np.sign(z)*pz*Rp*kp_2*J0*exp_electron
    alpha_parallel = lambda u: np.sqrt(alpha_x**2 + u**2)

    rp = lambda u: -alpha_parallel(u)/(alpha_parallel(u) - alfa_p)
      
    expo = lambda u: np.exp(-np.sqrt(alpha_x**2 + u**2)*k0*(2*zp + np.abs(b)))
    
    int_f_re_x = lambda u: np.real(rp(u)*expo(u)/alpha_parallel(u))
    int_f_im_x = lambda u: np.imag(rp(u)*expo(u)/alpha_parallel(u))
    
    INT_re_x,err = integrate.quad(int_f_re_x, cota_inf, cota_sup) 
    INT_im_x,err = integrate.quad(int_f_im_x, cota_inf, cota_sup) 
    
    INT_x = INT_re_x + 1j*INT_im_x
    
    
    int_f_re_y = lambda u: np.real(rp(u)*expo(u)*u/alpha_parallel(u))
    int_f_im_y = lambda u: np.imag(rp(u)*expo(u)*u/alpha_parallel(u))
    
    INT_re_y,err = integrate.quad(int_f_re_y, cota_inf, cota_sup) 
    INT_im_y,err = integrate.quad(int_f_im_y, cota_inf, cota_sup) 
    
    INT_y = INT_re_y + 1j*INT_im_y
        


    int_f_re_z = lambda u: np.real(rp(u)*expo(u))
    int_f_im_z = lambda u: np.imag(rp(u)*expo(u))
    
    INT_re_z,err = integrate.quad(int_f_re_z, cota_inf, cota_sup) 
    INT_im_z,err = integrate.quad(int_f_im_z, cota_inf, cota_sup) 
    
    INT_z = INT_re_z + 1j*INT_im_z


    arg = np.abs(b)*omegac*int_v
    K1 = special.kn(1,arg)
    K0 = special.kn(0,arg)    
    
    
    px = alffa_eff_x*1j*omegac*int_v*(K0 - INT_x)
    
    py = alffa_eff_y*1j*(2*1j*omegac*int_v*K1 - k0*INT_y)
    
    pz = alffa_eff_z*(-omegac*int_v*K1 + k0*INT_z )
    
    return px,py,pz








#     rp = lambda u: alfa_p/(alpha_parallel(u) - alfa_p)

def dipole_moment_pole_aprox_resonance_v2(omegac,epsi_silica,d_nano_film,d_thickness_disk_nano,D_disk_nano,int_v,b,zp):     
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
    


    rtaself_x, rtaself_y, rtaself_z  =  green_self_pole_aprox_v2(omegac,epsi_silica,d_nano_film,zp)

    alfa_parallel = polarizability_parallel(E,d_thickness_disk_nano,D_disk_nano,epsi_silica)
    alfa_perp = polarizability_perp(E,d_thickness_disk_nano,D_disk_nano,epsi_silica)

    
    alffa_eff_x = (1/alfa_parallel -  rtaself_x)**(-1)
    alffa_eff_y = alffa_eff_x
    alffa_eff_z = (1/alfa_perp -  rtaself_z)**(-1)


    epsi_x = epsilon_x(E)
    epsi_HBN_par = epsi_x
    
    d_micro = d_nano_film*1e-3
    alfa_p = epsi_silica(E)*2/(omegac*d_micro*(epsi_HBN_par-1))

    cota_inf = 0.01/k0
    cota_sup = 600/k0   
    
    alpha_x = int_v    
    
#    term6 = np.sign(z)*pz*Rp*kp_2*J0*exp_electron
    alpha_parallel = lambda u: np.sqrt(alpha_x**2 + u**2)

    rp = lambda u: -alfa_p/(alpha_parallel(u) - alfa_p)
      
    expo = lambda u: np.exp(-np.sqrt(alpha_x**2 + u**2)*k0*(2*zp + np.abs(b)))
    
    int_f_re_x = lambda u: np.real(rp(u)*expo(u)/alpha_parallel(u))
    int_f_im_x = lambda u: np.imag(rp(u)*expo(u)/alpha_parallel(u))
    
    INT_re_x,err = integrate.quad(int_f_re_x, cota_inf, cota_sup) 
    INT_im_x,err = integrate.quad(int_f_im_x, cota_inf, cota_sup) 
    
    INT_x = INT_re_x + 1j*INT_im_x
    
    
    int_f_re_y = lambda u: np.real(rp(u)*expo(u)*u/alpha_parallel(u))
    int_f_im_y = lambda u: np.imag(rp(u)*expo(u)*u/alpha_parallel(u))
    
    INT_re_y,err = integrate.quad(int_f_re_y, cota_inf, cota_sup) 
    INT_im_y,err = integrate.quad(int_f_im_y, cota_inf, cota_sup) 
    
    INT_y = INT_re_y + 1j*INT_im_y
        


    int_f_re_z = lambda u: np.real(rp(u)*expo(u))
    int_f_im_z = lambda u: np.imag(rp(u)*expo(u))
    
    INT_re_z,err = integrate.quad(int_f_re_z, cota_inf, cota_sup) 
    INT_im_z,err = integrate.quad(int_f_im_z, cota_inf, cota_sup) 
    
    INT_z = INT_re_z + 1j*INT_im_z


    arg = np.abs(b)*omegac*int_v
    K1 = special.kn(1,arg)
    K0 = special.kn(0,arg)    
    
    
    px = alffa_eff_x*1j*omegac*int_v*(K0 - INT_x)
    
    py = alffa_eff_y*1j*(2*1j*omegac*int_v*K1 - k0*INT_y)
    
    pz = alffa_eff_z*(-omegac*int_v*K1 + k0*INT_z )
    
    return px,py,pz


#%%
    