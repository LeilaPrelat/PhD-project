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

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/potential_field/potential_and_electric_field/potential_and_electric_field_with_dipole_moment_formula','')
#print('Importar modulos necesarios para este codigo')


try:
    sys.path.insert(1, path_basic)
    from dipole_moment import dipole_moment_anav2_resonance,dipole_moment_anav2_for_decay_rate_resonance,dipole_moment_anav2_for_decay_rate_resonance_dir
except ModuleNotFoundError:
    print('dipole_moment.py no se encuentra en ' + path_basic)

try:
    sys.path.insert(1, path_basic)
    from green_self_image import green_self_num,green_self_ana2,green_self_num_integral_inside_light_cone ## la parte imaginaria analitica del Green self image funciona bastante bien
    
except ModuleNotFoundError:
    print('graphene_sigma.py no se encuentra en ' + path_basic)

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

def EELS_film_ana_f(omegac,epsi1,epsi2,hbmu,hbgama,int_v,b,zp):     
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
    -------9
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
#    k1_2 = k1**2
    px_v,py_v,pz_v = dipole_moment_anav2_for_decay_rate_resonance(omegac,epsi1,epsi2,hbmu,hbgama,int_v,b,zp)
#    print(px_v,py_v,pz_v)
#    px_tot_2 = np.abs(px_v)**2 + np.abs(py_v)**2 + np.abs(pz_v)**2 
    rtaself_x1, rtaself_y1, rtaself_z1  =  green_self_ana2(omegac,epsi1,epsi2,hbmu,hbgama,zp)
    rtaself_x2, rtaself_y2, rtaself_z2  =  green_self_num_integral_inside_light_cone(omegac,epsi1,epsi2,hbmu,hbgama,zp)
    
    rtaself_x, rtaself_y, rtaself_z  = rtaself_x1 - rtaself_x2, rtaself_y1 - rtaself_y2, rtaself_z1 - rtaself_z2
    
    
    Green_self = rtaself_x*(np.abs(px_v)**2) + rtaself_y*(np.abs(py_v)**2)  + rtaself_z*(np.abs(pz_v)**2)
#    print(rtaself_x)

#    kp = np.real(kp)   ### neglecting the imaginary part of kp para la formula de EELS (la del campo)
#    alfa_p = np.real(alfa_p)  ### neglecting the imaginary part of kp para la formula de EELS (la del campo)
 
###################################################################################################

    cte_aux = alfac*(int_v**2)/(2*(np.pi**2)*c)
    
    cte_aux = cte_aux*1e9 ### cambiar unidades
    
    return cte_aux*np.imag(Green_self) 

def EELS_film_ana_f_sin_integrar(omegac,epsi1,epsi2,hbmu,hbgama,int_v,b,zp):     
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
#    k1_2 = k1**2
    px_v,py_v,pz_v = dipole_moment_anav2_resonance(omegac,epsi1,epsi2,hbmu,hbgama,int_v,b,zp)
#    print(px_v,py_v,pz_v)
#    px_tot_2 = np.abs(px_v)**2 + np.abs(py_v)**2 + np.abs(pz_v)**2 
    rtaself_x, rtaself_y, rtaself_z  =  green_self_num_integral_inside_light_cone(omegac,epsi1,epsi2,hbmu,hbgama,zp)
    
    Green_self = rtaself_x*(np.abs(px_v)**2) + rtaself_y*(np.abs(py_v)**2)  + rtaself_z*(np.abs(pz_v)**2)
#    print(rtaself_x)

#    kp = np.real(kp)   ### neglecting the imaginary part of kp para la formula de EELS (la del campo)
#    alfa_p = np.real(alfa_p)  ### neglecting the imaginary part of kp para la formula de EELS (la del campo)
 
###################################################################################################

    cte_aux = alfac*(int_v**2)/(2*(np.pi**2)*c)
    
    cte_aux = cte_aux*1e9 ### cambiar unidades
    
    return cte_aux*np.imag(Green_self) 


def EELS_film_ana_f_div_gamma0(omegac,epsi1,epsi2,hbmu,hbgama,int_v,b,zp):     
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
#    k1_2 = k1**2
    px_v,py_v,pz_v = dipole_moment_anav2_for_decay_rate_resonance(omegac,epsi1,epsi2,hbmu,hbgama,int_v,b,zp)
#    print(px_v,py_v,pz_v)
#    px_tot_2 = np.abs(px_v)**2 + np.abs(py_v)**2 + np.abs(pz_v)**2 
    rtaself_x1, rtaself_y1, rtaself_z1  =  green_self_ana2(omegac,epsi1,epsi2,hbmu,hbgama,zp)
    rtaself_x2, rtaself_y2, rtaself_z2  =  green_self_num_integral_inside_light_cone(omegac,epsi1,epsi2,hbmu,hbgama,zp)
    
    rtaself_x, rtaself_y, rtaself_z  = rtaself_x1 - rtaself_x2, rtaself_y1 - rtaself_y2, rtaself_z1 - rtaself_z2
    
    Green_self = rtaself_x*(np.abs(px_v)**2) + rtaself_y*(np.abs(py_v)**2)  + rtaself_z*(np.abs(pz_v)**2)
#    print(rtaself_x)

#    kp = np.real(kp)   ### neglecting the imaginary part of kp para la formula de EELS (la del campo)
#    alfa_p = np.real(alfa_p)  ### neglecting the imaginary part of kp para la formula de EELS (la del campo)
 
###################################################################################################
    v = c/int_v
    cte_aux = 1/(2*np.pi*v)
    
#    cte_aux = cte_aux*1e9 ### cambiar unidades

    gamma = np.sqrt(1 - (int_v)**(-2))**(-1)
#    gamma = 1
    alpha = -3*epsi1/(4*1j*k1**3)

    arg = np.abs(b)*omegac*int_v/gamma
    K1 = special.kn(1,arg)
    K0 = special.kn(0,arg)
   

    factor_gamma0 = (2*omegac*int_v/(v*gamma))**2
    gamma0 = factor_gamma0*(K0**2/gamma**2 + K1**2)*np.imag(alpha)/np.pi  ## decay rate de 1 dipolo # pero sin el "e/hbar" se cancela con el momento dipolar^2
    
#    gamma0 = 1
    rta_v1 = cte_aux*np.imag(Green_self)/gamma0
    
    alpha = 3*epsi1/(2*k1**3)
    return rta_v1
#    
#    rta_v2 = np.imag(Green_self)*(v*gamma)**2/(8*np.pi*(K0**2 + K1**2)*alpha)
#
#    return rta_v2



def EELS_film_ana_f_div_gamma0_v2(omegac,epsi1,epsi2,hbmu,hbgama,int_v,b,zp):     
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
#    k1_2 = k1**2
    px_v,py_v,pz_v = dipole_moment_anav2_for_decay_rate_resonance(omegac,epsi1,epsi2,hbmu,hbgama,int_v,b,zp) # multiplicar por e/(2*pi*v)
    

#    print(px_v,py_v,pz_v)
#    px_tot_2 = np.abs(px_v)**2 + np.abs(py_v)**2 + np.abs(pz_v)**2 
    rtaself_x1, rtaself_y1, rtaself_z1  =  green_self_ana2(omegac,epsi1,epsi2,hbmu,hbgama,zp)
    rtaself_x2, rtaself_y2, rtaself_z2  =  green_self_num_integral_inside_light_cone(omegac,epsi1,epsi2,hbmu,hbgama,zp)
    
    rtaself_x, rtaself_y, rtaself_z  = rtaself_x1 - rtaself_x2, rtaself_y1 - rtaself_y2, rtaself_z1 - rtaself_z2
    
    Green_self = rtaself_x*(np.abs(px_v)**2) + rtaself_y*(np.abs(py_v)**2)  + rtaself_z*(np.abs(pz_v)**2)

    usar_dif_p = 1
    usar_mismo_p = 0

    if usar_dif_p == 1:  
        
        px_dir,py_dir,pz_dir = dipole_moment_anav2_for_decay_rate_resonance_dir(omegac,int_v,b,zp)
        
        denominador = np.abs(px_dir)**2 +  np.abs(py_dir)**2 +  np.abs(pz_dir)**2

        rta = np.imag(Green_self*2*(omegac**3)/denominador)    
    
    if usar_mismo_p == 1: ## aca los momentos p se cancelan  
        
        px_dir,py_dir,pz_dir = dipole_moment_anav2_for_decay_rate_resonance(omegac,epsi1,epsi2,hbmu,hbgama,int_v,b,zp) # multiplicar por e/(2*pi*v)

        denominador = np.abs(px_dir)**2 +  np.abs(py_dir)**2 +  np.abs(pz_dir)**2

        alffa_eff_x = 1j*(2*omegac**3/(3*epsi1) +  np.imag(rtaself_x))**(-1)
        alffa_eff_y = 1j*(2*omegac**3/(3*epsi1) +  np.imag(rtaself_y))**(-1)
        alffa_eff_z = 1j*(2*omegac**3/(3*epsi1) +  np.imag(rtaself_z))**(-1)
    
        alfa_eff = np.abs(alffa_eff_x)**2 + np.abs(alffa_eff_y)**2 + np.abs(alffa_eff_z)**2 

    
        rta = np.imag(Green_self*alfa_eff*2*(omegac**3)/denominador )

    return rta 



def EELS_film_ana_f_div_gamma0_simpler(omegac,epsi1,epsi2,hbmu,hbgama,zp):     
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

    k = omegac
    k_3 = k**3

    rtaself_x1, rtaself_y1, rtaself_z1  =  green_self_ana2(omegac,epsi1,epsi2,hbmu,hbgama,zp)
    rtaself_x2, rtaself_y2, rtaself_z2  =  green_self_num_integral_inside_light_cone(omegac,epsi1,epsi2,hbmu,hbgama,zp)
    
    rtaself_x, rtaself_y, rtaself_z  = rtaself_x1 - rtaself_x2, rtaself_y1 - rtaself_y2, rtaself_z1 - rtaself_z2
    
    
    
    rta = 3*(rtaself_x + rtaself_y + rtaself_z)/(2*k_3*np.sqrt(epsi1))
   
    
    return np.imag(rta)


#%%

def EELS_dir_ana_f(omegac,epsi1,epsi2,hbmu,hbgama,int_v,b,zp):     
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
    k_3 = k1**3

#    k1_2 = (k0*cte1)**2
 #   n_v1 = int_v/cte1
#    k1_2 = k1**2

    px_v,py_v,pz_v = dipole_moment_anav2_resonance(omegac,epsi1,epsi2,hbmu,hbgama,int_v,b,zp)

    tot = np.abs(px_v)**2 + np.abs(py_v)**2  + np.abs(pz_v)**2
 
###################################################################################################
    
    
    cte_aux = k_3*alfac*(int_v**2)/(3*(np.pi**2)*c)

    cte_aux = cte_aux*1e9 ## cambiar unidades
    
    return cte_aux*tot 


#%%
    

def EELS_parallel_ana_f(omegac,epsi1,epsi2,hbmu,hbgama,L_nano,int_v,b,zp):     
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
#    k1_2 = k1**2

#    print(rtaself_x)

#    kp = np.real(kp)   ### neglecting the imaginary part of kp para la formula de EELS (la del campo)
#    alfa_p = np.real(alfa_p)  ### neglecting the imaginary part of kp para la formula de EELS (la del campo)
 
###################################################################################################

    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    Rp = 2*epsi1/(epsi1 + epsi2)
    alfa_p = 1j*(epsi1 + epsi2)/(cond*cte1)
    kp = alfa_p*k1
    
    expo = np.exp(-2*kp*(zp-b))
    
    rta = Rp*kp*expo/(np.sqrt(kp**2 - (int_v*omegac)**2))
    
    L_micro = L_nano*1e-3
    
    cte_aux = 1e9 ### cambiar unidades
    
    return np.real(rta)*4*alfac*(int_v**2)*L_micro*cte_aux/c





