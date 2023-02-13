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
path_constants =  path_basic.replace('/potential_field/potential_and_electric_field_with_dipole_moment_formula','')
#print('Importar modulos necesarios para este codigo')


try:
    sys.path.insert(1, path_basic)
    from dipole_moment_dif_sign import dipole_moment_pole_aprox_resonance_v1,dipole_moment_pole_aprox_resonance_v1_for_decay_rate,dipole_moment_pole_aprox_resonance_v2, dipole_moment_num_resonance,dipole_moment_anav1_for_decay_rate_resonance,dipole_moment_anav2_for_decay_rate_resonance, dipole_moment_ana_for_decay_rate_resonance_dir
except ModuleNotFoundError:
    print('dipole_moment.py no se encuentra en ' + path_basic)


try:
    sys.path.insert(1, path_constants)
    from green_self_image_dif_sign import green_self_ana_v1,green_self_ana_v2,green_self_num_integral_inside_light_cone,green_self_pole_aprox_v1,green_self_num
except ModuleNotFoundError:
    print('green_self_image_dif_sign.py no se encuentra en ' + path_constants)

 
try:
    sys.path.insert(1, path_constants)
    from constants import constantes
except ModuleNotFoundError:
    print('constants.py no se encuentra en ' + path_constants)

pi,hb,c,alfac,mu1,mu2 = constantes()
aux = c*hb

#%%

def EELS_film_ana_f(omegac,epsi_silica,d_nano_film,d_thickness_disk_nano,D_disk_nano,int_v,b,zp):     
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
#
#    E = omegac*aux
##    k0 = omegac #=omega/c
#    # x_y = ky/k0
#    n1 = epsi1*mu1
#    cte1 = np.sqrt(n1)
#    k1 = omegac*cte1
#    k1_2 = (k0*cte1)**2
 #   n_v1 = int_v/cte1
#    k1_2 = k1**2
    px_v,py_v,pz_v = dipole_moment_anav1_for_decay_rate_resonance(omegac,epsi_silica,d_nano_film,d_thickness_disk_nano,D_disk_nano,int_v,b,zp)


#    print(px_v,py_v,pz_v)
#    px_tot_2 = np.abs(px_v)**2 + np.abs(py_v)**2 + np.abs(pz_v)**2 
    rtaself_x1, rtaself_y1, rtaself_z1  =  green_self_ana_v1(omegac,epsi_silica,d_nano_film,zp)
    rtaself_x2, rtaself_y2, rtaself_z2  =  green_self_num_integral_inside_light_cone(omegac,epsi_silica,d_nano_film,zp)

    rtaself_x, rtaself_y, rtaself_z  =  rtaself_x1 - rtaself_x2, rtaself_y1 - rtaself_y2, rtaself_z1 - rtaself_z2
#    rtaself_x, rtaself_y, rtaself_z  =  rtaself_x1 , rtaself_y1 , rtaself_z1 
    
#    print(rtaself_x)
    Green_self = rtaself_x*(np.abs(px_v)**2) + rtaself_y*(np.abs(py_v)**2)  + rtaself_z*(np.abs(pz_v)**2)
#    print(rtaself_x)

#    kp = np.real(kp)   ### neglecting the imaginary part of kp para la formula de EELS (la del campo)
#    alfa_p = np.real(alfa_p)  ### neglecting the imaginary part of kp para la formula de EELS (la del campo)
 
###################################################################################################

    cte_aux = alfac*(int_v**2)/(2*(np.pi**2)*(c))
    
    cte_aux = cte_aux*1e9 ### cambiar unidades
    
    return cte_aux*np.imag(Green_self) 


def EELS_film_ana_f_div_gamma0(omegac,epsi_silica,d_nano_film,d_thickness_disk_nano,D_disk_nano,int_v,b,zp):     
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
#    n1 = epsi1*mu1
#    cte1 = np.sqrt(n1)
    k1 = omegac
#    k1_2 = (k0*cte1)**2
 #   n_v1 = int_v/cte1
#    k1_2 = k1**2
#    px_v,py_v,pz_v = dipole_moment_anav1_for_decay_rate_resonance(omegac,epsi_silica,d_nano,int_v,b,zp)
#    print(px_v,py_v,pz_v)
#    px_tot_2 = np.abs(px_v)**2 + np.abs(py_v)**2 + np.abs(pz_v)**2 
    
    px_v,py_v,pz_v = dipole_moment_pole_aprox_resonance_v1(omegac,epsi_silica,d_nano_film,d_thickness_disk_nano,D_disk_nano,int_v,b,zp)
#    print(px_v,py_v,pz_v)
#    px_tot_2 = np.abs(px_v)**2 + np.abs(py_v)**2 + np.abs(pz_v)**2 
    rtaself_x1, rtaself_y1, rtaself_z1  =  green_self_pole_aprox_v1(omegac,epsi_silica,d_nano_film,zp)

#    rtaself_x1, rtaself_y1, rtaself_z1  =  green_self_ana_v1(omegac,epsi_silica,d_nano,zp)
    rtaself_x2, rtaself_y2, rtaself_z2  =  green_self_num_integral_inside_light_cone(omegac,epsi_silica,d_nano_film,zp)

    rtaself_x, rtaself_y, rtaself_z  =  rtaself_x1 - rtaself_x2, rtaself_y1 - rtaself_y2, rtaself_z1 - rtaself_z2
    
    Green_self = rtaself_x*(np.abs(px_v)**2) + rtaself_y*(np.abs(py_v)**2)  + rtaself_z*(np.abs(pz_v)**2)
#    print(rtaself_x)

#    kp = np.real(kp)   ### neglecting the imaginary part of kp para la formula de EELS (la del campo)
#    alfa_p = np.real(alfa_p)  ### neglecting the imaginary part of kp para la formula de EELS (la del campo)
 
###################################################################################################
    v = c/int_v
    cte_aux = 1/(2*np.pi*v)
    
#    cte_aux = cte_aux*1e9 ### cambiar unidades

    gamma = np.sqrt(1 - (int_v)**(-2))**(-1)
    alpha = -3*epsi_silica(E)/(4*1j*k1**3)

    arg = np.abs(b)*omegac*int_v/gamma
    K1 = special.kn(1,arg)
    K0 = special.kn(0,arg)
    
#    print(K1)
   

    factor_gamma0 = (2*omegac*int_v/(v*gamma))**2
    gamma0 = factor_gamma0*(K0**2/gamma**2 + K1**2)*np.imag(alpha)/np.pi  ## decay rate de 1 dipolo # pero sin el "e/hbar" se cancela con el momento dipolar^2
    
    return cte_aux*np.imag(Green_self)/gamma0




def EELS_film_ana_f_div_gamma0_v2(omegac,epsi_silica,d_nano_film,d_thickness_disk_nano,D_disk_nano,int_v,b,zp):     
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

#    k1_2 = (k0*cte1)**2
 #   n_v1 = int_v/cte1
#    k1_2 = k1**2
    px_v,py_v,pz_v = dipole_moment_pole_aprox_resonance_v1(omegac,epsi_silica,d_nano_film,d_thickness_disk_nano,D_disk_nano,int_v,b,zp)
#    print(px_v,py_v,pz_v)
#    px_tot_2 = np.abs(px_v)**2 + np.abs(py_v)**2 + np.abs(pz_v)**2 
    rtaself_x1, rtaself_y1, rtaself_z1  =  green_self_pole_aprox_v1(omegac,epsi_silica,d_nano_film,zp)
    rtaself_x2, rtaself_y2, rtaself_z2  =  green_self_num_integral_inside_light_cone(omegac,epsi_silica,d_nano_film,zp)

    rtaself_x, rtaself_y, rtaself_z  =  rtaself_x1 - rtaself_x2, rtaself_y1 - rtaself_y2, rtaself_z1 - rtaself_z2
    
    Green_self = rtaself_x*(np.abs(px_v)**2) + rtaself_y*(np.abs(py_v)**2)  + rtaself_z*(np.abs(pz_v)**2)
#    print(rtaself_x)

#    kp = np.real(kp)   ### neglecting the imaginary part of kp para la formula de EELS (la del campo)
#    alfa_p = np.real(alfa_p)  ### neglecting the imaginary part of kp para la formula de EELS (la del campo)
 
###################################################################################################
    usar_dif_p = 0
    usar_mismo_p = 1

    if usar_dif_p == 1:  
        
        px_dir,py_dir,pz_dir = dipole_moment_ana_for_decay_rate_resonance_dir(omegac,int_v,b,zp)
        
        denominador = np.abs(px_dir)**2 +  np.abs(py_dir)**2 +  np.abs(pz_dir)**2

        rta = np.imag(Green_self*2*(omegac**3)/denominador)    
    
    if usar_mismo_p == 1: ## aca los momentos p se cancelan  
        
        px_dir,py_dir,pz_dir = dipole_moment_pole_aprox_resonance_v1(omegac,epsi_silica,d_nano_film,d_thickness_disk_nano,D_disk_nano,int_v,b,zp) # multiplicar por e/(2*pi*v)

        denominador = np.abs(px_dir)**2 +  np.abs(py_dir)**2 +  np.abs(pz_dir)**2

        alffa_eff_x = 1j*(2*omegac**3/(3*epsi_silica(E)) +  np.imag(rtaself_x))**(-1)
        alffa_eff_y = 1j*(2*omegac**3/(3*epsi_silica(E)) +  np.imag(rtaself_y))**(-1)
        alffa_eff_z = 1j*(2*omegac**3/(3*epsi_silica(E)) +  np.imag(rtaself_z))**(-1)
    
        alfa_eff = np.abs(alffa_eff_x)**2 + np.abs(alffa_eff_y)**2 + np.abs(alffa_eff_z)**2 

    
        rta = np.imag(Green_self*alfa_eff*2*(omegac**3)/denominador )

    return rta



def EELS_film_ana_f_div_gamma0_v3(omegac,epsi_silica,d_nano_film,d_thickness_disk_nano,D_disk_nano,int_v,b,zp):     ## normalizando con el paper 149
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
#    k1_2 = (k0*cte1)**2
 #   n_v1 = int_v/cte1
#    k1_2 = k1**2
    px_v,py_v,pz_v = dipole_moment_pole_aprox_resonance_v1_for_decay_rate(omegac,epsi_silica,d_nano_film,d_thickness_disk_nano,D_disk_nano,int_v,b,zp) # multiplicar por e/(2*pi*v)
    

#    print(px_v,py_v,pz_v)
#    px_tot_2 = np.abs(px_v)**2 + np.abs(py_v)**2 + np.abs(pz_v)**2 
    rtaself_x1, rtaself_y1, rtaself_z1  =  green_self_pole_aprox_v1(omegac,epsi_silica,d_nano_film,zp)
    rtaself_x2, rtaself_y2, rtaself_z2  =  green_self_num_integral_inside_light_cone(omegac,epsi_silica,d_nano_film,zp)

    rtaself_x, rtaself_y, rtaself_z  =  rtaself_x1 - rtaself_x2, rtaself_y1 - rtaself_y2, rtaself_z1 - rtaself_z2
    
    Green_self = rtaself_x*(np.abs(px_v)**2) + rtaself_y*(np.abs(py_v)**2)  + rtaself_z*(np.abs(pz_v)**2)

    arg = np.abs(b)*omegac*int_v
    K1 = special.kn(1,arg)
    K0 = special.kn(0,arg)
   

    
    factor_K = K0**2 + K1**2  ## decay rate de 1 dipolo # pero sin el "e/hbar" se cancela con el momento dipolar^2

    
#    px_dir,py_dir,pz_dir = dipole_moment_anav2_for_decay_rate_resonance_dir(omegac,int_v,b,zp)        
#    denominador = np.abs(px_dir)**2 +  np.abs(py_dir)**2 +  np.abs(pz_dir)**2

    k_prima = omegac*np.sqrt(epsi_silica(E))
#    epsi1 = 1
#    k_prima = omegac*np.sqrt(epsi1)
#    print(k_prima_2/k_prima)
    
    factor_final = k_prima*2*np.pi/(12*np.pi*(int_v**2))

    rta = np.imag(Green_self)*factor_final/factor_K    
    

    return rta 



 



#%%
    
def EELS_film_num_f(omegac,epsi_silica,d_nano_film,d_thickness_disk_nano,D_disk_nano,int_v,b,zp):     
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

#    E = omegac*aux
##    k0 = omegac #=omega/c
#    # x_y = ky/k0
#    n1 = epsi1*mu1
#    cte1 = np.sqrt(n1)
#    k1 = omegac*cte1
#    k1_2 = (k0*cte1)**2
 #   n_v1 = int_v/cte1
#    k1_2 = k1**2
    px_v,py_v,pz_v = dipole_moment_num_resonance(omegac,epsi_silica,d_nano_film,d_thickness_disk_nano,D_disk_nano,int_v,b,zp)
#    print(px_v,py_v,pz_v)
#    px_tot_2 = np.abs(px_v)**2 + np.abs(py_v)**2 + np.abs(pz_v)**2 
    rtaself_x, rtaself_y, rtaself_z  =  green_self_num(omegac,epsi_silica,d_nano_film,zp)
    
    Green_self = rtaself_x*(np.abs(px_v)**2) + rtaself_y*(np.abs(py_v)**2)  + rtaself_z*(np.abs(pz_v)**2)
#    print(rtaself_x)

#    kp = np.real(kp)   ### neglecting the imaginary part of kp para la formula de EELS (la del campo)
#    alfa_p = np.real(alfa_p)  ### neglecting the imaginary part of kp para la formula de EELS (la del campo)
 
###################################################################################################

    cte_aux = alfac*(int_v**2)/(2*(np.pi**2)*(c))
    
    cte_aux = cte_aux*1e9 ### cambiar unidades
    
    return cte_aux*np.imag(Green_self) 


#%%


def EELS_film_pole_aprox_f(omegac,epsi_silica,d_nano_film,d_thickness_disk_nano,D_disk_nano,int_v,b,zp):     
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

#    E = omegac*aux
##    k0 = omegac #=omega/c
#    # x_y = ky/k0
#    n1 = epsi1*mu1
#    cte1 = np.sqrt(n1)
#    k1 = omegac*cte1
#    k1_2 = (k0*cte1)**2
 #   n_v1 = int_v/cte1
#    k1_2 = k1**2
    px_v,py_v,pz_v = dipole_moment_pole_aprox_resonance_v1_for_decay_rate(omegac,epsi_silica,d_nano_film,d_thickness_disk_nano,D_disk_nano,int_v,b,zp)
#    print(px_v,py_v,pz_v)
#    px_tot_2 = np.abs(px_v)**2 + np.abs(py_v)**2 + np.abs(pz_v)**2 
    rtaself_x1, rtaself_y1, rtaself_z1  =  green_self_pole_aprox_v1(omegac,epsi_silica,d_nano_film,zp)
    rtaself_x2, rtaself_y2, rtaself_z2  =  green_self_num_integral_inside_light_cone(omegac,epsi_silica,d_nano_film,zp)
    
    rtaself_x, rtaself_y, rtaself_z  =  rtaself_x1 - rtaself_x2, rtaself_y1 - rtaself_y2, rtaself_z1 - rtaself_z2
   
    Green_self = rtaself_x*(np.abs(px_v)**2) + rtaself_y*(np.abs(py_v)**2)  + rtaself_z*(np.abs(pz_v)**2)
#    print(rtaself_x)

#    kp = np.real(kp)   ### neglecting the imaginary part of kp para la formula de EELS (la del campo)
#    alfa_p = np.real(alfa_p)  ### neglecting the imaginary part of kp para la formula de EELS (la del campo)
 
###################################################################################################

    cte_aux = alfac*(int_v**2)/(2*(np.pi**2)*(c))
    
    cte_aux = cte_aux*1e9 ### cambiar unidades
    
    return cte_aux*np.imag(Green_self) 


#%%
    
