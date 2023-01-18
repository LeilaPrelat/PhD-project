#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila
"""
import numpy as np
import sys
import os 

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/potential_field/potential_and_electric_field/potential_and_electric_field_with_dipole_moment_formula','')
#print('Importar modulos necesarios para este codigo')


try:
    sys.path.insert(1, path_basic)
    from dipole_moment import dipole_moment_anav2_resonance, dipole_moment_pole_aprox_resonance, dipole_moment_num_resonance
except ModuleNotFoundError:
    print('dipole_moment.py no se encuentra en ' + path_basic)


try:
    sys.path.insert(1, path_basic)
    from green_self_image import green_self_ana2
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

def EELS_film_ana_f(omegac,epsi1,epsi3,d_nano,int_v,b,zp):     
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
    px_v,py_v,pz_v = dipole_moment_anav2_resonance(omegac,epsi1,epsi3,d_nano,int_v,b,zp)
#    print(px_v,py_v,pz_v)
#    px_tot_2 = np.abs(px_v)**2 + np.abs(py_v)**2 + np.abs(pz_v)**2 
    rtaself_x, rtaself_y, rtaself_z  =  green_self_ana2(omegac,epsi1,epsi3,d_nano,zp)
    
    Green_self_parallel  = rtaself_x*(np.abs(px_v)**2) + rtaself_y*(np.abs(py_v)**2)  
    Green_self_perp  =  rtaself_z*(np.abs(pz_v)**2)

#    print(rtaself_x)

#    kp = np.real(kp)   ### neglecting the imaginary part of kp para la formula de EELS (la del campo)
#    alfa_p = np.real(alfa_p)  ### neglecting the imaginary part of kp para la formula de EELS (la del campo)
 
###################################################################################################

    cte_aux = alfac*(int_v**2)/(2*(np.pi**2)*(c))
    
    cte_aux = cte_aux*1e6 ### cambiar unidades
    
    return cte_aux*np.imag(Green_self_parallel) ,  cte_aux*np.imag(Green_self_perp) 


#%%


def EELS_dir_ana_f(omegac,epsi1,epsi3,d_nano,int_v,b,zp):     
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

    px_v,py_v,pz_v = dipole_moment_anav2_resonance(omegac,epsi1,epsi3,d_nano,int_v,b,zp)

    tot_parallel = np.abs(px_v)**2 + np.abs(py_v)**2  
    tot_perp = np.abs(pz_v)**2
    
    
    cte_aux = k_3*alfac*(int_v**2)/(3*(np.pi**2)*c)

    
    cte_aux = cte_aux*1e6 ### cambiar unidades
    
    return cte_aux*tot_parallel,  cte_aux*tot_perp



