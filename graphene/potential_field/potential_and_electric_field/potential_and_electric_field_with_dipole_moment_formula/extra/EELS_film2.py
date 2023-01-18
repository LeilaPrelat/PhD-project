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
    sys.path.insert(1, path_constants)
    from graphene_sigma import sigma_DL
except ModuleNotFoundError:
    print('graphene_sigma.py no se encuentra en ' + path_basic)


try:
    sys.path.insert(1, path_constants)
    from dipole_moment import dipole_moment_anav2
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

def EELS_film_ana_f(omegac,epsi1,epsi2,hbmu,hbgama,int_v,px,py,pz,b,zp,phi,R,z,omega0,kappa_factor_omega0,kappa_r_factor):     
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
    px_v,py_v,pz_v = dipole_moment_anav2(omegac,epsi1,epsi2,hbmu,hbgama,int_v,px,py,pz,b,zp,omega0,kappa_factor_omega0,kappa_r_factor)
#    print(px_v,py_v,pz_v)
    
    
    
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    Rp = 2*epsi1/(epsi1 + epsi2)
    alfa_p = 1j*(epsi1 + epsi2)/(cond*cte1)
    
    kp = alfa_p*k1
    
#    kp = np.real(kp)   ### neglecting the imaginary part of kp para la formula de EELS (la del campo)
#    alfa_p = np.real(alfa_p)  ### neglecting the imaginary part of kp para la formula de EELS (la del campo)
    
    kp_3 = kp**3 
    kp_2 = kp**2
    
    z_dip_barra = k1*(2*zp - z)        
    exp_electron = np.exp(-alfa_p*z_dip_barra)
    
#    phi = np.arctan2(np.abs(y),np.abs(x))    
#    R = np.sqrt(x**2 + y**2)
    term_px_py = px_v*np.cos(phi) + py_v*np.sin(phi)    
#    J1 = special.jv(1,kp*R)
    H0 = special.hankel1(0,kp*R)
    H1 = special.hankel1(1,kp*R)
    H_minus1 = special.hankel1(-1,kp*R)
    H2 = special.hankel1(2,kp*R)
    


    Ex1 = Rp*np.pi*1j*kp_3*0.5*pz_v*np.cos(phi)*np.sign(z)*(H_minus1 - H1)
    Ex2 = Rp*np.pi*1j*kp_3*0.5*np.cos(phi)*term_px_py*(H0 - H2)  
    Ex3 = -Rp*np.pi*1j*kp_2*H1*(-px_v*np.sin(phi)**2 + py_v*np.sin(phi)/np.cos(phi) - py_v*np.sin(phi)**3/np.cos(phi))/R  
    
    Ex = (Ex1 + Ex2 + Ex3)*exp_electron
    
    Ez = z*Rp*np.pi*1j*kp_3*(-pz_v*np.sign(z)*H0 + term_px_py*H1)*exp_electron
#    print(pz_v*np.sign(z)*(H_minus1 - H1))
#    
#    print( term_px_py*(H0 - H2))
#    
    Ey1 = Rp*np.pi*1j*kp_3*0.5*pz_v*np.sin(phi)*np.sign(z)*(H_minus1 - H1)
    Ey2 = Rp*np.pi*1j*kp_3*0.5*np.sin(phi)*term_px_py*(H0 - H2)  
    Ey3 = -Rp*np.pi*1j*kp_2*H1*(px_v*np.sin(phi)*np.cos(phi) - py_v + py_v*np.sin(phi)**2)/R  
    
    Ey = (Ey1 + Ey2 + Ey3)*exp_electron
    
    
    term_f2 = Ex*np.conjugate(px_v) + Ey*np.conjugate(py_v) + Ez*np.conjugate(pz_v)
###################################################################################################

    return term_f2 


#%%