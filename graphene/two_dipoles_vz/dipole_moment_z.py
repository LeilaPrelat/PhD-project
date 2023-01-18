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
path_constants =  path_basic.replace('/two_dipoles_vz','')
#print('Importar modulos necesarios para este codigo')

try:
    sys.path.insert(1, path_basic)
    from green_self_image import green_self_num,green_self_ana2,green_self_pole_aprox
except ModuleNotFoundError:
    print('green_self_image.py no se encuentra en ' + path_basic)


try:
    sys.path.insert(1, path_basic)
    from green_tensor_direct_zz import green_tensor_ANA_QE
except ModuleNotFoundError:
    print('green_tensor_direct.py no se encuentra en ' + path_basic)


try:
    sys.path.insert(1, path_basic)
    from green_tensor_ref_zz_graphene import green_tensor_PP2, green_tensor_ref_pole_aprox
except ModuleNotFoundError:
    print('green_tensor_ref_xx_graphene.py no se encuentra en ' + path_basic)

try:
    sys.path.insert(1, path_basic)
    from fieldE_direct_z import Efield_ANA
except ModuleNotFoundError:
    print('fieldE_direct_z.py no se encuentra en ' + path_basic)


#try:
#    sys.path.insert(1, path_basic)
#    from fieldE_ref_z import Efield_ref_ana,Efield_ref_fresnel
#except ModuleNotFoundError:
#    print('fieldE_ref.py no se encuentra en ' + path_basic)     

try:
    sys.path.insert(1, path_constants)
    from constants import constantes
except ModuleNotFoundError:
    print('constants.py no se encuentra en ' + path_constants)

pi,hb,c,alfac,mu1,mu2 = constantes()
aux = c*hb


#%%

def alpha_function(epsi1,omegac,omega0,kappa_factor_omega0,kappa_r_factor):
    """
    Parameters
    ----------
    epsilon1 : permeabilidad electrica del medio 1
    omegac : frequencia in units of 1/micrometers 
    omega0 : resonant frequency 
    kappa_factor_omega0 : kappa (collision frequency) = kappa_factor*omega0
    kappa_r_factor : kappa_r = kappa_r_factor*kappa with kappa_r_factor<1 
    Returns
    -------
    lorentzian model for polarizabilty 
    """
    omega = omegac*c
    n1 = epsi1*mu1
    cte1 = np.sqrt(n1)
    k1 = cte1
    k1_3 = k1**3
    
    kappa = kappa_factor_omega0*omega0
    kappa_r = kappa_r_factor*kappa
    A = 3*kappa_r*0.25/k1_3
    
    den = (omega0 - omega)**2 + (kappa/2)**2
    num = omega0 - omega + 1j*kappa/2

    rta = A*num/den
    return rta

#%%

def dipole_moment_z_ana(omegac,epsi1,epsi2,hbmu,hbgama,int_v,zp,x1,z1,x2,z2,omega01,kappa_factor_omega01,kappa_r_factor1,omega02,kappa_factor_omega02,kappa_r_factor2):     
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

    alffa1 = alpha_function(epsi1,omegac,omega01,kappa_factor_omega01,kappa_r_factor1)
    alffa2 = alpha_function(epsi1,omegac,omega02,kappa_factor_omega02,kappa_r_factor2)
    
    rtaself_x, rtaself_y, rtaself_z  =  green_self_pole_aprox(omegac,epsi1,epsi2,hbmu,hbgama,zp)
    alffa_eff1_x = (1/alffa1 -  rtaself_x)**(-1)

    alffa_eff2_x = (1/alffa2 -  rtaself_x)**(-1)

    phi = 0
    
    R = np.abs(x1-x2)
    Gdir_12 = green_tensor_ANA_QE(epsi1,R,z1,z2)

    Gref_12 = green_tensor_PP2(omegac,epsi1,epsi2,hbmu,hbgama,R,phi,z1,z2,zp)

    R1 = np.abs(x1)
    Edir1 = Efield_ANA(omegac,epsi1,R1,z1,int_v)
    
#    Eref1 = Efield_ref_ana((megac,epsi1,epsi2,hbmu,hbgama,x1,z1,xe,zp,int_v)
    
    px = alffa_eff1_x*alffa_eff2_x*Edir1*(Gdir_12 + Gref_12 )

    
    return px

#%%

def dipole_moment_z_num(omegac,epsi1,epsi2,hbmu,hbgama,int_v,zp,x1,z1,x2,z2,omega01,kappa_factor_omega01,kappa_r_factor1,omega02,kappa_factor_omega02,kappa_r_factor2):     
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

    alffa1 = alpha_function(epsi1,omegac,omega01,kappa_factor_omega01,kappa_r_factor1)
    alffa2 = alpha_function(epsi1,omegac,omega02,kappa_factor_omega02,kappa_r_factor2)
    
    rtaself_x, rtaself_y, rtaself_z  =  green_self_pole_aprox(omegac,epsi1,epsi2,hbmu,hbgama,zp)
    alffa_eff1_x = (1/alffa1 -  rtaself_x)**(-1)

    alffa_eff2_x = (1/alffa2 -  rtaself_x)**(-1)


    phi = 0
    
    R = np.abs(x1-x2)
    Gdir_12 = green_tensor_ANA_QE(epsi1,R,z1,z2)

    Gref_12 = green_tensor_ref_pole_aprox(omegac,epsi1,epsi2,hbmu,hbgama,R,phi,z1,z2,zp)

    R1 = np.abs(x1)
    Edir1 = Efield_ANA(omegac,epsi1,R1,z1,int_v)

    px = alffa_eff1_x*alffa_eff2_x*Edir1*(Gdir_12 + Gref_12 )

    return px

#%%