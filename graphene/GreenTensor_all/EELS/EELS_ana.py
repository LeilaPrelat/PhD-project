#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila

campo externo directo analitico con la convencion de z hacia abajo
en el punto x  = 0, y = 0,z = 0
"""
import numpy as np
import sys
import os 
from scipy import special

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/EELS','')
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

def Efield_dir_ana_versionL(omegac,epsi1,int_v,b):
    """
    Parameters
    ----------
    omegac : omega/c = k0 en 1/micrometros    
    epsi1 : epsilon del medio de arriba del plano
    int_v : fraccion de la vel de la c que vale la velocidad uniforme del electron en x, si v = c/int
    b : electron position en z, b < 0 (dipolo en 0 y eje de z apunta hacia abajo)  
    Returns
    -------
    External field direct (green tensor direct)
    analitico luego de aplicar QE y de despreciar
    2 terminos (rosa y violeta)
    """

 #   charge_electron_c = 1j*(1.602*1e-19*c)*omegac/(2*np.pi)
    
    # x_y = ky/k0
    n1 = epsi1*mu1
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1
    k1_2 = k1**2

    arg = np.abs(b)*omegac*int_v
    K0 = special.kn(0,arg)

    cte_general = -(omegac*int_v)**2 + k1_2
    rta = cte_general*K0
    
    return rta


#%%
    

def Efield_dir_ana_version149(omegac,epsi1,int_v,b):
    """
    Parameters
    ----------
    omegac : omega/c = k0 en 1/micrometros    
    epsi1 : epsilon del medio de arriba del plano
    int_v : fraccion de la vel de la c que vale la velocidad uniforme del electron en x, si v = c/int
    b : electron position en z, b < 0 (dipolo en 0 y eje de z apunta hacia abajo)  
    Returns
    -------
    External field direct version paper 149
    """

#    charge_electron_c = 2*1j*(1.602*1e-19*c)*omegac
    
    # x_y = ky/k0
#    n1 = epsi1*mu1
#    cte1 = np.sqrt(n1)
#    k1 = omegac*cte1
#    k1_2 = k1**2
    v = c/int_v
    v2 = v**2
    
    arg = np.abs(b)*omegac*int_v
    K0 = special.kn(0,arg)
    # cte = omega/(v*np.abs(b))
    
    rta = K0/(v2*epsi1)
    
    return rta


#%%

def Efield_ref_ana1b(omegac,epsi1,epsi2,hbmu,hbgama,int_v,b,zp,z):
    """
    Parameters
    ----------
    omegac : omega/c = k0 en 1/micrometros    
    epsi1 : epsilon del medio de arriba del plano
    int_v : fraccion de la vel de la c que vale la velocidad uniforme del electron en x, si v = c/int
    b: electron position en z, b < 0 (dipolo en 0 y eje de z apunta hacia abajo) 
    x : punto x donde miramos el campo en micrones
    y : punto y donde miramos el campo en micrones
    z : punto z donde miramos el campo en micrones 
    Returns
    -------
    External field direct (green tensor direct)
    dividido por 1j*e/(2*pi)
    analitico opcion 1.b
    """

#    cte_aux = (1.602*1e-19*c)/int_v
    
    E = omegac*aux  
 #   omega = omegac*c
#    k0 = omegac
    n1 = epsi1*mu1
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1
    k1_2 = k1**2
    n_v1 = int_v/cte1
    pi2 = np.pi**2
    
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    Rp = 2*epsi1/(epsi1 + epsi2)
    alfa_p = 1j*(epsi1 + epsi2)/(cond*cte1)

    exp_electron = np.exp(-alfa_p*k1*(np.abs(z) + np.abs(b) + 2*zp))
    
    cte = pi2*k1_2*Rp*(alfa_p**(5/2))
    cte2 = 1/np.sqrt(alfa_p + n_v1) + 1/np.sqrt(alfa_p - n_v1) 
    cte3 = (1 - 1j)*np.sqrt(-1j)
    
    final_expression = cte*cte2*cte3*exp_electron
    
    return final_expression

#%%
    

def green_tensor_self_ana1_xx(omegac,epsi1,epsi2,hbmu,hbgama,xD,zD,zp):
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
    zD : coordenada z del dipolo 
    zp : posicion del plano (>0)
    Returns
    -------
    Gxx self (superficie)
    con z hacia abajo (convencion del paper)
    numerico con la aprox QE (no hay alpha_z)
    + p polarization approx (fresnel coefficient)
    solo para x,y,z = xD,yD,zD = 0,0,0
    Aprox analitica
    """

    E = omegac*aux
    k0 = omegac #=omega/c
    # x_y = ky/k0
    n1 = epsi1*mu1
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1
    k1_3 = (k0*cte1)**3
    
    z_dip_barra = k1*(np.abs(zD) + 2*zp + np.abs(zD))
    Rself = np.sqrt(2)*k1*xD

    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)

    alfa_p = (k0/k1)*1j*(epsi1 + epsi2)/cond
    alfap_3 = alfa_p**3
 
    Rp = 2*epsi1/(epsi1 + epsi2)
    
    expB = np.exp(-alfa_p*z_dip_barra) 
    
    
    arg = alfa_p*Rself
    J0 = special.jn(0,arg)
    J2 = special.jn(2,arg)
    
    
    final = np.pi*k1_3*Rp*alfap_3*(J0 + J2)
    
    return final*expB

#%%
    

def green_tensor_self_ana1_zz(omegac,epsi1,epsi2,hbmu,hbgama,xD,zD,zp):
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
    zD : coordenada z del dipolo 
    zp : posicion del plano (>0)
    Returns
    -------
    Gxx self (superficie)
    con z hacia abajo (convencion del paper)
    numerico con la aprox QE (no hay alpha_z)
    + p polarization approx (fresnel coefficient)
    solo para x,y,z = xD,yD,zD = 0,0,0
    Aprox analitica
    """

    E = omegac*aux
    k0 = omegac #=omega/c
    # x_y = ky/k0
    n1 = epsi1*mu1
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1
    k1_3 = (k0*cte1)**3
    
    z_dip_barra = k1*(np.abs(zD) + 2*zp + np.abs(zD))
    Rself = np.sqrt(2)*k1*xD

    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)

    alfa_p = (k0/k1)*1j*(epsi1 + epsi2)/cond
    alfap_3 = alfa_p**3
 
    Rp = 2*epsi1/(epsi1 + epsi2)
    
    expB = np.exp(-alfa_p*z_dip_barra) 
    
    
    arg = alfa_p*Rself
    J0 = special.jn(0,arg)

    
    final = 2*np.pi*1j*k1_3*Rp*alfap_3*J0 
    
    return final*expB

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
    
def EELS_f(omegac,epsi1,epsi2,hbmu,hbgama,z,xD,zD,b,zp,int_v,omega0,kappa_factor_omega0,kappa_r_factor,px,py,pz):
    
    
    alffa = alpha_function(epsi1,omegac,omega0,kappa_factor_omega0,kappa_r_factor)
    Gself_xx = green_tensor_self_ana1_xx(omegac,epsi1,epsi2,hbmu,hbgama,z,xD,zD,zp)
    Gself_zz = green_tensor_self_ana1_zz(omegac,epsi1,epsi2,hbmu,hbgama,z,xD,zD,zp)
    
    Gself = px*Gself_xx + py*Gself_xx + pz*Gself_zz
    
    Edir = Efield_dir_ana_versionL(omegac,epsi1,int_v,b)
    Eref = Efield_ref_ana1b(omegac,epsi1,epsi2,hbmu,hbgama,int_v,b,zp,z)
    
    Etot = (Edir + Eref)**2
    
    den = 1/alffa - Gself  
    
    pi2 = np.pi**2
    
    aux_final = alfac*c/(pi2*int_v)
    
    return aux_final*Etot/den

#%%

    