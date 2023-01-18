#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila


En esta version, los limites de integracion son alrededor de cada polo

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
path_ctes =  path_basic.replace('/' + 'E_field','')
#print('Importar modulos necesarios para este codigo')
path_load = path_ctes
try:
    sys.path.insert(1, path_ctes)
    from constants import constantes
except ModuleNotFoundError:
    print('constants.py no se encuentra en ' + path_ctes)


try:
    sys.path.insert(1, path_ctes)
    from graphene_sigma import sigma_DL
except ModuleNotFoundError:
    print('graphene_sigma.py no se encuentra en ' + path_basic)

pi,hb,c,alfac,mu1,mu2 = constantes()
aux = c*hb

epsi2 = 1 ### los polos fueron hallados para este valor de epsilon2 
epsi1 = 1

epsilon_pole = 1e-8

#%%

def Ex_polp_num(omegac,hbmu,hbgama,R,phi,z,zp,px,py,pz,u):     
    """    
    Parameters
    ----------
    omegac : omega/c = k0 en 1/micrometros    
    d
    R
    phi
    z : coordenada z
    zp : posicion del plano (>0)
    px : coordenada x del dipolo 
    py : coordenada y del dipolo
    pz : coordenada z del dipolo
    Returns
    -------
    Ex pol s
    """
    E = omegac*aux
    k1 = omegac
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)      
    
    alpha_z1 = 1j*np.sqrt(u**2 - 1)
        
    z_dip_barra = k1*(2*zp - z)
    
    exp_electron_f = np.exp(1j*alpha_z1(u)*z_dip_barra)
    J0 = special.jv(0,u*k1*R)  
    J1 = special.jv(1,u*k1*R)  
    J2 = special.jv(2,u*k1*R)
    
    
    rp_num = epsi2*1j*u - epsi1*1j*u - cond*u**2
    rp_den =  epsi2*1j*u + epsi1*1j*u - cond*u**2
    rp =  rp_num(u)/rp_den(u)   
    
    
    term_Ex_px =  alpha_z1(u)*(-J0(u) + J2(u)*np.cos(2*phi))
    term_Ex_py =  alpha_z1(u)*J2(u)*np.sin(2*phi)
    
    
    termEx_num_re_f = np.real(u*rp(u)*(px*term_Ex_px(u) + py*term_Ex_py(u) + 2*1j*pz*u*J1(u)*np.cos(phi))*exp_electron_f(u))
    termEx_num_im_f = np.imag(u*rp(u)*(px*term_Ex_px(u) + py*term_Ex_py(u) + 2*1j*pz*u*J1(u)*np.cos(phi))*exp_electron_f(u))
  
    

    termEx_num_int = termEx_num_re_f + 1j*termEx_num_im_f


    return 1j*0.5*termEx_num_int*(k1**3)  

#%%

def Ex_polp_pole_aprox(omegac,hbmu,hbgama,R,phi,z,zp,px,py,pz,u):     
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
    k1 = omegac
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)      
    
    alpha_z1 = 1j*np.sqrt(u**2 - 1)
        
    z_dip_barra = k1*(2*zp - z)
    
    exp_electron_f = np.exp(1j*alpha_z1(u)*z_dip_barra)
    J0 =  special.jv(0,u*k1*R)  
    J1 =   special.jv(1,u*k1*R)  
    J2 =   special.jv(2,u*k1*R)
    
    
    Rp = 2*epsi1/(epsi1 + epsi2)
    alfa_p = 1j*(epsi1 + epsi2)/(cond)
    rp = Rp*alfa_p/(u - alfa_p)  
    
    
    term_Ex_px =  alpha_z1(u)*(-J0(u) + J2(u)*np.cos(2*phi))
    term_Ex_py =  alpha_z1(u)*J2(u)*np.sin(2*phi)
    
    
    termEx_num_re_f = np.real(u*rp(u)*(px*term_Ex_px(u) + py*term_Ex_py(u) + 2*1j*pz*u*J1(u)*np.cos(phi))*exp_electron_f(u))
    termEx_num_im_f = np.imag(u*rp(u)*(px*term_Ex_px(u) + py*term_Ex_py(u) + 2*1j*pz*u*J1(u)*np.cos(phi))*exp_electron_f(u))
          

    termEx_num_int = termEx_num_re_f + 1j*termEx_num_im_f


    return 1j*0.5*termEx_num_int*(k1**3)  


#%%


#%%

def Ey_polp_num(omegac,hbmu,hbgama,R,phi,z,zp,px,py,pz,u):     
    """    
    Parameters
    ----------
    omegac : omega/c = k0 en 1/micrometros    
    d
    R
    phi
    z : coordenada z
    zp : posicion del plano (>0)
    px : coordenada x del dipolo 
    py : coordenada y del dipolo
    pz : coordenada z del dipolo
    Returns
    -------
    Ex pol s
    """

    E = omegac*aux
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    k1 = omegac
    
    
    alpha_z1 = 1j*np.sqrt(u**2 - 1)

        
    z_dip_barra = k1*(2*zp - z)
    
    exp_electron_f = np.exp(1j*alpha_z1(u)*z_dip_barra)
    J0 =  special.jv(0,u*k1*R)  
    J1 =   special.jv(1,u*k1*R)  
    J2 =   special.jv(2,u*k1*R)
    
    
    rp_num = epsi2*1j*u - epsi1*1j*u - cond*u**2
    rp_den =  epsi2*1j*u + epsi1*1j*u - cond*u**2
    rp =  rp_num(u)/rp_den(u)   
    
    
    term_Ex_px =  -alpha_z1(u)*(J0(u) + J2(u)*np.cos(2*phi))
    term_Ex_py = alpha_z1(u)*J2(u)*np.sin(2*phi)
    
    
    termEx_num_re_f =  np.real(u*rp(u)*(py*term_Ex_px(u) + px*term_Ex_py(u) + 2*1j*pz*u*J1(u)*np.cos(phi))*exp_electron_f(u))
    termEx_num_im_f =  np.imag(u*rp(u)*(py*term_Ex_px(u) + px*term_Ex_py(u) + 2*1j*pz*u*J1(u)*np.cos(phi))*exp_electron_f(u))

    termEx_num_int = termEx_num_re_f + 1j*termEx_num_im_f

    return 1j*0.5*termEx_num_int*(k1**3)  

#%%

def Ey_polp_pole_aprox(omegac,hbmu,hbgama,R,phi,z,zp,px,py,pz,u):     
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
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    k1 = omegac
    
    alpha_z1 =1j*np.sqrt(u**2 - 1)
        
    z_dip_barra = k1*(2*zp - z)
    
    exp_electron_f = np.exp(1j*alpha_z1(u)*z_dip_barra)
    J0 =  special.jv(0,u*k1*R)  
    J1 = special.jv(1,u*k1*R)
    J2 =   special.jv(2,u*k1*R)

    Rp = 2*epsi1/(epsi1 + epsi2)
    alfa_p = 1j*(epsi1 + epsi2)/(cond)
    rp =  Rp*alfa_p/(u - alfa_p)  
    
    term_Ex_px =  -alpha_z1(u)*(J0(u) + J2(u)*np.cos(2*phi))
    term_Ex_py =  alpha_z1(u)*J2(u)*np.sin(2*phi)
    
    
    termEx_num_re_f = np.real(u*rp(u)*(py*term_Ex_px(u) + px*term_Ex_py(u) + 2*1j*pz*u*J1(u)*np.cos(phi))*exp_electron_f(u))
    termEx_num_im_f =  np.imag(u*rp(u)*(py*term_Ex_px(u) + px*term_Ex_py(u) + 2*1j*pz*u*J1(u)*np.cos(phi))*exp_electron_f(u))


    termEy_pole_aprox_int = termEx_num_re_f + 1j*termEx_num_im_f

    return 1j*0.5*termEy_pole_aprox_int*(k1**3)  


#%%
    
def Ez_polp_pole_aprox(omegac,hbmu,hbgama,R,phi,z,zp,px,py,pz,u):     
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
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    k1 = omegac
    
    alpha_z1 =  1j*np.sqrt(u**2 - 1)
        
    z_dip_barra = k1*(2*zp - z)
    
    exp_electron_f = np.exp(1j*alpha_z1(u)*z_dip_barra)
    J0 =   special.jv(0,u*k1*R)  
    J1 =  special.jv(1,u*k1*R)

    Rp = 2*epsi1/(epsi1 + epsi2)
    alfa_p = 1j*(epsi1 + epsi2)/(cond)
    rp = Rp*alfa_p/(u - alfa_p)   ## u es alfa_parallel y la funcion depende de k_parallel
    
    term_Ex = px*np.cos(phi) + py*np.sin(phi) 
    
    
    termEx_num_re_f = np.real(u**2*rp(u)*(J1(u)*term_Ex + 1j*pz*J0(u)*u/alpha_z1(u))*exp_electron_f(u))
    termEx_num_im_f =np.imag(u**2*rp(u)*(J1(u)*term_Ex + 1j*pz*J0(u)*u/alpha_z1(u))*exp_electron_f(u))


    termEz_pole_aprox_int = termEx_num_re_f + 1j*termEx_num_im_f
    return termEz_pole_aprox_int*(k1**3)  


#%%


def Ez_polp_num(omegac,hbmu,hbgama,R,phi,z,zp,px,py,pz,u):     
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
    k1 = omegac
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)   
    
    
    alpha_z1 = 1j*np.sqrt(u**2 - 1)     
    z_dip_barra = k1*(2*zp - z)
    
    exp_electron_f =  np.exp(1j*alpha_z1(u)*z_dip_barra)
    J0 = special.jv(0,u*k1*R)  
    J1 =  special.jv(1,u*k1*R)

    rp_num = epsi2*1j*u - epsi1*1j*u - cond*u**2
    rp_den = epsi2*1j*u + epsi1*1j*u - cond*u**2
    rp = rp_num(u)/rp_den(u) 

    
    term_Ex = px*np.cos(phi) + py*np.sin(phi) 
    
    
    termEx_num_re_f = np.real(u**2*rp(u)*(J1(u)*term_Ex + 1j*pz*J0(u)*u/alpha_z1(u))*exp_electron_f(u))
    termEx_num_im_f =  np.imag(u**2*rp(u)*(J1(u)*term_Ex + 1j*pz*J0(u)*u/alpha_z1(u))*exp_electron_f(u))
    

    termEy_pole_aprox_int = termEx_num_re_f + 1j*termEx_num_im_f


    return termEy_pole_aprox_int*(k1**3)  


#%%
    

 