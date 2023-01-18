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
path_constants =  path_basic.replace('/potential_field/potential_and_electric_field/potential_and_electric_field_simple_version','')
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
limitt = 60

#%%

def G1_ana(omegac,epsi1,epsi2,hbmu,hbgama,R,z,zp):     
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
    k1_2 = k1**2
 #   n_v1 = int_v/cte1
#    k1_2 = k1**2
   
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    Rp = 2*epsi1/(epsi1 + epsi2)
    alfa_p = 1j*(epsi1 + epsi2)/(cond*cte1)
    kp = alfa_p*k1
    kp_2 = kp**2 
    
    z_dip_barra = k1*(-z + 2*zp)        
    exp_electron = np.exp(-alfa_p*z_dip_barra)
     

    H0 = special.hankel1(0,kp*R)
    
    term1 = Rp*np.pi*1j*kp_2*H0*exp_electron

    
###################################################################################################

    return term1 


#%%

def G1_num(omegac,epsi1,epsi2,hbmu,hbgama,R,z,zp):     
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
    k1_2 = k1**2
 #   n_v1 = int_v/cte1
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)  

    
#    cota_inf = 0.01
#    cota_sup = 1000/k1
    
    z_dip_barra = 2*zp - z
    
    exp_electron_f = lambda u: np.exp(-u*k1*z_dip_barra)
    J0 =  lambda u : special.jv(0,u*k1*R)


    rp_num = lambda u: epsi2*1j*u - epsi1*1j*u - cte1*cond*u**2
    rp_den = lambda u: epsi2*1j*u + epsi1*1j*u - cte1*cond*u**2
    rp = lambda u: rp_num(u)/rp_den(u)   
    
    term2_re_f = lambda u: np.real(J0(u)*u*rp(u)*exp_electron_f(u))
    term2_im_f = lambda u: np.imag(J0(u)*u*rp(u)*exp_electron_f(u))
  
    
    term2_int0 = 0
    
    cota_alpha_y = 600
    nmax = int(cota_alpha_y*R*omegac/np.pi - 0.5)
    for n in range(0,nmax+1): 

        cota_sup = (2*n+1)*np.pi*0.5/(omegac*R)
        if n == 0:
            cota_inf = 0
        else:
            cota_inf = (2*(n-1)+1)*np.pi*0.5/(omegac*R)
        
        term2_int_re,err = integrate.quad(term2_re_f, cota_inf, cota_sup, limit = limitt) 
        term2_int_im,err = integrate.quad(term2_im_f, cota_inf, cota_sup, limit = limitt) 
    
        term2_int = term2_int_re + 1j*term2_int_im
    
        term2_int0 = term2_int0 + term2_int
    
    return term2_int0*k1_2 

#%%

def G1_pole_aprox(omegac,epsi1,epsi2,hbmu,hbgama,R,z,zp):     
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
    k1_2 = k1**2
 #   n_v1 = int_v/cte1
#    k1_2 = k1**2
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    Rp = 2*epsi1/(epsi1 + epsi2)
    alfa_p = 1j*(epsi1 + epsi2)/(cond*cte1)
    kp = alfa_p*k1
    
    cota_inf = 0.01/k1
    cota_sup = 600/k1
    
    z_dip_barra = 2*zp - z
    
    exp_electron_f = lambda u: np.exp(-u*k1*z_dip_barra)
    J0 = lambda u : special.jv(0,u*k1*R)

    rp = lambda u: Rp*alfa_p/(u - alfa_p)   
    
    term2_re_f = lambda u: np.real(J0(u)*u*rp(u)*exp_electron_f(u))
    term2_im_f = lambda u: np.imag(J0(u)*u*rp(u)*exp_electron_f(u))
    
    term2_int_re,err = integrate.quad(term2_re_f, cota_inf, cota_sup, limit = limitt) 
    term2_int_im,err = integrate.quad(term2_im_f, cota_inf, cota_sup, limit = limitt) 

    term2_int = term2_int_re + 1j*term2_int_im
    
    
    term2_final = term2_int
    
    return term2_final*k1_2


#%%
    
def G2_ana(omegac,epsi1,epsi2,hbmu,hbgama,R,z,zp):     
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
    k1_2 = k1**2
 #   n_v1 = int_v/cte1
#    k1_2 = k1**2
   
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    Rp = 2*epsi1/(epsi1 + epsi2)
    alfa_p = 1j*(epsi1 + epsi2)/(cond*cte1)
    kp = alfa_p*k1
    kp_2 = kp**2 
    
    z_dip_barra = k1*(-z + 2*zp)        
    exp_electron = np.exp(-alfa_p*z_dip_barra)
     

    H1 = special.hankel1(1,kp*R)
    
    term1 = Rp*np.pi*1j*kp_2*H1*exp_electron

    
###################################################################################################

    return term1 



#%%

def G2_num(omegac,epsi1,epsi2,hbmu,hbgama,R,z,zp):     
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
    k1_2 = k1**2
 #   n_v1 = int_v/cte1
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)


#    cota_inf = 0.01
#    cota_sup = 1000/k1
    
    z_dip_barra = 2*zp - z
    
    exp_electron_f = lambda u: np.exp(-u*k1*z_dip_barra)
    J1 =  lambda u : special.jv(1,u*k1*R)


    rp_num = lambda u: epsi2*1j*u - epsi1*1j*u - cte1*cond*u**2
    rp_den = lambda u: epsi2*1j*u + epsi1*1j*u - cte1*cond*u**2
    rp = lambda u: rp_num(u)/rp_den(u)   
    
    term2_re_f = lambda u: np.real(J1(u)*u*rp(u)*exp_electron_f(u))
    term2_im_f = lambda u: np.imag(J1(u)*u*rp(u)*exp_electron_f(u))


    term2_int0 = 0
    
    cota_alpha_y = 600
    nmax = int(cota_alpha_y*R*omegac/np.pi - 0.5)
    for n in range(0,nmax+1): 

        cota_sup = (2*n+1)*np.pi*0.5/(omegac*R)
        if n == 0:
            cota_inf = 0
        else:
            cota_inf = (2*(n-1)+1)*np.pi*0.5/(omegac*R)
 
        term2_int_re,err = integrate.quad(term2_re_f, cota_inf, cota_sup, limit = limitt) 
        term2_int_im,err = integrate.quad(term2_im_f, cota_inf, cota_sup, limit = limitt) 
    
        term2_int = term2_int_re + 1j*term2_int_im
    
        term2_int0 = term2_int0 + term2_int
    
    
    return term2_int0*k1_2 

#%%

def G2_pole_aprox(omegac,epsi1,epsi2,hbmu,hbgama,R,z,zp):     
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
    k1_2 = k1**2
 #   n_v1 = int_v/cte1

    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    Rp = 2*epsi1/(epsi1 + epsi2)
    alfa_p = 1j*(epsi1 + epsi2)/(cond*cte1)
    kp = alfa_p*k1
    
    cota_inf = 0.01/k1
    cota_sup = 600/k1
    
    z_dip_barra = 2*zp - z
    
    exp_electron_f = lambda u: np.exp(-u*k1*z_dip_barra)
    J1 = lambda u : special.jv(1,u*k1*R)

    rp = lambda u: Rp*alfa_p/(u - alfa_p)   
    
    term2_re_f = lambda u: np.real(J1(u)*u*rp(u)*exp_electron_f(u))
    term2_im_f = lambda u: np.imag(J1(u)*u*rp(u)*exp_electron_f(u))
    
    term2_int_PP_re,err = integrate.quad(term2_re_f, cota_inf, cota_sup, limit = limitt) 
    term2_int_PP_im,err = integrate.quad(term2_im_f, cota_inf, cota_sup, limit = limitt) 

    termG2_int = term2_int_PP_re + 1j*term2_int_PP_im
    
    
    term2_final = termG2_int
    
    return term2_final*k1_2  



#%%
    

def potential_final_ana(omegac,epsi1,epsi2,hbmu,hbgama,phi,R,z,zp,px,py,pz):
    
    G1 =  G1_ana(omegac,epsi1,epsi2,hbmu,hbgama,R,z,zp)
    
    G2 =  G2_ana(omegac,epsi1,epsi2,hbmu,hbgama,R,z,zp)
    
    
    term0_3 = (np.abs(z)**2 + R**2)**(3/2)
    term0_1 = (np.abs(z)**2 + R**2)**(1/2)
    
    px_py_term = px*np.cos(phi) + py*np.sin(phi)
    
    term1 = -px_py_term*(G2 + (np.abs(z)**2/term0_3 - 1/term0_1)/R)
    
    term2 = -pz*np.sign(z)*(G1 + np.abs(z)/term0_3)
    
    return term1 + term2


#%%s


def potential_final_pole_aprox(omegac,epsi1,epsi2,hbmu,hbgama,phi,R,z,zp,px,py,pz):
    
    G1 =  G1_pole_aprox(omegac,epsi1,epsi2,hbmu,hbgama,R,z,zp)
    
    G2 =  G2_pole_aprox(omegac,epsi1,epsi2,hbmu,hbgama,R,z,zp)
    
    
    term0_3 = (np.abs(z)**2 + R**2)**(3/2)
    term0_1 = (np.abs(z)**2 + R**2)**(1/2)
    
    px_py_term = px*np.cos(phi) + py*np.sin(phi)
    
    term1 = -px_py_term*(G2 + (np.abs(z)**2/term0_3 - 1/term0_1)/R)
    
    term2 = -pz*np.sign(z)*(G1 + np.abs(z)/term0_3)
    
    return term1 + term2

 #%%s   
 
def potential_final_num(omegac,epsi1,epsi2,hbmu,hbgama,phi,R,z,zp,px,py,pz):
    
    G1 =  G1_num(omegac,epsi1,epsi2,hbmu,hbgama,R,z,zp)
    
    G2 =  G2_num(omegac,epsi1,epsi2,hbmu,hbgama,R,z,zp)
    
    term0_3 = (np.abs(z)**2 + R**2)**(3/2)
    term0_1 = (np.abs(z)**2 + R**2)**(1/2)
    
    px_py_term = px*np.cos(phi) + py*np.sin(phi)
    
    term1 = -px_py_term*(G2 + (np.abs(z)**2/term0_3 - 1/term0_1)/R)
    
    term2 = -pz*np.sign(z)*(G1 + np.abs(z)/term0_3)
    
    return term1 + term2 
 
#%%