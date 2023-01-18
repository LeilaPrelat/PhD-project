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
from scipy import optimize

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

limitt = 200

#%%

#plot_vs_R = 1
#plot_vs_E = 0

#def zeroEx_num(omegac,hbmu,hbgama,R,phi,z,zp,px,py,pz):
#
#    k1 = omegac
#    
#    def function(omegac,R,u):
#    
#        alpha_z1 = 1j*np.sqrt(u**2 - 1)
#        
#        J0 =  special.jv(0,u*k1*R)  
#        J1 =  special.jv(1,u*k1*R)  
#        J2 =  special.jv(2,u*k1*R)
#        
#        eq = alpha_z1*(px*(J2(u)*np.cos(2*phi)-J0(u)) + py*J2(u)*np.sin(2*phi)) + 2*1j*pz*u*np.cos(phi)*J1
#    
#        return eq
#    
#    
#    deltta = 50
#    def raiz(omegac,R,u):
#        if plot_vs_R == 1:
#            R1 = R - deltta
#            R2 = R + deltta
#        
#            eq1 = function(omegac,R1,u)
#            eq2 = function(omegac,R2,u)
#            
#            if eq1*eq2 < 0:
#                
#                def function_raiz(u):
#                    return raiz(omegac,R,u)
#                
#                init_guess = np.mean([R1,R2])
#                sol = optimize.root(function_raiz, [init_guess])
#            
#                listy = np.array([function(omegac,R1,sol.x), function(omegac,R2,sol.x)])
#                listy = np.sorted(listy)
#                if sol.x in listy:
#                    return sol.x
#            
#        elif plot_vs_E == 1:
#        
#            omegac1 = omegac - deltta
#            omegac2 = omegac + deltta
#            
#            eq1 = function(omegac1,R,u)
#            eq2 = function(omegac2,R,u)
#            
#            if eq1*eq2 < 0:
#
#                def function_raiz(u):
#                    return raiz(omegac,R,u)       
#                
#                init_guess = np.mean([omegac1,omegac2]) 
#                sol = optimize.root(function_raiz, [init_guess])
#                
#                listy = np.array([function(omegac1,R,sol.x), function(omegac2,R,sol.x)])
#                listy = np.sorted(listy)
#                if sol.x in listy:
#                    return sol.x
            
               
#%%

def Ex_polp_num(omegac,hbmu,hbgama,R,phi,z,zp,px,py,pz):     
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
    
    alpha_z1 =  lambda u : np.sqrt(1 - u**2) if (1 > np.abs(u)) else 1j*np.sqrt(u**2 -1)   
        
    z_dip_barra = k1*(2*zp - z)
    
    exp_electron_f = lambda u: np.exp(1j*alpha_z1(u)*z_dip_barra)
    J0 =  lambda u : special.jv(0,u*k1*R)  
    J1 =  lambda u : special.jv(1,u*k1*R)  
    J2 =  lambda u : special.jv(2,u*k1*R)
    
    
    rp_num = lambda u: epsi2*1j*u - epsi1*1j*u - cond*u**2
    rp_den = lambda u: epsi2*1j*u + epsi1*1j*u - cond*u**2
    rp = lambda u: rp_num(u)/rp_den(u)   
    
    
    term_Ex_px = lambda u: alpha_z1(u)*(J2(u)*np.cos(2*phi) - J0(u))
    term_Ex_py = lambda u: alpha_z1(u)*J2(u)*np.sin(2*phi)
    
    
    termEx_num_re_f = lambda u: np.real(u*rp(u)*(px*term_Ex_px(u) + py*term_Ex_py(u) + 2*1j*pz*u*J1(u)*np.cos(phi))*exp_electron_f(u))
    termEx_num_im_f = lambda u: np.imag(u*rp(u)*(px*term_Ex_px(u) + py*term_Ex_py(u) + 2*1j*pz*u*J1(u)*np.cos(phi))*exp_electron_f(u))
  
#    cota_inf = 1
#    cota_sup = 600*k1
    
    term2_int0 = 0
    
    cota_alpha_y = 600
    nmax = int(cota_alpha_y*R*omegac/np.pi - 0.5)
    for n in range(0,nmax+1): 

        cota_sup = (2*n+1)*np.pi*0.5/(omegac*R)
        if n == 0:
            cota_inf = 0
        else:
            cota_inf = (2*(n-1)+1)*np.pi*0.5/(omegac*R)
 
        term2_int_re,err = integrate.quad(termEx_num_re_f, cota_inf, cota_sup, limit = 50) 
        term2_int_im,err = integrate.quad(termEx_num_im_f, cota_inf, cota_sup, limit = 50) 
    
        term2_int = term2_int_re + 1j*term2_int_im
    
        term2_int0 = term2_int0 + term2_int        

    
#    termEx_num_int_re,err = integrate.quad(termEx_num_re_f, cota_inf, cota_sup,limit = limitt) 
#    termEx_num_int_im,err = integrate.quad(termEx_num_im_f, cota_inf, cota_sup,limit = limitt) 
#
#    termEx_num_int = termEx_num_int_re + 1j*termEx_num_int_im


    return 1j*0.5*term2_int0*(k1**3)  

#%%

def Ex_polp_pole_aprox(omegac,hbmu,hbgama,R,phi,z,zp,px,py,pz):     
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
    
    alpha_z1 =  lambda u : np.sqrt(1 - u**2) if (1 > np.abs(u)) else 1j*np.sqrt(u**2 -1)   
        
    z_dip_barra = k1*(2*zp - z)
    
    exp_electron_f = lambda u: np.exp(1j*alpha_z1(u)*z_dip_barra)
    J0 =  lambda u : special.jv(0,u*k1*R)  
    J1 =  lambda u : special.jv(1,u*k1*R)  
    J2 =  lambda u : special.jv(2,u*k1*R)
    
    
    Rp = 2*epsi1/(epsi1 + epsi2)
    alfa_p = 1j*(epsi1 + epsi2)/(cond)
    rp = lambda u: Rp*alfa_p/(u - alfa_p)  
    
    
    term_Ex_px = lambda u: alpha_z1(u)*( J2(u)*np.cos(2*phi) - J0(u))
    term_Ex_py = lambda u: alpha_z1(u)*J2(u)*np.sin(2*phi)
    
    
    termEx_num_re_f = lambda u: np.real(u*rp(u)*(px*term_Ex_px(u) + py*term_Ex_py(u) + 2*1j*pz*u*J1(u)*np.cos(phi))*exp_electron_f(u))
    termEx_num_im_f = lambda u: np.imag(u*rp(u)*(px*term_Ex_px(u) + py*term_Ex_py(u) + 2*1j*pz*u*J1(u)*np.cos(phi))*exp_electron_f(u))
    
    cota_inf = 0.01
    cota_sup = 600
        

    
    termEx_num_int_re,err = integrate.quad(termEx_num_re_f, cota_inf, cota_sup) 
    termEx_num_int_im,err = integrate.quad(termEx_num_im_f, cota_inf, cota_sup) 

    termEx_num_int = termEx_num_int_re + 1j*termEx_num_int_im


    return 1j*0.5*termEx_num_int*(k1**3)  


#%%
    

def Ex_polp_ana(omegac,hbmu,hbgama,R,phi,z,zp,px,py,pz):     
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
    k0 = omegac
   

    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    Rp = 2*epsi1/(epsi1 + epsi2)
    alfa_p = 1j*(epsi1 + epsi2)/(cond)
    kp = alfa_p*k0

    
    alpha_z1 = 1j*np.sqrt(alfa_p**2 - 1)
        
    z_dip_barra = 2*zp - z
    
    exp_electron_f = np.exp(1j*alpha_z1*k0*z_dip_barra)
    H0 = special.hankel1(0,kp*R)  
    H1 = special.hankel1(1,kp*R)  
    H2 = special.hankel1(2,kp*R)

    rs = Rp*alfa_p

    term_Ex_px = alpha_z1*(-H0 + H2*np.cos(2*phi))
    term_Ex_py = alpha_z1*H2*np.sin(2*phi)        

    termEx = alfa_p*rs*(px*term_Ex_px + py*term_Ex_py + 2*1j*pz*H1*np.cos(phi)*alfa_p)*exp_electron_f

        
    return -np.pi*0.5*termEx*(k0**3)  
 

#%%
    
def Ey_polp_ana(omegac,hbmu,hbgama,R,phi,z,zp,px,py,pz):     
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
    k0 = omegac

    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    Rp = 2*epsi1/(epsi1 + epsi2)
    alfa_p = 1j*(epsi1 + epsi2)/(cond)
    kp = alfa_p*k0

    
    alpha_z1 = np.sqrt(1 - alfa_p**2) if (1 > np.abs(alfa_p)) else 1j*np.sqrt(alfa_p**2 -1)   
        
    z_dip_barra = 2*zp - z
    
    exp_electron_f = np.exp(1j*alpha_z1*k0*z_dip_barra)
    H0 = special.hankel1(0,kp*R)  
    H1 = special.hankel1(1,kp*R)  
    H2 = special.hankel1(2,kp*R)

    rp = Rp*alfa_p
    
    term_Ex_px = -alpha_z1*(H0 + H2*np.cos(2*phi))
    term_Ex_py = alpha_z1*H2*np.sin(2*phi)     
    

    termEx = alfa_p*rp*(py*term_Ex_px + px*term_Ex_py + 2*1j*pz*H1*np.cos(phi)*alfa_p)*exp_electron_f
        
    return -np.pi*0.5*termEx*(k0**3)  
 

#%%

def Ey_polp_num(omegac,hbmu,hbgama,R,phi,z,zp,px,py,pz):     
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
    
    
    alpha_z1 =  lambda u : np.sqrt(1 - u**2) if (1 > np.abs(u)) else 1j*np.sqrt(u**2 -1)   

        
    z_dip_barra = k1*(2*zp - z)
    
    exp_electron_f = lambda u: np.exp(1j*alpha_z1(u)*z_dip_barra)
    J0 =  lambda u : special.jv(0,u*k1*R)  
    J1 =  lambda u : special.jv(1,u*k1*R)  
    J2 =  lambda u : special.jv(2,u*k1*R)
    
    
    rp_num = lambda u: epsi2*1j*u - epsi1*1j*u - cond*u**2
    rp_den = lambda u: epsi2*1j*u + epsi1*1j*u - cond*u**2
    rp = lambda u: rp_num(u)/rp_den(u)   
    
    
    term_Ex_px = lambda u: -alpha_z1(u)*(J0(u) + J2(u)*np.cos(2*phi))
    term_Ex_py = lambda u: alpha_z1(u)*J2(u)*np.sin(2*phi)
    
    
    termEx_num_re_f = lambda u: np.real(u*rp(u)*(py*term_Ex_px(u) + px*term_Ex_py(u) + 2*1j*pz*u*J1(u)*np.cos(phi))*exp_electron_f(u))
    termEx_num_im_f = lambda u: np.imag(u*rp(u)*(py*term_Ex_px(u) + px*term_Ex_py(u) + 2*1j*pz*u*J1(u)*np.cos(phi))*exp_electron_f(u))

    
    cota_inf = 0.01
    cota_sup = 600
            
    termEy_num_int_re,err = integrate.quad(termEx_num_re_f, cota_inf, cota_sup,limit = limitt) 
    termEy_num_int_im,err = integrate.quad(termEx_num_im_f, cota_inf, cota_sup,limit = limitt) 

    termEx_num_int = termEy_num_int_re + 1j*termEy_num_int_im

    return 1j*0.5*termEx_num_int*(k1**3)  

#%%

def Ey_polp_pole_aprox(omegac,hbmu,hbgama,R,phi,z,zp,px,py,pz):     
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
    
    alpha_z1 =  lambda u : np.sqrt(1 - u**2) if (1 > np.abs(u)) else 1j*np.sqrt(u**2 -1)   
        
    z_dip_barra = k1*(2*zp - z)
    
    exp_electron_f = lambda u: np.exp(1j*alpha_z1(u)*z_dip_barra)
    J0 =  lambda u : special.jv(0,u*k1*R)  
    J1 =  lambda u : special.jv(1,u*k1*R)
    J2 =  lambda u : special.jv(2,u*k1*R)

    Rp = 2*epsi1/(epsi1 + epsi2)
    alfa_p = 1j*(epsi1 + epsi2)/(cond)
    rp = lambda u: Rp*alfa_p/(u - alfa_p)  
    
    term_Ex_px = lambda u: -alpha_z1(u)*(J0(u) + J2(u)*np.cos(2*phi))
    term_Ex_py = lambda u: alpha_z1(u)*J2(u)*np.sin(2*phi)
    
    
    termEx_num_re_f = lambda u: np.real(u*rp(u)*(py*term_Ex_px(u) + px*term_Ex_py(u) + 2*1j*pz*u*J1(u)*np.cos(phi))*exp_electron_f(u))
    termEx_num_im_f = lambda u: np.imag(u*rp(u)*(py*term_Ex_px(u) + px*term_Ex_py(u) + 2*1j*pz*u*J1(u)*np.cos(phi))*exp_electron_f(u))
    
    cota_inf = 0.01
    cota_sup = 600
        
    termEy_pole_aprox_int_re,err = integrate.quad(termEx_num_re_f, cota_inf, cota_sup) 
    termEy_pole_aprox_int_im,err = integrate.quad(termEx_num_im_f, cota_inf, cota_sup) 

    termEy_pole_aprox_int = termEy_pole_aprox_int_re + 1j*termEy_pole_aprox_int_im

    return 1j*0.5*termEy_pole_aprox_int*(k1**3)  


#%%
    
def Ez_polp_pole_aprox(omegac,hbmu,hbgama,R,phi,z,zp,px,py,pz):     
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
    
    alpha_z1 =  lambda u : np.sqrt(1 - u**2) if (1 > np.abs(u)) else 1j*np.sqrt(u**2 -1)   
        
    z_dip_barra = k1*(2*zp - z)
    
    exp_electron_f = lambda u: np.exp(1j*alpha_z1(u)*z_dip_barra)
    J0 =  lambda u : special.jv(0,u*k1*R)  
    J1 =  lambda u : special.jv(1,u*k1*R)

    Rp = 2*epsi1/(epsi1 + epsi2)
    alfa_p = 1j*(epsi1 + epsi2)/(cond)
    rp = lambda u: Rp*alfa_p/(u - alfa_p)   ## u es alfa_parallel y la funcion depende de k_parallel
    
    term_Ex = px*np.cos(phi) + py*np.sin(phi) 
    
    
    termEx_num_re_f = lambda u: np.real(u**2*rp(u)*(J1(u)*term_Ex + 1j*pz*J0(u)*u/alpha_z1(u))*exp_electron_f(u))
    termEx_num_im_f = lambda u: np.imag(u**2*rp(u)*(J1(u)*term_Ex + 1j*pz*J0(u)*u/alpha_z1(u))*exp_electron_f(u))
    
    
    cota_inf = 0.01
    cota_sup = 600
        
    termEz_pole_aprox_int_re,err = integrate.quad(termEx_num_re_f, cota_inf, cota_sup) 
    termEz_pole_aprox_int_im,err = integrate.quad(termEx_num_im_f, cota_inf, cota_sup) 

    termEz_pole_aprox_int = termEz_pole_aprox_int_re + 1j*termEz_pole_aprox_int_im
    return termEz_pole_aprox_int*(k1**3)  


#%%


def Ez_polp_num(omegac,hbmu,hbgama,R,phi,z,zp,px,py,pz):     
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
    
    
    alpha_z1 =  lambda u : np.sqrt(1 - u**2) if (1 > np.abs(u)) else 1j*np.sqrt(u**2 -1)   
    z_dip_barra = k1*(2*zp - z)
    
    exp_electron_f = lambda u: np.exp(1j*alpha_z1(u)*z_dip_barra)
    J0 =  lambda u : special.jv(0,u*k1*R)  
    J1 =  lambda u : special.jv(1,u*k1*R)

    rp_num = lambda u: epsi2*1j*u - epsi1*1j*u - cond*u**2
    rp_den = lambda u: epsi2*1j*u + epsi1*1j*u - cond*u**2
    rp = lambda u: rp_num(u)/rp_den(u) 

    
    term_Ex = px*np.cos(phi) + py*np.sin(phi) 
    
    
    termEx_num_re_f = lambda u: np.real(u**2*rp(u)*(J1(u)*term_Ex + 1j*pz*J0(u)*u/alpha_z1(u))*exp_electron_f(u))
    termEx_num_im_f = lambda u: np.imag(u**2*rp(u)*(J1(u)*term_Ex + 1j*pz*J0(u)*u/alpha_z1(u))*exp_electron_f(u))
    
    cota_inf = 0.01
    cota_sup = 600
        
    termEy_pole_aprox_int_re,err = integrate.quad(termEx_num_re_f, cota_inf, cota_sup,limit = limitt) 
    termEy_pole_aprox_int_im,err = integrate.quad(termEx_num_im_f, cota_inf, cota_sup,limit = limitt) 

    termEy_pole_aprox_int = termEy_pole_aprox_int_re + 1j*termEy_pole_aprox_int_im


    return termEy_pole_aprox_int*(k1**3)  


#%%
    

    
def Ez_polp_ana(omegac,hbmu,hbgama,R,phi,z,zp,px,py,pz):     
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
    k0 = omegac
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)   
    Rp = 2*epsi1/(epsi1 + epsi2)
    alfa_p = 1j*(epsi1 + epsi2)/(cond)
    kp = alfa_p*k0
    
    
    alpha_z1 = np.sqrt(1 - alfa_p**2) if (1 > np.abs(alfa_p)) else 1j*np.sqrt(alfa_p**2 -1)   
        
    z_dip_barra = 2*zp - z
        
    exp_electron_f = np.exp(1j*alpha_z1*k0*z_dip_barra)
    H0 =  special.hankel1(0,kp*R)  
    H1 = special.hankel1(1,kp*R)  

    rp = Rp*alfa_p  ## qj_polo*Rp
    
    term_Ex = px*np.cos(phi) + py*np.sin(phi)   
        
    termEx = alfa_p**2*rp*(H1*term_Ex + 1j*pz*H0*alfa_p/alpha_z1)*exp_electron_f
        
    return 1j*np.pi*termEx*(k0**3)  
 