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
from scipy.interpolate import interp1d

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

pi,hb,c,alfac,mu1,mu2,mu3 = constantes()
aux = c*hb

epsilon2 = 12 ### los polos fueron hallados para este valor de epsilon2 

cota_inf_alpha_parallel = 1.05
cota_sup_alpha_parallel = np.sqrt(epsilon2) - 0.05


epsilon_pole = 1e-4

#%%

def function_poles(type_sol,value,d):
    
    os.chdir(path_load)
    tabla_TM1A = np.loadtxt('sol_%s_d%inm.txt' %(type_sol,d*1e3), delimiter='\t', skiprows=1)
    tabla_TM1A = np.transpose(tabla_TM1A)
    [listx,listy] = tabla_TM1A
    
    if np.min(listx) <= value <= np.max(listx):
    
        return float(interp1d(listx, listy)(value))   

def coef_fresnel_TE_pole_aprox(omega_omegaWG,d,k_parallel):
    eta = 1
    omegaWG = (np.pi/d)/(np.sqrt(epsilon2-1))  # ya divido por c 
    k0 = omega_omegaWG*omegaWG
    q = k_parallel/k0

    list_all_modes = ['TE1A','TE1B','TE2A','TE2B']
    list_poles = []
    for type_of_mode in list_all_modes:  
        rta = function_poles(type_of_mode,omega_omegaWG,d)
        if rta != None:
            list_poles.append(rta)

#    if number_of_poles != len(list_poles):
#        raise TypeError('TM: N(omega) != [omega/omegaWG]')
        
    if k_parallel**2 <= epsilon2*(k0**2):    
        chi_prima = d*np.sqrt(epsilon2*k0**2 - k_parallel**2)
    else:
        chi_prima = d*1j*np.sqrt(k_parallel**2 - epsilon2*k0**2)

    if k_parallel**2 <= k0**2 :
        chi = d*1j*np.sqrt(k0**2 - k_parallel**2)
    else:
        chi = d*np.sqrt(k_parallel**2 - k0**2)
    
    term0 = 0
    for j in range(len(list_poles)):
        polo_j = list_poles[j]
#            print(polo_j)
        
        den =  (polo_j*d)**2
        
        Rj_aux = 2*chi_prima*chi/den
        
        term1 = (chi**2 + chi_prima**2)/(chi*chi_prima)
        term2 =  ((eta*chi)**2 + chi_prima**2)/(eta*chi_prima)
        
        Rj_term_final = 1/(2*term1 + term2)
        
        Rj = Rj_aux*Rj_term_final
                
        qj_polo = polo_j/k0
        
        if q - qj_polo != 0:
        
            term0 = term0 + Rj*qj_polo/(q - qj_polo - 1j*epsilon_pole)  
        
    return term0


def coef_fresnel_TE_ana(omega_omegaWG,d):
    eta = 1
    omegaWG = (np.pi/d)/(np.sqrt(epsilon2-1))  # ya divido por c 
    k0 = omega_omegaWG*omegaWG

    list_all_modes = ['TE1A','TE1B','TE2A','TE2B']
    list_poles = []
    for type_of_mode in list_all_modes:  
        rta = function_poles(type_of_mode,omega_omegaWG,d)
        if rta != None:
            list_poles.append(rta)

#    if number_of_poles != len(list_poles):
#        raise TypeError('TM: N(omega) != [omega/omegaWG]')
    
    term0 = 0
    for j in range(len(list_poles)):
        polo_j = list_poles[j]
#            print(polo_j)
        
        if polo_j**2 <= epsilon2*(k0**2):    
            chi_prima = d*np.sqrt(epsilon2*k0**2 - polo_j**2)
        else:
            chi_prima = d*1j*np.sqrt(polo_j**2 - epsilon2*k0**2)
    
        if polo_j**2 <= k0**2 :
            chi = d*1j*np.sqrt(k0**2 - polo_j**2)
        else:
            chi = d*np.sqrt(polo_j**2 - k0**2)
        
        den =  (polo_j*d)**2
        
        Rj_aux = 2*chi_prima*chi/den
        
        term1 = (chi**2 + chi_prima**2)/(chi*chi_prima)
        term2 =  ((eta*chi)**2 + chi_prima**2)/(eta*chi_prima)
        
        Rj_term_final = 1/(2*term1 + term2)
        
        Rj = Rj_aux*Rj_term_final
                
        qj_polo = polo_j/k0
        
        term0 = term0 + Rj*qj_polo
        
    return term0

#%%

def Ex_pols_num(omega_omegaWG,d,R,phi,z,zp,px,py,pz):     
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
    omegaWG = (np.pi/d)/(np.sqrt(epsilon2-1)) 
    omegac = omega_omegaWG*omegaWG
    k1 = omegac
    
    
    alpha_z1 = lambda u: 1j*np.sqrt(u**2 - 1)
    alpha_z2 = lambda u: np.sqrt(epsilon2 - u**2)
        
    z_dip_barra = k1*(2*zp - z)
    
    exp_electron_f = lambda u: np.exp(1j*alpha_z1(u)*z_dip_barra)
    J0 =  lambda u : special.jv(0,u*k1*R)  
    J2 =  lambda u : special.jv(2,u*k1*R)

    r12_s = lambda u: (alpha_z1(u) - alpha_z2(u))/(alpha_z1(u) + alpha_z2(u)) 
    r21_s = lambda u: (alpha_z2(u) - alpha_z1(u))/(alpha_z1(u) + alpha_z2(u))

    r23_s = lambda u: (alpha_z2(u) - alpha_z1(u))/(alpha_z2(u) + alpha_z1(u))
    t12_s = lambda u: 2*alpha_z1(u)/(alpha_z1(u) + alpha_z2(u))

    t21_s = lambda u: 2*alpha_z2(u)/(alpha_z1(u) + alpha_z2(u))

    expo_coef = lambda u: np.exp(2*1j*k1*alpha_z2(u)*d)
        
    rs = lambda u: r12_s(u) + t12_s(u)*t21_s(u)*r23_s(u)*expo_coef(u)/(1 - r21_s(u)*r23_s(u)*expo_coef(u))
    
    term_Ex_px = lambda u: J0(u) + J2(u)*np.cos(2*phi)
    term_Ex_py = lambda u: J2(u)*np.sin(2*phi)
    
    
    termEx_num_re_f = lambda u: np.real(u*rs(u)*(px*term_Ex_px(u) + py*term_Ex_py(u))*exp_electron_f(u)/alpha_z1(u))
    termEx_num_im_f = lambda u: np.imag(u*rs(u)*(px*term_Ex_px(u) + py*term_Ex_py(u))*exp_electron_f(u)/alpha_z1(u))
  
    
    termEx_num_int_re,err = integrate.quad(termEx_num_re_f, cota_inf_alpha_parallel, cota_sup_alpha_parallel) 
    termEx_num_int_im,err = integrate.quad(termEx_num_im_f, cota_inf_alpha_parallel, cota_sup_alpha_parallel) 

    termEx_num_int = termEx_num_int_re + 1j*termEx_num_int_im

    return 1j*0.5*termEx_num_int*(k1**3)  

#%%

def Ex_pols_pole_aprox(omega_omegaWG,d,R,phi,z,zp,px,py,pz):     
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

    omegaWG = (np.pi/d)/(np.sqrt(epsilon2-1))  ##saco el c porque divido por c en omegac 
    omegac = omega_omegaWG*omegaWG
    k1 = omegac
    
    
    alpha_z1 = lambda u: 1j*np.sqrt(u**2 - 1)
        
    z_dip_barra = k1*(2*zp - z)
    
    exp_electron_f = lambda u: np.exp(1j*alpha_z1(u)*z_dip_barra)
    J0 =  lambda u : special.jv(0,u*k1*R)  
    J2 =  lambda u : special.jv(2,u*k1*R)

    rs = lambda u: coef_fresnel_TE_pole_aprox(omega_omegaWG,d,u*k1)  ## u es alfa_parallel y la funcion depende de k_parallel
    
    term_Ex_px = lambda u: J0(u) + J2(u)*np.cos(2*phi)
    term_Ex_py = lambda u: J2(u)*np.sin(2*phi)
    
    
    termEx_num_re_f = lambda u: np.real(u*rs(u)*(px*term_Ex_px(u) + py*term_Ex_py(u))*exp_electron_f(u)/alpha_z1(u))
    termEx_num_im_f = lambda u: np.imag(u*rs(u)*(px*term_Ex_px(u) + py*term_Ex_py(u))*exp_electron_f(u)/alpha_z1(u))

    
    termEx_pole_aprox_int_re,err = integrate.quad(termEx_num_re_f, cota_inf_alpha_parallel, cota_sup_alpha_parallel) 
    termEx_pole_aprox_int_im,err = integrate.quad(termEx_num_im_f, cota_inf_alpha_parallel, cota_sup_alpha_parallel) 

    termEx_pole_aprox_int = termEx_pole_aprox_int_re + 1j*termEx_pole_aprox_int_im

    return 1j*0.5*termEx_pole_aprox_int*(k1**3)  


#%%
    

def Ex_pols_ana(omega_omegaWG,d,R,phi,z,zp,px,py,pz):     
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

    omegaWG = (np.pi/d)/(np.sqrt(epsilon2-1))  ##saco el c porque divido por c en omegac 
    omegac = omega_omegaWG*omegaWG
    k1 = omegac


    list_all_modes = ['TE1A','TE1B','TE2A','TE2B']
    list_poles = []
    for type_of_mode in list_all_modes:  
        rta = function_poles(type_of_mode,omega_omegaWG,d)
        if rta != None:
            list_poles.append(rta)
            
    term0 = 0
    for j in range(len(list_poles)):
        polo_j = list_poles[j]
        qj_polo = polo_j/omegac
        
        alpha_z1 = 1j*np.sqrt(qj_polo**2 - 1)
            
        z_dip_barra = 2*zp - z
        
        exp_electron_f = np.exp(1j*alpha_z1*k1*z_dip_barra)  ## simetria
        H0 =  special.hankel1(0,polo_j*R)  
        H2 =   special.hankel1(2,polo_j*R)
    
        rs = np.real(coef_fresnel_TE_ana(omega_omegaWG,d))  ## qj_polo*Rp
        
        term_Ex_px = H0 + H2*np.cos(2*phi)
        term_Ex_py = H2*np.sin(2*phi)
        

        termEx = qj_polo*rs*(px*term_Ex_px + py*term_Ex_py)*exp_electron_f/alpha_z1
        term0 = term0 + termEx
        
    return -np.pi*0.5*term0*(k1**3)  
 

#%%
    
def Ey_pols_ana(omega_omegaWG,d,R,phi,z,zp,px,py,pz):     
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

    omegaWG = (np.pi/d)/(np.sqrt(epsilon2-1))  ##saco el c porque divido por c en omegac 
    omegac = omega_omegaWG*omegaWG
    k1 = omegac


    list_all_modes = ['TE1A','TE1B','TE2A','TE2B']
    list_poles = []
    for type_of_mode in list_all_modes:  
        rta = function_poles(type_of_mode,omega_omegaWG,d)
        if rta != None:
            list_poles.append(rta)
            
    term0 = 0
    for j in range(len(list_poles)):
        polo_j = list_poles[j]
        qj_polo = polo_j/omegac
        
        alpha_z1 = 1j*np.sqrt(qj_polo**2 - 1)
            
        z_dip_barra = 2*zp - z
        
        exp_electron_f = np.exp(1j*alpha_z1*k1*z_dip_barra)  ## simetria
        H0 =  special.hankel1(0,polo_j*R)  
        H2 =   special.hankel1(2,polo_j*R)
    
        rs = np.real(coef_fresnel_TE_ana(omega_omegaWG,d))  ## qj_polo*Rp
        
        term_Ex_py = H0 - H2*np.cos(2*phi)
        term_Ex_px = H2*np.sin(2*phi)
        

        termEx = qj_polo*rs*(px*term_Ex_px + py*term_Ex_py)*exp_electron_f/alpha_z1
        term0 = term0 + termEx
        
    return -np.pi*0.5*term0*(k1**3)  
 

#%%

def Ey_pols_num(omega_omegaWG,d,R,phi,z,zp,px,py,pz):     
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
    omegaWG = (np.pi/d)/(np.sqrt(epsilon2-1)) 
    omegac = omega_omegaWG*omegaWG
    k1 = omegac
    
    
    alpha_z1 = lambda u: 1j*np.sqrt(u**2 - 1)
    alpha_z2 = lambda u: np.sqrt(epsilon2 - u**2)
        
    z_dip_barra = k1*(2*zp - z)
    
    exp_electron_f = lambda u: np.exp(1j*alpha_z1(u)*z_dip_barra)
    J0 =  lambda u : special.jv(0,u*k1*R)  
    J2 =  lambda u : special.jv(2,u*k1*R)

    r12_s = lambda u: (alpha_z1(u) - alpha_z2(u))/(alpha_z1(u) + alpha_z2(u)) 
    r21_s = lambda u: (alpha_z2(u) - alpha_z1(u))/(alpha_z1(u) + alpha_z2(u))

    r23_s = lambda u: (alpha_z2(u) - alpha_z1(u))/(alpha_z2(u) + alpha_z1(u))
    t12_s = lambda u: 2*alpha_z1(u)/(alpha_z1(u) + alpha_z2(u))

    t21_s = lambda u: 2*alpha_z2(u)/(alpha_z1(u) + alpha_z2(u))

    expo_coef = lambda u: np.exp(2*1j*k1*alpha_z2(u)*d)
        
    rs = lambda u: r12_s(u) + t12_s(u)*t21_s(u)*r23_s(u)*expo_coef(u)/(1 - r21_s(u)*r23_s(u)*expo_coef(u))
    
    term_Ex_py = lambda u: J0(u) - J2(u)*np.cos(2*phi)
    term_Ex_px = lambda u: J2(u)*np.sin(2*phi)
    
    
    termEx_num_re_f = lambda u: np.real(u*rs(u)*(px*term_Ex_px(u) + py*term_Ex_py(u))*exp_electron_f(u)/alpha_z1(u))
    termEx_num_im_f = lambda u: np.imag(u*rs(u)*(px*term_Ex_px(u) + py*term_Ex_py(u))*exp_electron_f(u)/alpha_z1(u))
  

    
    termEy_num_int_re,err = integrate.quad(termEx_num_re_f, cota_inf_alpha_parallel, cota_sup_alpha_parallel) 
    termEy_num_int_im,err = integrate.quad(termEx_num_im_f, cota_inf_alpha_parallel, cota_sup_alpha_parallel) 

    termEx_num_int = termEy_num_int_re + 1j*termEy_num_int_im

    return 1j*0.5*termEx_num_int*(k1**3)  

#%%

def Ey_pols_pole_aprox(omega_omegaWG,d,R,phi,z,zp,px,py,pz):     
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
    omegaWG = (np.pi/d)/(np.sqrt(epsilon2-1))  ##saco el c porque divido por c en omegac 

    omegaWG = (np.pi/d)/(np.sqrt(epsilon2-1))  ##saco el c porque divido por c en omegac 
    omegac = omega_omegaWG*omegaWG
    k1 = omegac
    
    
    alpha_z1 = lambda u: 1j*np.sqrt(u**2 - 1)
        
    z_dip_barra = k1*(2*zp - z)
    
    exp_electron_f = lambda u: np.exp(1j*alpha_z1(u)*z_dip_barra)
    J0 =  lambda u : special.jv(0,u*k1*R)  
    J2 =  lambda u : special.jv(2,u*k1*R)

    rs = lambda u: coef_fresnel_TE_pole_aprox(omega_omegaWG,d,u*k1)  ## u es alfa_parallel y la funcion depende de k_parallel
    
    term_Ex_py = lambda u: J0(u) - J2(u)*np.cos(2*phi)
    term_Ex_px = lambda u: J2(u)*np.sin(2*phi)
    
    
    termEx_num_re_f = lambda u: np.real(u*rs(u)*(px*term_Ex_px(u) + py*term_Ex_py(u))*exp_electron_f(u)/alpha_z1(u))
    termEx_num_im_f = lambda u: np.imag(u*rs(u)*(px*term_Ex_px(u) + py*term_Ex_py(u))*exp_electron_f(u)/alpha_z1(u))
  

    termEy_pole_aprox_int_re,err = integrate.quad(termEx_num_re_f, cota_inf_alpha_parallel, cota_sup_alpha_parallel) 
    termEy_pole_aprox_int_im,err = integrate.quad(termEx_num_im_f, cota_inf_alpha_parallel, cota_sup_alpha_parallel) 

    termEy_pole_aprox_int = termEy_pole_aprox_int_re + 1j*termEy_pole_aprox_int_im

    return 1j*0.5*termEy_pole_aprox_int*(k1**3)   


#%%
    

