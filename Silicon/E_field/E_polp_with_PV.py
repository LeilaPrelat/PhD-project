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

cota_inf_alpha_parallel = 1.001
cota_sup_alpha_parallel = np.sqrt(epsilon2) - 0.001

#%%

def function_poles(type_sol,value,d):
    
    os.chdir(path_load)
    tabla_TM1A = np.loadtxt('sol_%s_d%inm.txt' %(type_sol,d*1e3), delimiter='\t', skiprows=1)
    tabla_TM1A = np.transpose(tabla_TM1A)
    [omega_omegaWG,k_parallel] = tabla_TM1A
    
    if np.min(omega_omegaWG) <= value <= np.max(omega_omegaWG):
    
        return float(interp1d(omega_omegaWG, k_parallel)(value))   
    
def Ex_PV_polp_pole_aprox(omega_omegaWG,d,R,phi,z,zp,px,py,pz,epsilon_pole):     
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
    eta = epsilon2
    omegaWG = (np.pi/d)/(np.sqrt(epsilon2-1))  ##saco el c porque divido por c en omegac 
    omegac = omega_omegaWG*omegaWG
    k1 = omegac
    
    
    alpha_z1 = lambda u: 1j*np.sqrt(u**2 - 1)
        
    z_dip_barra = k1*(2*zp - z)
    
    exp_electron_f = lambda u: np.exp(1j*alpha_z1(u)*z_dip_barra)
    J0 =  lambda u : special.jv(0,u*k1*R)  
    J1 =  lambda u : special.jv(1,u*k1*R)  
    J2 =  lambda u : special.jv(2,u*k1*R)

    
    term_Ex_px = lambda u: alpha_z1(u)*(-J0(u) + J2(u)*np.cos(2*phi))
    term_Ex_py = lambda u: alpha_z1(u)*J2(u)*np.sin(2*phi)



    list_all_modes = ['TM1A','TM1B','TM2A','TM2B']
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
        qj_polo  = polo_j/k1
        
        den =  (polo_j*d)**2

        chi_prima =  d*np.sqrt(epsilon2 - qj_polo**2)*k1
        chi =   d*np.sqrt(qj_polo**2 - 1)*k1
    
        
        Rj_aux  =  2*chi_prima*chi/den
        
        term1  = (chi**2 + chi_prima**2)/(chi*chi_prima)
        term2  = ((eta*chi)**2 + chi_prima**2)/(eta*chi_prima)
        
        Rj_term_final = 1/(2*term1 + term2)
        
        Rj  = Rj_aux*Rj_term_final
        
        rp = lambda u: Rj*qj_polo/(u - qj_polo - 1j*epsilon_pole)
    
    
        termEx_num_re_f = lambda u: np.real(u*(px*term_Ex_px(u) + py*term_Ex_py(u) + 2*1j*pz*u*J1(u)*np.cos(phi))*exp_electron_f(u))
        termEx_num_im_f = lambda u: np.imag(u*(px*term_Ex_px(u) + py*term_Ex_py(u) + 2*1j*pz*u*J1(u)*np.cos(phi))*exp_electron_f(u))


        termEx_num_re_f0 = np.real(qj_polo*(px*term_Ex_px(qj_polo) + py*term_Ex_py(qj_polo) + 2*1j*pz*qj_polo*J1(qj_polo)*np.cos(phi))*exp_electron_f(qj_polo))
        termEx_num_im_f0 = np.imag(qj_polo*(px*term_Ex_px(qj_polo) + py*term_Ex_py(qj_polo) + 2*1j*pz*qj_polo*J1(qj_polo)*np.cos(phi))*exp_electron_f(qj_polo))
    
        termEx_num_re_f1 =  lambda u: np.real((termEx_num_re_f(u) - termEx_num_re_f0)*rp(u))
        termEx_num_im_f1 =  lambda u: np.imag((termEx_num_im_f(u) - termEx_num_im_f0)*rp(u))

        termEx_pole_aprox_int_re,err = integrate.quad(termEx_num_re_f1, cota_inf_alpha_parallel, cota_sup_alpha_parallel) 
        termEx_pole_aprox_int_im,err = integrate.quad(termEx_num_im_f1, cota_inf_alpha_parallel, cota_sup_alpha_parallel) 

        termEx_pole_aprox_int = termEx_pole_aprox_int_re + 1j*termEx_pole_aprox_int_im

        term0 = term0 + termEx_pole_aprox_int
    
    return term0*(k1**3)  

#%%

def Ey_PV_polp_pole_aprox(omega_omegaWG,d,R,phi,z,zp,px,py,pz,epsilon_pole):     
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
    eta = epsilon2
    omegaWG = (np.pi/d)/(np.sqrt(epsilon2-1))  ##saco el c porque divido por c en omegac 
    omegac = omega_omegaWG*omegaWG
    k1 = omegac
    
    
    alpha_z1 = lambda u: 1j*np.sqrt(u**2 - 1)
        
    z_dip_barra = k1*(2*zp - z)
    
    exp_electron_f = lambda u: np.exp(1j*alpha_z1(u)*z_dip_barra)
    J0 =  lambda u : special.jv(0,u*k1*R)  
    J1 =  lambda u : special.jv(1,u*k1*R)
    J2 =  lambda u : special.jv(2,u*k1*R)

    
    term_Ex_px = lambda u: -alpha_z1(u)*(J0(u) + J2(u)*np.cos(2*phi))
    term_Ex_py = lambda u: alpha_z1(u)*J2(u)*np.sin(2*phi)


    list_all_modes = ['TM1A','TM1B','TM2A','TM2B']
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
        qj_polo  = polo_j/k1        
        den =  (polo_j*d)**2

        chi_prima =  d*np.sqrt(epsilon2 - qj_polo**2)*k1
        
        chi =  d*np.sqrt(qj_polo**2 - 1)*k1
    
        
        Rj_aux  =  2*chi_prima*chi/den
        
        term1  = (chi**2 + chi_prima**2)/(chi*chi_prima)
        term2  =  ((eta*chi)**2 + chi_prima**2)/(eta*chi_prima)
        
        Rj_term_final  =  1/(2*term1 + term2)
        
        Rj  =  Rj_aux*Rj_term_final
        

        rp = lambda u: Rj*qj_polo/(u - qj_polo - 1j*epsilon_pole)
    
        termEx_num_re_f = lambda u: np.real(u*(py*term_Ex_px(u) + px*term_Ex_py(u) + 2*1j*pz*u*J1(u)*np.cos(phi))*exp_electron_f(u))
        termEx_num_im_f = lambda u: np.imag(u*(py*term_Ex_px(u) + px*term_Ex_py(u) + 2*1j*pz*u*J1(u)*np.cos(phi))*exp_electron_f(u))
 

        termEx_num_re_f0 = np.real(qj_polo*(py*term_Ex_px(qj_polo) + px*term_Ex_py(qj_polo) + 2*1j*pz*qj_polo*J1(qj_polo)*np.cos(phi))*exp_electron_f(qj_polo))
        termEx_num_im_f0 = np.imag(qj_polo*(py*term_Ex_px(qj_polo) + px*term_Ex_py(qj_polo) + 2*1j*pz*qj_polo*J1(qj_polo)*np.cos(phi))*exp_electron_f(qj_polo))
    
        termEx_num_re_f1 =  lambda u: np.real((termEx_num_re_f(u) - termEx_num_re_f0)*rp(u))
        termEx_num_im_f1 =  lambda u: np.imag((termEx_num_im_f(u) - termEx_num_im_f0)*rp(u))
        
           
        termEy_pole_aprox_int_re,err = integrate.quad(termEx_num_re_f1, cota_inf_alpha_parallel, cota_sup_alpha_parallel) 
        termEy_pole_aprox_int_im,err = integrate.quad(termEx_num_im_f1, cota_inf_alpha_parallel, cota_sup_alpha_parallel) 
    
        termEy_pole_aprox_int = termEy_pole_aprox_int_re + 1j*termEy_pole_aprox_int_im

        term0 = term0 + termEy_pole_aprox_int

    return term0*(k1**3)  


#%%
    
def Ez_PV_polp_pole_aprox(omega_omegaWG,d,R,phi,z,zp,px,py,pz,epsilon_pole):     
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
    eta = epsilon2
    omegaWG = (np.pi/d)/(np.sqrt(epsilon2-1))  ##saco el c porque divido por c en omegac 
    omegac = omega_omegaWG*omegaWG
    k1 = omegac
    
    
    alpha_z1 = lambda u: 1j*np.sqrt(u**2 - 1)
        
    z_dip_barra = k1*(2*zp - z)
    
    exp_electron_f = lambda u: np.exp(1j*alpha_z1(u)*z_dip_barra)
    J0 =  lambda u : special.jv(0,u*k1*R)  
    J1 =  lambda u : special.jv(1,u*k1*R)

    
    term_Ex = px*np.cos(phi) + py*np.sin(phi) 

    list_all_modes = ['TM1A','TM1B','TM2A','TM2B']
    list_poles = []
    for type_of_mode in list_all_modes:  
        rta = function_poles(type_of_mode,omega_omegaWG,d)
        if rta != None:
            list_poles.append(rta)    

    term0 = 0
    for j in range(len(list_poles)):
        polo_j = list_poles[j]
#            print(polo_j)
        qj_polo  = polo_j/k1        
        den =  (polo_j*d)**2

        chi_prima =  d*np.sqrt(epsilon2 - qj_polo**2)*k1
        chi =  d*np.sqrt(qj_polo**2 - 1)*k1
    
        
        Rj_aux  = 2*chi_prima*chi/den
        
        term1  = (chi**2 + chi_prima**2)/(chi*chi_prima)
        term2  =  ((eta*chi)**2 + chi_prima**2)/(eta*chi_prima)
        
        Rj_term_final  =  1/(2*term1 + term2)
        
        Rj  =  Rj_aux*Rj_term_final
        

        rp = lambda u: Rj*qj_polo/(u - qj_polo - 1j*epsilon_pole)
    
        termEx_num_re_f = lambda u: np.real(u**2*(J1(u)*term_Ex + 1j*pz*J0(u)*u/alpha_z1(u))*exp_electron_f(u))
        termEx_num_im_f = lambda u: np.imag(u**2*(J1(u)*term_Ex + 1j*pz*J0(u)*u/alpha_z1(u))*exp_electron_f(u))
    
        termEx_num_re_f0 = np.real(qj_polo**2*(J1(qj_polo)*term_Ex + 1j*pz*J0(qj_polo)*qj_polo/alpha_z1(qj_polo))*exp_electron_f(qj_polo))
        termEx_num_im_f0 = np.imag(qj_polo**2*(J1(qj_polo)*term_Ex + 1j*pz*J0(qj_polo)*qj_polo/alpha_z1(qj_polo))*exp_electron_f(qj_polo))    
    
        termEx_num_re_f1 =  lambda u: np.real((termEx_num_re_f(u) - termEx_num_re_f0)*rp(u))
        termEx_num_im_f1 =  lambda u: np.imag((termEx_num_im_f(u) - termEx_num_im_f0)*rp(u))

        
        termEy_pole_aprox_int_re,err = integrate.quad(termEx_num_re_f1, cota_inf_alpha_parallel, cota_sup_alpha_parallel) 
        termEy_pole_aprox_int_im,err = integrate.quad(termEx_num_im_f1, cota_inf_alpha_parallel, cota_sup_alpha_parallel) 

        termEy_pole_aprox_int = termEy_pole_aprox_int_re + 1j*termEy_pole_aprox_int_im
        
        term0 = term0 + termEy_pole_aprox_int
        
    return term0*(k1**3)  


#%%

def Ex_polp_num(omega_omegaWG,d,R,phi,z,zp,px,py,pz):     
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
    J1 =  lambda u : special.jv(1,u*k1*R)  
    J2 =  lambda u : special.jv(2,u*k1*R)
    
    
    r12_p = lambda u: (alpha_z1(u)*epsilon2 - alpha_z2(u))/(alpha_z1(u)*epsilon2 + alpha_z2(u))
    r21_p = lambda u: (alpha_z2(u) - alpha_z1(u)*epsilon2)/(alpha_z1(u)*epsilon2 +  alpha_z2(u))

    r23_p = lambda u: (alpha_z2(u) - alpha_z1(u)*epsilon2)/(alpha_z2(u) + alpha_z1(u)*epsilon2)

    t12_p = lambda u: 2*alpha_z1(u)*np.sqrt(epsilon2)/(alpha_z1(u)*epsilon2 + alpha_z2(u))
    t21_p = lambda u: 2*alpha_z1(u)*np.sqrt(epsilon2)/(alpha_z1(u)*epsilon2 + alpha_z2(u))

    expo_coef = lambda u: np.exp(2*1j*alpha_z2(u)*k1*d)

    rp = lambda u: r12_p(u) + t12_p(u)*t21_p(u)*r23_p(u)*expo_coef(u)/(1 - r21_p(u)*r23_p(u)*expo_coef(u))
    
    term_Ex_px = lambda u: alpha_z1(u)*(-J0(u) + J2(u)*np.cos(2*phi))
    term_Ex_py = lambda u: alpha_z1(u)*J2(u)*np.sin(2*phi)
    
    
    termEx_num_re_f = lambda u: np.real(u*rp(u)*(px*term_Ex_px(u) + py*term_Ex_py(u) + 2*1j*pz*u*J1(u)*np.cos(phi))*exp_electron_f(u))
    termEx_num_im_f = lambda u: np.imag(u*rp(u)*(px*term_Ex_px(u) + py*term_Ex_py(u) + 2*1j*pz*u*J1(u)*np.cos(phi))*exp_electron_f(u))
    
    termEx_num_int_re,err = integrate.quad(termEx_num_re_f, cota_inf_alpha_parallel, cota_sup_alpha_parallel) 
    termEx_num_int_im,err = integrate.quad(termEx_num_im_f, cota_inf_alpha_parallel, cota_sup_alpha_parallel) 

    termEx_num_int = termEx_num_int_re + 1j*termEx_num_int_im


    return 1j*0.5*termEx_num_int*(k1**3)  


#%%
    

#def Ex_polp_pole_aprox(omega_omegaWG,d,R,phi,z,zp,px,py,pz):     
#    """    
#    Parameters
#    ----------
#    omegac : omega/c = k0 en 1/micrometros    
#    epsi1 : epsilon del medio de arriba del plano
#    epsi2 : epsilon del medio de abajo del plano
#    hbmu : chemical potential in eV  
#    hbgama : collision frequency in eV
#    z : coordenada z
#    xD : coordenada x del dipolo 
#    yD : coordenada y del dipolo
#    zD : coordenada z del dipolo 
#    zp : posicion del plano (>0)
#    px : coordenada x del dipolo 
#    py : coordenada y del dipolo
#    pz : coordenada z del dipolo
#    Returns
#    -------
#    formula del potencial electric con QE approximation, rp con 
#    aproximacion del polo y con aprox de principal value para las integrales
#    con rp
#    """
#    eta = epsilon2
#    omegaWG = (np.pi/d)/(np.sqrt(epsilon2-1))  ##saco el c porque divido por c en omegac 
#    omegac = omega_omegaWG*omegaWG
#    k1 = omegac
#    
#    
#    alpha_z1 = lambda u: 1j*np.sqrt(u**2 - 1)
#        
#    z_dip_barra = k1*(2*zp - z)
#    
#    exp_electron_f = lambda u: np.exp(1j*alpha_z1(u)*z_dip_barra)
#    J0 =  lambda u : special.jv(0,u*k1*R)  
#    J1 =  lambda u : special.jv(1,u*k1*R)  
#    J2 =  lambda u : special.jv(2,u*k1*R)
#
#    term_Ex_px = lambda u: alpha_z1(u)*(-J0(u) + J2(u)*np.cos(2*phi))
#    term_Ex_py = lambda u: alpha_z1(u)*J2(u)*np.sin(2*phi)
#        
#    list_all_modes = ['TM1A','TM1B','TM2A','TM2B']
#    list_poles = []
#    for type_of_mode in list_all_modes:  
#        rta = function_poles(type_of_mode,omega_omegaWG,d)
#        if rta != None:
#            list_poles.append(rta)
#            
#    term0 = lambda u: 0
#    for j in range(len(list_poles)):
#        polo_j = list_poles[j]
##            print(polo_j)
#        qj_polo  = polo_j/k1
#        den =  (polo_j*d)**2
#
#        chi_prima = d*np.sqrt(epsilon2 - qj_polo**2)*k1
#        chi = d*np.sqrt(qj_polo**2 - 1)*k1
#        
#        Rj_aux  = 2*chi_prima*chi/den
#        
#        term1  = (chi**2 + chi_prima**2)/(chi*chi_prima)
#        term2  =  ((eta*chi)**2 + chi_prima**2)/(eta*chi_prima)
#        
#        Rj_term_final  =  1/(2*term1 + term2)
#        
#        Rj  =  Rj_aux*Rj_term_final
#        
#        rp = lambda u: Rj*qj_polo/(u - qj_polo - 1j*epsilon_pole) 
#        
#        term0 = lambda u: term0(u) + rp(u)
##    rp = lambda u: coef_fresnel_TM_pole_aprox(omega_omegaWG,d,u)  ## u es alfa_parallel y la funcion depende de k_parallel
#
#        
#    termEx_num_re_f = lambda u: np.real(u*term0(u)*(px*term_Ex_px(u) + py*term_Ex_py(u) + 2*1j*pz*u*J1(u)*np.cos(phi))*exp_electron_f(u))
#    termEx_num_im_f = lambda u: np.imag(u*term0(u)*(px*term_Ex_px(u) + py*term_Ex_py(u) + 2*1j*pz*u*J1(u)*np.cos(phi))*exp_electron_f(u))
#
#        
#    termEx_pole_aprox_int_re,err = integrate.quad(termEx_num_re_f, cota_inf_alpha_parallel, cota_sup_alpha_parallel) 
#    termEx_pole_aprox_int_im,err = integrate.quad(termEx_num_im_f, cota_inf_alpha_parallel, cota_sup_alpha_parallel) 
#
#    termEx_pole_aprox_int = termEx_pole_aprox_int_re + 1j*termEx_pole_aprox_int_im
#
#    PV_term = Ex_PV_polp_pole_aprox(omega_omegaWG,d,R,phi,z,zp,px,py,pz)
#    return 1j*0.5*(np.pi*1j*termEx_pole_aprox_int*(k1**3) + PV_term)


#%%
    

def Ex_polp_ana_and_PV(omega_omegaWG,d,R,phi,z,zp,px,py,pz,epsilon_pole):     
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
    eta = epsilon2
    omegaWG = (np.pi/d)/(np.sqrt(epsilon2-1))  ##saco el c porque divido por c en omegac 
    omegac = omega_omegaWG*omegaWG
    k0 = omegac


    list_all_modes = ['TM1A','TM1B','TM2A','TM2B']
    list_poles = []
    for type_of_mode in list_all_modes:  
        rta = function_poles(type_of_mode,omega_omegaWG,d)
        if rta != None:
            list_poles.append(rta)
            
    term0 = 0
    for j in range(len(list_poles)):
        polo_j = list_poles[j]
        qj_polo = polo_j/omegac
        
#        print(polo_j)

        chi_prima = d*np.sqrt(epsilon2*k0**2 - polo_j**2)
        chi = d*np.sqrt(polo_j**2 - k0**2)
            
        den =  (polo_j*d)**2
        
        Rjaux = 2*chi_prima*chi/den
        
        term1 = (chi**2 + chi_prima**2)/(chi*chi_prima)
        term2 =  ((eta*chi)**2 + chi_prima**2)/(eta*chi_prima)
        
        Rjterm_final = 1/(2*term1 + term2)
        
        Rj = Rjaux*Rjterm_final
        
        alpha_z1 = 1j*np.sqrt(qj_polo**2 - 1)
            
        z_dip_barra = 2*zp - z
        
        exp_electron_f = np.exp(1j*alpha_z1*k0*z_dip_barra)
        J0 = special.jv(0,polo_j*R)  
        J1 = special.jv(1,polo_j*R)  
        J2 = special.jv(2,polo_j*R)
    
        rs = Rj*qj_polo

        term_Ex_px = alpha_z1*(-J0 + J2*np.cos(2*phi))
        term_Ex_py = alpha_z1*J2*np.sin(2*phi)        

        termEx = qj_polo*rs*(px*term_Ex_px + py*term_Ex_py + 2*1j*pz*J1*np.cos(phi)*qj_polo)*exp_electron_f
        term0 = term0 + termEx

    PV_term = Ex_PV_polp_pole_aprox(omega_omegaWG,d,R,phi,z,zp,px,py,pz,epsilon_pole)
    return 1j*0.5*(np.pi*1j*term0*(k0**3) + PV_term)
        

#%%
    
def Ey_polp_ana_and_PV(omega_omegaWG,d,R,phi,z,zp,px,py,pz,epsilon_pole):     
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
    eta = epsilon2
    omegaWG = (np.pi/d)/(np.sqrt(epsilon2-1))  ##saco el c porque divido por c en omegac 
    omegac = omega_omegaWG*omegaWG
    k0 = omegac


    list_all_modes = ['TM1A','TM1B','TM2A','TM2B']
    list_poles = []
    for type_of_mode in list_all_modes:  
        rta = function_poles(type_of_mode,omega_omegaWG,d)
        if rta != None:
            list_poles.append(rta)
            
    term0 = 0
    for j in range(len(list_poles)):
        polo_j = list_poles[j]
        qj_polo = polo_j/omegac

        chi_prima = d*np.sqrt(epsilon2*k0**2 - polo_j**2)
        chi = d*np.sqrt(polo_j**2 - k0**2)
            
        den =  (polo_j*d)**2
        
        Rjaux = 2*chi_prima*chi/den
        
        term1 = (chi**2 + chi_prima**2)/(chi*chi_prima)
        term2 =  ((eta*chi)**2 + chi_prima**2)/(eta*chi_prima)
        
        Rjterm_final = 1/(2*term1 + term2)
        
        Rj = Rjaux*Rjterm_final
        
        alpha_z1 = 1j*np.sqrt(qj_polo**2 - 1)
            
        z_dip_barra = 2*zp - z
        
        exp_electron_f = np.exp(1j*alpha_z1*k0*z_dip_barra)
        J0 =  special.jv(0,polo_j*R)  
        J1 = special.jv(1,polo_j*R)  
        J2 =   special.jv(2,polo_j*R)
    
        rp = Rj*qj_polo  ## qj_polo*Rp
        
        term_Ex_px = -alpha_z1*(J0 + J2*np.cos(2*phi))
        term_Ex_py = alpha_z1*J2*np.sin(2*phi)     
        

        termEx = qj_polo*rp*(py*term_Ex_px + px*term_Ex_py + 2*1j*pz*J1*np.cos(phi)*qj_polo)*exp_electron_f
        term0 = term0 + termEx
        
    PV_term = Ey_PV_polp_pole_aprox(omega_omegaWG,d,R,phi,z,zp,px,py,pz,epsilon_pole)
    return 1j*0.5*(np.pi*1j*term0*(k0**3) + PV_term)

 

#%%

def Ey_polp_num(omega_omegaWG,d,R,phi,z,zp,px,py,pz,):     
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
    J1 =  lambda u : special.jv(1,u*k1*R)
    J2 =  lambda u : special.jv(2,u*k1*R)

    r12_p = lambda u: (alpha_z1(u)*epsilon2 - alpha_z2(u))/(alpha_z1(u)*epsilon2 + alpha_z2(u))
    r21_p = lambda u: (alpha_z2(u) - alpha_z1(u)*epsilon2)/(alpha_z1(u)*epsilon2 +  alpha_z2(u))

    r23_p = lambda u: (alpha_z2(u) - alpha_z1(u)*epsilon2)/(alpha_z2(u) + alpha_z1(u)*epsilon2)

    t12_p = lambda u: 2*alpha_z1(u)*np.sqrt(epsilon2)/(alpha_z1(u)*epsilon2 + alpha_z2(u))
    t21_p = lambda u: 2*alpha_z1(u)*np.sqrt(epsilon2)/(alpha_z1(u)*epsilon2 + alpha_z2(u))

    expo_coef = lambda u: np.exp(2*1j*alpha_z2(u)*k1*d)

    rp = lambda u: r12_p(u) + t12_p(u)*t21_p(u)*r23_p(u)*expo_coef(u)/(1 - r21_p(u)*r23_p(u)*expo_coef(u))
    
    
    term_Ex_px = lambda u: -alpha_z1(u)*(J0(u) + J2(u)*np.cos(2*phi))
    term_Ex_py = lambda u: alpha_z1(u)*J2(u)*np.sin(2*phi)
    
    
    termEx_num_re_f = lambda u: np.real(u*rp(u)*(py*term_Ex_px(u) + px*term_Ex_py(u) + 2*1j*pz*u*J1(u)*np.cos(phi))*exp_electron_f(u))
    termEx_num_im_f = lambda u: np.imag(u*rp(u)*(py*term_Ex_px(u) + px*term_Ex_py(u) + 2*1j*pz*u*J1(u)*np.cos(phi))*exp_electron_f(u))
  
            
    termEy_num_int_re,err = integrate.quad(termEx_num_re_f, cota_inf_alpha_parallel, cota_sup_alpha_parallel) 
    termEy_num_int_im,err = integrate.quad(termEx_num_im_f, cota_inf_alpha_parallel, cota_sup_alpha_parallel) 

    termEx_num_int = termEy_num_int_re + 1j*termEy_num_int_im

    return 1j*0.5*termEx_num_int*(k1**3)  

#%%

#def Ey_polp_pole_aprox(omega_omegaWG,d,R,phi,z,zp,px,py,pz):     
#    """    
#    Parameters
#    ----------
#    omegac : omega/c = k0 en 1/micrometros    
#    epsi1 : epsilon del medio de arriba del plano
#    epsi2 : epsilon del medio de abajo del plano
#    hbmu : chemical potential in eV  
#    hbgama : collision frequency in eV
#    z : coordenada z
#    xD : coordenada x del dipolo 
#    yD : coordenada y del dipolo
#    zD : coordenada z del dipolo 
#    zp : posicion del plano (>0)
#    px : coordenada x del dipolo 
#    py : coordenada y del dipolo
#    pz : coordenada z del dipolo
#    Returns
#    -------
#    formula del potencial electric con QE approximation, rp con 
#    aproximacion del polo y con aprox de principal value para las integrales
#    con rp
#    """
#    eta = epsilon2
#    omegaWG = (np.pi/d)/(np.sqrt(epsilon2-1))  ##saco el c porque divido por c en omegac 
#    omegac = omega_omegaWG*omegaWG
#    k1 = omegac
#    
#    
#    alpha_z1 = lambda u: 1j*np.sqrt(u**2 - 1)
#        
#    z_dip_barra = k1*(2*zp - z)
#    
#    exp_electron_f = lambda u: np.exp(1j*alpha_z1(u)*z_dip_barra)
#    J0 =  lambda u : special.jv(0,u*k1*R)  
#    J1 =  lambda u : special.jv(1,u*k1*R)
#    J2 =  lambda u : special.jv(2,u*k1*R)
#
##    rp = lambda u: coef_fresnel_TM_pole_aprox(omega_omegaWG,d,u)  ## u es alfa_parallel y la funcion depende de k_parallel
#    
#    term_Ex_px = lambda u: -alpha_z1(u)*(J0(u) + J2(u)*np.cos(2*phi))
#    term_Ex_py = lambda u: alpha_z1(u)*J2(u)*np.sin(2*phi)
#    
#    list_all_modes = ['TM1A','TM1B','TM2A','TM2B']
#    list_poles = []
#    for type_of_mode in list_all_modes:  
#        rta = function_poles(type_of_mode,omega_omegaWG,d)
#        if rta != None:
#            list_poles.append(rta)
#            
#    term0 = lambda u: 0
#    for j in range(len(list_poles)):
#        polo_j = list_poles[j]
##            print(polo_j)
#        qj_polo  = polo_j/k1
#        
#        den =  (polo_j*d)**2
#
#        chi_prima = d*np.sqrt(epsilon2 - qj_polo**2)*k1
#        chi = d*np.sqrt(qj_polo**2 - 1)*k1
#        
#        Rj_aux  = 2*chi_prima*chi/den
#        
#        term1  = (chi**2 + chi_prima**2)/(chi*chi_prima)
#        term2  =  ((eta*chi)**2 + chi_prima**2)/(eta*chi_prima)
#        
#        Rj_term_final  =  1/(2*term1 + term2)
#        
#        Rj  =  Rj_aux*Rj_term_final
#          
#        rp = lambda u: Rj*qj_polo/(u - qj_polo - 1j*epsilon_pole) 
#        
#        term0 = lambda u: term0(u) + rp(u)
#        
#        
#            
#    termEx_num_re_f = lambda u: np.real(u*term0(u)*(py*term_Ex_px(u) + px*term_Ex_py(u) + 2*1j*pz*u*J1(u)*np.cos(phi))*exp_electron_f(u))
#    termEx_num_im_f = lambda u: np.imag(u*term0(u)*(py*term_Ex_px(u) + px*term_Ex_py(u) + 2*1j*pz*u*J1(u)*np.cos(phi))*exp_electron_f(u))
#
# 
#    termEy_pole_aprox_int_re,err = integrate.quad(termEx_num_re_f, cota_inf_alpha_parallel, cota_sup_alpha_parallel) 
#    termEy_pole_aprox_int_im,err = integrate.quad(termEx_num_im_f, cota_inf_alpha_parallel, cota_sup_alpha_parallel) 
#
#    termEy_pole_aprox_int = termEy_pole_aprox_int_re + 1j*termEy_pole_aprox_int_im
#
#       
#    PV_term = Ey_PV_polp_pole_aprox(omega_omegaWG,d,R,phi,z,zp,px,py,pz)
#    
#    return 1j*0.5*(np.pi*1j*termEy_pole_aprox_int*(k1**3) + PV_term)


#%%
    
#def Ez_polp_pole_aprox(omega_omegaWG,d,R,phi,z,zp,px,py,pz):     
#    """    
#    Parameters
#    ----------
#    omegac : omega/c = k0 en 1/micrometros    
#    epsi1 : epsilon del medio de arriba del plano
#    epsi2 : epsilon del medio de abajo del plano
#    hbmu : chemical potential in eV  
#    hbgama : collision frequency in eV
#    z : coordenada z
#    xD : coordenada x del dipolo 
#    yD : coordenada y del dipolo
#    zD : coordenada z del dipolo 
#    zp : posicion del plano (>0)
#    px : coordenada x del dipolo 
#    py : coordenada y del dipolo
#    pz : coordenada z del dipolo
#    Returns
#    -------
#    formula del potencial electric con QE approximation, rp con 
#    aproximacion del polo y con aprox de principal value para las integrales
#    con rp
#    """
#    eta = epsilon2
#    omegaWG = (np.pi/d)/(np.sqrt(epsilon2-1))  ##saco el c porque divido por c en omegac 
#    omegac = omega_omegaWG*omegaWG
#    k1 = omegac
#    
#    
#    alpha_z1 = lambda u: 1j*np.sqrt(u**2 - 1)
#        
#    z_dip_barra = k1*(2*zp - z)
#    
#    exp_electron_f = lambda u: np.exp(1j*alpha_z1(u)*z_dip_barra)
#    J0 =  lambda u : special.jv(0,u*k1*R)  
#    J1 =  lambda u : special.jv(1,u*k1*R)
##
##    rp = lambda u: coef_fresnel_TM_pole_aprox(omega_omegaWG,d,u)  ## u es alfa_parallel y la funcion depende de k_parallel
#    
#    term_Ex = px*np.cos(phi) + py*np.sin(phi) 
#    
#    list_all_modes = ['TM1A','TM1B','TM2A','TM2B']
#    list_poles = []
#    for type_of_mode in list_all_modes:  
#        rta = function_poles(type_of_mode,omega_omegaWG,d)
#        if rta != None:
#            list_poles.append(rta)
#            
#    term0 = 0
#    for j in range(len(list_poles)):
#        polo_j = list_poles[j]
##            print(polo_j)
#        qj_polo  = polo_j/k1        
#        den =  (polo_j*d)**2
#
#        chi_prima =  d*np.sqrt(epsilon2 - polo_j**2)*k1
#        chi = d*np.sqrt(polo_j**2 - 1)*k1
#        
#        Rj_aux  =  2*chi_prima*chi/den
#        
#        term1  = (chi**2 + chi_prima**2)/(chi*chi_prima)
#        term2  =  ((eta*chi)**2 + chi_prima**2)/(eta*chi_prima)
#        
#        Rj_term_final  = 1/(2*term1 + term2)
#        
#        Rj  = Rj_aux*Rj_term_final
#        
#        rp = lambda u: Rj*qj_polo/(u - qj_polo - 1j*epsilon_pole) 
#            
#        
#        termEx_num_re_f = lambda u: np.real(u**2*rp(u)*(J1(u)*term_Ex + 1j*pz*J0(u)*u/alpha_z1(u))*exp_electron_f(u))
#        termEx_num_im_f = lambda u: np.imag(u**2*rp(u)*(J1(u)*term_Ex + 1j*pz*J0(u)*u/alpha_z1(u))*exp_electron_f(u))
#        
#            
#        termEy_pole_aprox_int_re,err = integrate.quad(termEx_num_re_f, cota_inf_alpha_parallel, cota_sup_alpha_parallel) 
#        termEy_pole_aprox_int_im,err = integrate.quad(termEx_num_im_f, cota_inf_alpha_parallel, cota_sup_alpha_parallel) 
#    
#        termEy_pole_aprox_int = termEy_pole_aprox_int_re + 1j*termEy_pole_aprox_int_im
#        
#        term0 = term0 + termEy_pole_aprox_int
#        
#    PV_term = Ez_PV_polp_pole_aprox(omega_omegaWG,d,R,phi,z,zp,px,py,pz)        
#    return np.pi*1j*term0*(k1**3) + PV_term


#%%


def Ez_polp_num(omega_omegaWG,d,R,phi,z,zp,px,py,pz):     
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
    alpha_z2 = lambda u: np.sqrt(epsilon2 - u**2)        
    z_dip_barra = k1*(2*zp - z)
    
    exp_electron_f = lambda u: np.exp(1j*alpha_z1(u)*z_dip_barra)
    J0 =  lambda u : special.jv(0,u*k1*R)  
    J1 =  lambda u : special.jv(1,u*k1*R)

    r12_p = lambda u: (alpha_z1(u)*epsilon2 - alpha_z2(u))/(alpha_z1(u)*epsilon2 + alpha_z2(u))
    r21_p = lambda u: (alpha_z2(u) - alpha_z1(u)*epsilon2)/(alpha_z1(u)*epsilon2 +  alpha_z2(u))

    r23_p = lambda u: (alpha_z2(u) - alpha_z1(u)*epsilon2)/(alpha_z2(u) + alpha_z1(u)*epsilon2)

    t12_p = lambda u: 2*alpha_z1(u)*np.sqrt(epsilon2)/(alpha_z1(u)*epsilon2 + alpha_z2(u))
    t21_p = lambda u: 2*alpha_z1(u)*np.sqrt(epsilon2)/(alpha_z1(u)*epsilon2 + alpha_z2(u))

    expo_coef = lambda u: np.exp(2*1j*alpha_z2(u)*k1*d)

    rp = lambda u: r12_p(u) + t12_p(u)*t21_p(u)*r23_p(u)*expo_coef(u)/(1 - r21_p(u)*r23_p(u)*expo_coef(u))
    
    term_Ex = px*np.cos(phi) + py*np.sin(phi) 
    
    
    termEx_num_re_f = lambda u: np.real(u**2*rp(u)*(J1(u)*term_Ex + 1j*pz*J0(u)*u/alpha_z1(u))*exp_electron_f(u))
    termEx_num_im_f = lambda u: np.imag(u**2*rp(u)*(J1(u)*term_Ex + 1j*pz*J0(u)*u/alpha_z1(u))*exp_electron_f(u))
    

        
    termEy_pole_aprox_int_re,err = integrate.quad(termEx_num_re_f, cota_inf_alpha_parallel, cota_sup_alpha_parallel) 
    termEy_pole_aprox_int_im,err = integrate.quad(termEx_num_im_f, cota_inf_alpha_parallel, cota_sup_alpha_parallel) 

    termEy_pole_aprox_int = termEy_pole_aprox_int_re + 1j*termEy_pole_aprox_int_im


    return termEy_pole_aprox_int*(k1**3)  


#%%
    

    
def Ez_polp_ana_and_PV(omega_omegaWG,d,R,phi,z,zp,px,py,pz,epsilon_pole):     
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
    eta = epsilon2
    omegaWG = (np.pi/d)/(np.sqrt(epsilon2-1))  ##saco el c porque divido por c en omegac 
    omegac = omega_omegaWG*omegaWG
    k0 = omegac


    list_all_modes = ['TM1A','TM1B','TM2A','TM2B']
    list_poles = []
    for type_of_mode in list_all_modes:  
        rta = function_poles(type_of_mode,omega_omegaWG,d)
        if rta != None:
            list_poles.append(rta)
            
    term0 = 0
    for j in range(len(list_poles)):
        polo_j = list_poles[j]
        qj_polo = polo_j/omegac

        chi_prima = d*np.sqrt(epsilon2*k0**2 - polo_j**2)
        chi = d*np.sqrt(polo_j**2 - k0**2)
            
        den =  (polo_j*d)**2
        
        Rjaux = 2*chi_prima*chi/den
        
        term1 = (chi**2 + chi_prima**2)/(chi*chi_prima)
        term2 =  ((eta*chi)**2 + chi_prima**2)/(eta*chi_prima)
        
        Rjterm_final = 1/(2*term1 + term2)
        
        Rj = Rjaux*Rjterm_final
        
        alpha_z1 = 1j*np.sqrt(qj_polo**2 - 1)
            
        z_dip_barra = 2*zp - z
        
        
        exp_electron_f = np.exp(1j*alpha_z1*k0*z_dip_barra)
        J0 =  special.jv(0,polo_j*R)  
        J1 = special.jv(1,polo_j*R)  
    
        rp = Rj*qj_polo  ## qj_polo*Rp
        
    
        term_Ex = px*np.cos(phi) + py*np.sin(phi)   
        
        
        termEx = qj_polo**2*rp*(J1*term_Ex + 1j*pz*J0*qj_polo/alpha_z1)*exp_electron_f
        term0 = term0 + termEx
 

    PV_term = Ez_PV_polp_pole_aprox(omega_omegaWG,d,R,phi,z,zp,px,py,pz,epsilon_pole)
    
    return np.pi*1j*term0*(k0**3) + PV_term
      
 