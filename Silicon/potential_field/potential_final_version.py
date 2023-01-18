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
path_ctes =  path_basic.replace('/' + 'potential_field','')
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

#%%

def function_poles(type_sol,value,d):
    
    os.chdir(path_load)
    tabla_TM1A = np.loadtxt('sol_%s_d%inm.txt' %(type_sol,d*1e3), delimiter='\t', skiprows=1)
    tabla_TM1A = np.transpose(tabla_TM1A)
    [listx,listy] = tabla_TM1A
    
    if np.min(listx) <= value <= np.max(listx):
    
        return float(interp1d(listx, listy)(value))


def coef_fresnel_TM_pole_aprox(omega_omegaWG,d,k_parallel):
    eta = epsilon2
    omegaWG = (np.pi/d)/(np.sqrt(epsilon2-1))  # ya divido por c 
    k0 = omega_omegaWG*omegaWG
    q = k_parallel/k0

    list_all_modes = ['TM1A','TM1B','TM2A','TM2B']
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
        term2 =  (eta*chi)**2 + chi_prima**2/(eta*chi_prima)
        
        Rj_term_final = 1/(2*term1 + term2)
        
        Rj = Rj_aux*Rj_term_final
                
        qj_polo = polo_j/k0
        
        if q - qj_polo != 0:
        
            term0 = term0 + Rj*qj_polo/(q - qj_polo)  
        
    return term0
   

def coef_fresnel_TE_pole_aprox(omega_omegaWG,d,k_parallel):
    eta = epsilon2
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
        term2 =  (eta*chi)**2 + chi_prima**2/(eta*chi_prima)
        
        Rj_term_final = 1/(2*term1 + term2)
        
        Rj = Rj_aux*Rj_term_final
                
        qj_polo = polo_j/k0
        
        if q - qj_polo != 0:
        
            term0 = term0 + Rj*qj_polo/(q - qj_polo)  
        
    return term0

#%%

def G1_num(omega_omegaWG,d,R,z,zp):     
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
    omegaWG = (np.pi/d)/(np.sqrt(epsilon2-1)) 
    omegac = omega_omegaWG*omegaWG
    k1 = omegac
    
    
    z_dip_barra = k1*(2*zp - z)
    
    exp_electron_f = lambda u: np.exp(-u*k1*z_dip_barra)
    J0 =  lambda u : special.jv(0,u*k1*R)

    
    kz1 = lambda u: np.sqrt(omegac**2 - u**2) if (omegac > u) else 1j*np.sqrt(u**2 - omegac**2)
    kz2 = lambda u: np.sqrt(epsilon2*omegac**2 - u**2) if (omegac*np.sqrt(epsilon2) > u) else 1j*np.sqrt(u**2 - epsilon2*omegac**2)
    
    r12_s = lambda u: (kz1(u) - kz2(u))/(kz1(u) + kz2(u))
    r12_p = lambda u: (kz1(u)*epsilon2 - kz2(u))/(kz1(u)*epsilon2 + kz2(u))

    r21_s = lambda u: (kz2(u) - kz1(u))/(kz1(u) + kz2(u))
    r21_p = lambda u: (kz2(u) - kz1(u)*epsilon2)/(kz1(u)*epsilon2 + kz2(u))


    r23_s = lambda u: (kz2(u) - kz1(u))/(kz2(u) + kz1(u))
    r23_p = lambda u: (kz2(u) - kz1(u)*epsilon2)/(kz2(u) + kz1(u)*epsilon2)

    t12_s = lambda u: 2*kz1(u)/(kz1(u) + kz2(u))
    t12_p = lambda u: 2*kz1(u)*np.sqrt(epsilon2)/(kz1(u) + kz2(u))

    t21_s = lambda u: 2*kz2(u)/(kz1(u) + kz2(u))
    t21_p = lambda u: 2*kz2(u)*np.sqrt(epsilon2)/(kz1(u) + kz2(u))

    expo_coef = lambda u: np.exp(2*1j*kz2(u)*d)
        
    rs = lambda u: r12_s(u) + t12_s(u)*t21_s(u)*r23_s(u)*expo_coef(u)/(1 - r21_s(u)*r23_s(u)*expo_coef(u))
    
    rp = lambda u: r12_p(u) + t12_p(u)*t21_p(u)*r23_p(u)*expo_coef(u)/(1 - r21_p(u)*r23_p(u)*expo_coef(u))
    
    termG1_num_re_f = lambda u: np.real(J0(u)*u*(rp(u*k1) + rs(u*k1))*exp_electron_f(u))
    termG1_num_im_f = lambda u: np.imag(J0(u)*u*(rp(u*k1) + rs(u*k1))*exp_electron_f(u))
  
    termG1_num_int0 = 0
    
    cota_alpha_y = 600*k1
    nmax = int(cota_alpha_y*R*omegac/np.pi - 0.5)
    for n in range(0,nmax+1): 

        cota_sup = (2*n+1)*np.pi*0.5/(omegac*R)
        if n == 0:
            cota_inf = 0
        else:
            cota_inf = (2*(n-1)+1)*np.pi*0.5/(omegac*R)
        
        termG1_num_int_re,err = integrate.quad(termG1_num_re_f, cota_inf, cota_sup) 
        termG1_num_int_im,err = integrate.quad(termG1_num_im_f, cota_inf, cota_sup) 
    
        termG1_num_int = termG1_num_int_re + 1j*termG1_num_int_im
    
        termG1_num_int0 = termG1_num_int0 + termG1_num_int
    
    return termG1_num_int0*(k1**2)  

#%%

def G1_pole_aprox(omega_omegaWG,d,R,z,zp):     
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
    
    cota_inf = 1
    cota_sup = 600*k1
    
    z_dip_barra = k1*(2*zp - z)
    
    exp_electron_f = lambda u: np.exp(-u*z_dip_barra)
    J0 = lambda u : special.jv(0,u*k1*R)

    rp = lambda u: coef_fresnel_TM_pole_aprox(omega_omegaWG,d,u*k1)
    rs = lambda u: coef_fresnel_TE_pole_aprox(omega_omegaWG,d,u*k1)
    
    
    termG1_pole_re_f = lambda u: np.real(J0(u)*u*(rp(u) + rs(u))*exp_electron_f(u))
    termG1_pole_im_f = lambda u: np.imag(J0(u)*u*(rp(u) + rs(u))*exp_electron_f(u))
    
    termG1_pole_int_re,err = integrate.quad(termG1_pole_re_f, cota_inf, cota_sup) 
    termG1_pole_int_im,err = integrate.quad(termG1_pole_im_f, cota_inf, cota_sup) 

    termG1_pole_int = termG1_pole_int_re + 1j*termG1_pole_int_im
    
    return termG1_pole_int*(k1**2)  



#%%
    

def G1_ana(omega_omegaWG,d,R,z,zp):     
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
    k0 = omegac
   

    
    z_dip_barra = k1*(-z + 2*zp)        
    exp_electron = np.exp(-alfa_p*z_dip_barra)
     

    H0 = special.hankel1(0,kp*R)
    
    term1 = Rp*np.pi*1j*kp_2*H0*exp_electron


    def coef_fresnel_pole_aprox_ana(omega_omegaWG,d,eta):
        omegaWG = (np.pi*c/d)/(np.sqrt(epsilon2-1)) 

        list_all_modes = ['TM1A','TM1B','TM2A','TM2B','TE1A','TE1B','TE2A','TE2B']
    
        list_polos = []
        for type_of_mode in list_all_modes:  
            rta = function_poles(type_of_mode,omega_omegaWG,d)
            if rta != None:
                list_polos.append(rta)
    
#
#        q = k_parallel/omegac            

        
        term0 = 0
        for j in range(len(list_polos)):
            polo_j = list_polos[j]
            
            
            
            if polo_j**2 <= epsilon2*(k0**2):    
                chi_prima = d*np.sqrt(epsilon2*k0**2 - polo_j**2)
            else:
                chi_prima = d*1j*np.sqrt(polo_j**2 - epsilon2*k0**2)
        
            if polo_j**2 <= k0**2 :
                chi = d*1j*np.sqrt(k0**2 - polo_j**2)
            else:
                chi = d*np.sqrt(polo_j**2 - k0**2)            
                
#            print(polo_j)
            
            den =  (polo_j*d)**2
            
            Rj_aux = 2*chi_prima*chi/den
            
            term1 = (chi**2 + chi_prima**2)/(chi*chi_prima)
            term2 =  (eta*chi)**2 + chi_prima**2/(eta*chi_prima)
            
            Rj_term_final = 1/(2*term1 + term2)
            
            Rj = Rj_aux*Rj_term_final
                    
            qj_polo = polo_j/k0
            
            term0 = term0 + Rj*qj_polo 
            
        return term0
    
###################################################################################################

    return term1/k1_2 


#%%
    
def G2_ana(omega_omegaWG,d,R,z,zp):     
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

    omegaWG = (np.pi*c/d)/(np.sqrt(epsilon2-1)) 
    omega = omega_omegaWG*omegaWG
    omegac = omega/c
   
 
    
    z_dip_barra = k1*(-z + 2*zp)        
    exp_electron = np.exp(-alfa_p*z_dip_barra)
     

    H1 = special.hankel1(1,kp*R)
    
    term1 = Rp*np.pi*1j*kp_2*H1*exp_electron

    
###################################################################################################

    return term1/k1_2 

#%%

def G2_num(omega_omegaWG,d,R,z,zp):     
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
    

    z_dip_barra = 2*zp - z
    
    exp_electron_f = lambda u: np.exp(-u*k1*z_dip_barra)
    J1 =  lambda u : special.jv(1,u*k1*R)

    kz1 = lambda u: np.sqrt(omegac**2 - u**2) if (omegac > u) else 1j*np.sqrt(u**2 - omegac**2)
    kz2 = lambda u: np.sqrt(epsilon2*omegac**2 - u**2) if (omegac*np.sqrt(epsilon2) > u) else 1j*np.sqrt(u**2 - epsilon2*omegac**2)
    
    
    r12_s = lambda u: (kz1(u) - kz2(u))/(kz1(u) + kz2(u))
    r12_p = lambda u: (kz1(u)*epsilon2 - kz2(u))/(kz1(u)*epsilon2 + kz2(u))

    r21_s = lambda u: (kz2(u) - kz1(u))/(kz1(u) + kz2(u))
    r21_p = lambda u: (kz2(u) - kz1(u)*epsilon2)/(kz1(u)*epsilon2 + kz2(u))


    r23_s = lambda u: (kz2(u) - kz1(u))/(kz2(u) + kz1(u))
    r23_p = lambda u: (kz2(u) - kz1(u)*epsilon2)/(kz2(u) + kz1(u)*epsilon2)

    t12_s = lambda u: 2*kz1(u)/(kz1(u) + kz2(u))
    t12_p = lambda u: 2*kz1(u)*np.sqrt(epsilon2)/(kz1(u) + kz2(u))

    t21_s = lambda u: 2*kz2(u)/(kz1(u) + kz2(u))
    t21_p = lambda u: 2*kz2(u)*np.sqrt(epsilon2)/(kz1(u) + kz2(u))

    expo_coef = lambda u: np.exp(2*1j*kz2(u)*d)
        
    rs = lambda u: r12_s(u) + t12_s(u)*t21_s(u)*r23_s(u)*expo_coef(u)/(1 - r21_s(u)*r23_s(u)*expo_coef(u))
    
    rp = lambda u: r12_p(u) + t12_p(u)*t21_p(u)*r23_p(u)*expo_coef(u)/(1 - r21_p(u)*r23_p(u)*expo_coef(u))
    
    termG2_num_re_f = lambda u: np.real(J1(u)*u*(rp(u*k1) + rs(u*k1))*exp_electron_f(u))
    termG2_num_im_f = lambda u: np.imag(J1(u)*u*(rp(u*k1) + rs(u*k1))*exp_electron_f(u))
  
    termG2_num_int0 = 0
    
    cota_alpha_y = 600*k1
    nmax = int(cota_alpha_y*R*omegac/np.pi - 0.5)
    for n in range(0,nmax+1): 

        cota_sup = (2*n+1)*np.pi*0.5/(omegac*R)
        if n == 0:
            cota_inf = 0
        else:
            cota_inf = (2*(n-1)+1)*np.pi*0.5/(omegac*R)
        
        termG2_num_int_re,err = integrate.quad(termG2_num_re_f, cota_inf, cota_sup) 
        termG2_num_int_im,err = integrate.quad(termG2_num_im_f, cota_inf, cota_sup) 
    
        termG2_num_int = termG2_num_int_re + 1j*termG2_num_int_im
    
        termG2_num_int0 = termG2_num_int0 + termG2_num_int
    
    return termG2_num_int0*(k1**2)  

#%%

def G2_pole_aprox(omega_omegaWG,d,R,z,zp):     
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
    
    cota_inf = 1
    cota_sup = 600*k1
    
    z_dip_barra = 2*zp - z
    
    exp_electron_f = lambda u: np.exp(-u*k1*z_dip_barra)
    J1 = lambda u : special.jv(1,u*k1*R)

    rp = lambda u: coef_fresnel_TM_pole_aprox(omega_omegaWG,d,u*k1)
    rs = lambda u: coef_fresnel_TE_pole_aprox(omega_omegaWG,d,u*k1)
    
    
    termG2_pole_re_f = lambda u: np.real(J1(u)*u*(rp(u) + rs(u))*exp_electron_f(u))
    termG2_pole_im_f = lambda u: np.imag(J1(u)*u*(rp(u) + rs(u))*exp_electron_f(u))
    
    termG2_pole_int_re,err = integrate.quad(termG2_pole_re_f, cota_inf, cota_sup) 
    termG2_pole_int_im,err = integrate.quad(termG2_pole_im_f, cota_inf, cota_sup) 

    termG2_pole_int = termG2_pole_int_re + 1j*termG2_pole_int_im
    
    return termG2_pole_int*(k1**2)  


#%%
    

def potential_final_ana(omega_omegaWG,d,phi,R,z,zp,px,py,pz):
    
    G1 =  G1_ana(omega_omegaWG,d,R,z,zp)
    
    G2 =  G2_ana(omega_omegaWG,d,R,z,zp)
    
    
    term0_3 = (np.abs(z)**2 + R**2)**(3/2)
    term0_1 = (np.abs(z)**2 + R**2)**(1/2)
    
    px_py_term = px*np.cos(phi) + py*np.sin(phi)
    
    term1 = -px_py_term*(G2 + (np.abs(z)**2/term0_3 - 1/term0_1)/R)
    
    term2 = -pz*np.sign(z)*(G1 + np.abs(z)/term0_3)
    
    return term1 + term2


#%%s


def potential_final_pole_aprox(omega_omegaWG,d,phi,R,z,zp,px,py,pz):
    
    G1 =  G1_pole_aprox(omega_omegaWG,d,R,z,zp)
    
    G2 =  G2_pole_aprox(omega_omegaWG,d,R,z,zp)
    
    
    term0_3 = (np.abs(z)**2 + R**2)**(3/2)
    term0_1 = (np.abs(z)**2 + R**2)**(1/2)
    
    px_py_term = px*np.cos(phi) + py*np.sin(phi)
    
    term1 = -px_py_term*(G2 + (np.abs(z)**2/term0_3 - 1/term0_1)/R)
    
    term2 = -pz*np.sign(z)*(G1 + np.abs(z)/term0_3)
    
    return term1 + term2

 #%%s   
 
def potential_final_num(omega_omegaWG,d,phi,R,z,zp,px,py,pz):
    
    G1 =  G1_num(omega_omegaWG,d,R,z,zp)
    
    G2 =  G2_num(omega_omegaWG,d,R,z,zp)
    
    term0_3 = (np.abs(z)**2 + R**2)**(3/2)
    term0_1 = (np.abs(z)**2 + R**2)**(1/2)
    
    px_py_term = px*np.cos(phi) + py*np.sin(phi)
    
    term1 = -px_py_term*(G2 + (np.abs(z)**2/term0_3 - 1/term0_1)/R)
    
    term2 = -pz*np.sign(z)*(G1 + np.abs(z)/term0_3)
    
    return term1 + term2 
 
#%%