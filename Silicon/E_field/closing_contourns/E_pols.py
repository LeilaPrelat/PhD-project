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
path_ctes =  path_basic.replace('/' + 'E_field/closing_contourns','')
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

#cota_inf_alpha_parallel = 1.001
#cota_sup_alpha_parallel = np.sqrt(epsilon2) - 0.001 

epsilon_pole = 1e-10 ### ayuda a la convergencia. parece que con 1e-3 funciona sin problema
limitt = 70


#%%

def function_poles(type_sol,value,d):
    
    os.chdir(path_load)
    tabla_TM1A = np.loadtxt('sol_%s_d%inm.txt' %(type_sol,d*1e3), delimiter='\t', skiprows=1)
    tabla_TM1A = np.transpose(tabla_TM1A)
    [omega_omegaWG,k_parallel] = tabla_TM1A
    
    if np.min(omega_omegaWG) <= value <= np.max(omega_omegaWG):
    
        return float(interp1d(omega_omegaWG, k_parallel)(value))   
    

def integration_limits(omega_omegaWG,d):

    def k_parallel_air(omega):
    
        return omega/c
    
    def k_parallel_medium(omega,epsilon2):
    
        return omega*np.sqrt(epsilon2)/c

    omegaWG = (np.pi*c/d)/(np.sqrt(epsilon2 - 1)) 
    omega = omega_omegaWG*omegaWG
    k1 = omega/c
    
#    kair =  k_parallel_air(omega)
#    kmedium = k_parallel_medium(omega,epsilon2)

    list_all_modes = ['TE1A','TE1B','TE2A','TE2B']
    list_poles = []
    for type_of_mode in list_all_modes:  
        rta = function_poles(type_of_mode,omega_omegaWG,d)
        if rta != None:
            list_poles.append(rta)
    
    list_poles = sorted(list_poles)
    
    inf_limits = []
    sup_limits = []
    
    for pole in list_poles:
        deltta = 0.2*k1
        if pole - deltta > k1:
            inf_limits.append(pole - deltta)
        else:
            inf_limits.append(k1 + 0.01) ## el integrando diverge en kz = 0 (k parallel = k)
        

        if pole + deltta < np.sqrt(epsilon2)*k1:
            sup_limits.append(pole + deltta)
        else:
            sup_limits.append(np.sqrt(epsilon2)*k1)
 
#        print(pole/k1, inf_limits/k1, sup_limits/k1)
    
    return [inf_limits,sup_limits]

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
    
    
    alpha_z1 = lambda u: np.sqrt(1 - u**2) if (1 > np.abs(u)) else 1j*np.sqrt(u**2 -1)   
    alpha_z2 = lambda u: np.sqrt(epsilon2 - u**2) if (np.sqrt(epsilon2) > np.abs(u)) else 1j*np.sqrt(u**2 - epsilon2)   
        
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
  
    term0_re = 0
    
    [inf_limits,sup_limits] = integration_limits(omega_omegaWG,d)
    for j in range(len(sup_limits)):
    
        
        inf_limits_int = inf_limits[j]/k1
        sup_limits_int = sup_limits[j]/k1
        
        
        termEx_num_int_re,err = integrate.quad(termEx_num_re_f, inf_limits_int, sup_limits_int, limit=limitt) 
        termEx_num_int_im,err = integrate.quad(termEx_num_im_f, inf_limits_int, sup_limits_int, limit=limitt) 

        termEx_num_int = termEx_num_int_re + 1j*termEx_num_int_im
        
        term0_re = term0_re + termEx_num_int
#


#    sin_electron_f = lambda u: np.sin(alpha_z1(u)*z_dip_barra)    
#    termEx_num_re_f2 = lambda u: -u*np.imag(rs(u))*(px*term_Ex_px(u) + py*term_Ex_py(u))*sin_electron_f(u)/alpha_z1(u)
#    
    
    
    cota_inf_div =  0.001
    cota_sup_div = 0.999
    listy_re = []
    listx = np.linspace(cota_inf_div,cota_sup_div,3000)
    for u in listx:
        
        listy_re.append(termEx_num_re_f(u))
    
    zeros_re = []
    for j in range(len(listy_re)-1):
        a = listy_re[j]
        b = listy_re[j+1]
        
        if a*b < 0 :
            
            value_zero = np.mean([listx[j],listx[j+1]])
            zeros_re.append(value_zero)



    listy_im = []
    for u in listx:
        
        listy_im.append(termEx_num_im_f(u))
    
    zeros_im = []
    for j in range(len(listy_im)-1):
        a = listy_im[j]
        b = listy_im[j+1]
        
        if a*b < 0 :
            
            value_zero = np.mean([listx[j],listx[j+1]])
            zeros_im.append(value_zero)





    if len(zeros_re) == 0:
#        print('no hay zeros')
        termEx_num_int_re3,err = integrate.quad(termEx_num_re_f, cota_inf_div , cota_sup_div, limit = 50)   

    elif len(zeros_re) == 1 :
#        print('hay 1 zeros')
#        print(zeros[0])
        
        termEx_num_int_re3_1,err = integrate.quad(termEx_num_re_f, cota_inf_div , zeros_re[0], limit = 5)
        
        termEx_num_int_re3_2,err = integrate.quad(termEx_num_re_f, zeros_re[0] , cota_sup_div, limit = 5)
        
        termEx_num_int_re3 = termEx_num_int_re3_1 + termEx_num_int_re3_2
 

    elif len(zeros_re) == 2 :
#        print('hay 2 zeros')
#        print(zeros[0])
        
        termEx_num_int_re3_1,err = integrate.quad(termEx_num_re_f, cota_inf_div , zeros_re[0], limit = 5)
        
        termEx_num_int_re3_2,err = integrate.quad(termEx_num_re_f, zeros_re[0] , zeros_re[1], limit = 5)
        
        termEx_num_int_re3_3,err = integrate.quad(termEx_num_re_f, zeros_re[1] , cota_sup_div, limit = 5)
        
        termEx_num_int_re3 = termEx_num_int_re3_1 + termEx_num_int_re3_2 + termEx_num_int_re3_3
       
    else  :
#        print('hay muchos zeros')
        term0_im = 0
        termEx_num_int_re3_1, err = integrate.quad(termEx_num_re_f, cota_inf_div, zeros_re[0], limit = 5)
        for j in range(len(zeros_re)-1):
            termEx_num_int_re3_2, err = integrate.quad(termEx_num_re_f, zeros_re[j], zeros_re[j+1], limit = 5)
            
            term0_im = term0_im + termEx_num_int_re3_2
        
        termEx_num_int_re3_3,err =  integrate.quad(termEx_num_re_f, zeros_re[-1], cota_sup_div, limit = 5)
        termEx_num_int_re3 = termEx_num_int_re3_1 + term0_im + termEx_num_int_re3_3 
        #    termEx_num_int_im3,err = integrate.quad(termEx_num_im_f, 0 , 0.99, limit = 50)   
        



    if len(zeros_im) == 0:
#        print('no hay zeros')
        termEx_num_int_im3,err = integrate.quad(termEx_num_im_f, cota_inf_div , cota_sup_div, limit = 50)   

    elif len(zeros_im) == 1 :
#        print('hay 1 zeros')
#        print(zeros[0])
        
        termEx_num_int_im3_1,err = integrate.quad(termEx_num_im_f, cota_inf_div , zeros_im[0], limit = 5)
        
        termEx_num_int_im3_2,err = integrate.quad(termEx_num_im_f, zeros_im[0] , cota_sup_div, limit = 5)
        
        termEx_num_int_im3 = termEx_num_int_im3_1 + termEx_num_int_im3_2
 

    elif len(zeros_im) == 2 :
#        print('hay 2 zeros')
#        print(zeros[0])
        
        termEx_num_int_im3_1,err = integrate.quad(termEx_num_im_f, cota_inf_div , zeros_im[0], limit = 5)
        
        termEx_num_int_im3_2,err = integrate.quad(termEx_num_im_f, zeros_im[0] , zeros_im[1], limit = 5)
        
        termEx_num_int_im3_3,err = integrate.quad(termEx_num_im_f, zeros_im[1] , cota_sup_div, limit = 5)
        
        termEx_num_int_im3 = termEx_num_int_im3_1 + termEx_num_int_im3_2 + termEx_num_int_im3_3
       
    else  :
#        print('hay muchos zeros')
        term0_im = 0
        termEx_num_int_im3_1, err = integrate.quad(termEx_num_im_f, cota_inf_div, zeros_im[0], limit = 5)
        for j in range(len(zeros_im)-1):
            termEx_num_int_im3_2, err = integrate.quad(termEx_num_im_f, zeros_im[j], zeros_im[j+1], limit = 5)
            
            term0_im = term0_im + termEx_num_int_im3_2
        
        termEx_num_int_im3_3,err =  integrate.quad(termEx_num_im_f, zeros_im[-1], cota_sup_div, limit = 5)
        termEx_num_int_im3 = termEx_num_int_im3_1 + term0_im + termEx_num_int_im3_3 

        
    termEx_num_int_3 = termEx_num_int_re3 +1j*termEx_num_int_im3
        
    
    return 1j*0.5*(term0_re + termEx_num_int_3 )*(k1**3)  


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
    eta = 1
    omegaWG = (np.pi/d)/(np.sqrt(epsilon2-1))  ##saco el c porque divido por c en omegac 
    omegac = omega_omegaWG*omegaWG
    k1 = omegac
    
    
    alpha_z1 = lambda u: np.sqrt(1 - u**2) if (1 > u) else 1j*np.sqrt(u**2 -1)   
#    alpha_z2 = lambda u: np.sqrt(epsilon2 - u**2) if (np.sqrt(epsilon2) > u) else 1j*np.sqrt(u**2 - epsilon2)   
        
    z_dip_barra = k1*(2*zp - z)
    
    exp_electron_f = lambda u: np.exp(1j*alpha_z1(u)*z_dip_barra)
    J0 =  lambda u : special.jv(0,u*k1*R)  
    J2 =  lambda u : special.jv(2,u*k1*R)

    term_Ex_px = lambda u: J0(u) + J2(u)*np.cos(2*phi)
    term_Ex_py = lambda u: J2(u)*np.sin(2*phi)


    list_all_modes = ['TE1A','TE1B','TE2A','TE2B']
    list_poles = []
    for type_of_mode in list_all_modes:  
        rta = function_poles(type_of_mode,omega_omegaWG,d)
        if rta != None:
            list_poles.append(rta)
    
    [inf_limits,sup_limits] = integration_limits(omega_omegaWG,d)
    term0_re = 0
    
    for j in range(len(list_poles)):

        polo_j = list_poles[j]
#            print(polo_j)
        qj_polo = polo_j/k1      
        den =  (polo_j*d)**2

        chi_prima = d*np.sqrt(epsilon2 - qj_polo**2)*k1
        chi = d*np.sqrt(qj_polo**2 - 1)*k1        
        den =  (polo_j*d)**2
        
        Rj_aux = 2*chi_prima*chi/den
        
        term1 = (chi**2 + chi_prima**2)/(chi*chi_prima)
        term2 =  ((eta*chi)**2 + chi_prima**2)/(eta*chi_prima)
        
        Rj_term_final = 1/(2*term1 + term2)
        
        Rj = Rj_aux*Rj_term_final
                
        
        rp = lambda u: Rj*qj_polo/(u - qj_polo - 1j*epsilon_pole) 

        
#        print(qj_polo)

        inf_limits_int = inf_limits[j]/k1
        sup_limits_int = sup_limits[j]/k1
        
        termEx_num_re_f = lambda u: np.real(u*rp(u)*(px*term_Ex_px(u) + py*term_Ex_py(u))*exp_electron_f(u)/alpha_z1(u))
        termEx_num_im_f = lambda u: np.imag(u*rp(u)*(px*term_Ex_px(u) + py*term_Ex_py(u))*exp_electron_f(u)/alpha_z1(u))

        termEx_num_int_im,err = integrate.quad(termEx_num_im_f, inf_limits_int, sup_limits_int, limit=limitt) 
 
        termEx_num_int = 1j*termEx_num_int_im
        
        term0_re = term0_re + termEx_num_int       

    return 1j*0.5*term0_re*(k1**3)  


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
    eta = 1
    omegaWG = (np.pi/d)/(np.sqrt(epsilon2-1))  ##saco el c porque divido por c en omegac 
    omegac = omega_omegaWG*omegaWG
    k0 = omegac


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
        J0 =  special.hankel1(0,polo_j*R)  
        J2 =   special.hankel1(2,polo_j*R)
        
#        print(polo_j*R,J0,J2)
        
#        J0 = np.real(J0)
#        J2 = np.real(J2)
#    
        rs = Rj*qj_polo  ## qj_polo*Rp
        
        term_Ex_px = J0 + J2*np.cos(2*phi)
        term_Ex_py = J2*np.sin(2*phi)
        

        termEx = qj_polo*rs*(px*term_Ex_px + py*term_Ex_py)*exp_electron_f/alpha_z1
        term0 = term0 + termEx
        

    return 1j*0.5*np.pi*1j*term0*(k0**3)

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
    eta = 1
    omegaWG = (np.pi/d)/(np.sqrt(epsilon2-1))  ##saco el c porque divido por c en omegac 
    omegac = omega_omegaWG*omegaWG
    k0 = omegac


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
        J0 =  special.hankel1(0,polo_j*R)  
        J2 =   special.hankel1(2,polo_j*R)
    
        rs = Rj*qj_polo  ## qj_polo*Rp
        
        term_Ex_py = J0 - J2*np.cos(2*phi)
        term_Ex_px = J2*np.sin(2*phi)
        

        termEx = qj_polo*rs*(px*term_Ex_px + py*term_Ex_py)*exp_electron_f/alpha_z1
        term0 = term0 + termEx
        
    return 1j*0.5*np.pi*1j*term0*(k0**3)
 

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
    
    alpha_z1 = lambda u: np.sqrt(1 - u**2) if (1 > u) else 1j*np.sqrt(u**2 -1)   
    alpha_z2 = lambda u: np.sqrt(epsilon2 - u**2) if (np.sqrt(epsilon2) > u) else 1j*np.sqrt(u**2 - epsilon2)   
        
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
        
    term0_re = 0
    
    [inf_limits,sup_limits] = integration_limits(omega_omegaWG,d)
    for j in range(len(sup_limits)):
    
        
        inf_limits_int = inf_limits[j]/k1
        sup_limits_int = sup_limits[j]/k1
        

        termEy_num_int_re,err = integrate.quad(termEx_num_re_f, inf_limits_int, sup_limits_int,  limit=limitt) 
        termEy_num_int_im,err = integrate.quad(termEx_num_im_f, inf_limits_int, sup_limits_int,  limit=limitt) 

        termEx_num_int = termEy_num_int_re + 1j*termEy_num_int_im
        
        term0_re = term0_re + termEx_num_int

    termEx_num_int_re3,err = integrate.quad(termEx_num_re_f,  0.001 , 0.999, limit = 5)   
    termEx_num_int_im3,err = integrate.quad(termEx_num_im_f,  0.001 , 0.999, limit = 5)   
    termEx_num_int_3 = termEx_num_int_re3  + 1j*termEx_num_int_im3

    return 1j*0.5*(term0_re + termEx_num_int_3)*(k1**3)  

#%%
#
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
    eta = 1

    omegaWG = (np.pi/d)/(np.sqrt(epsilon2-1))  ##saco el c porque divido por c en omegac 
    omegac = omega_omegaWG*omegaWG
    k1 = omegac
    
    
    alpha_z1 = lambda u: np.sqrt(1 - u**2) if (1 > u) else 1j*np.sqrt(u**2 -1)   
    alpha_z2 = lambda u: np.sqrt(epsilon2 - u**2) if (np.sqrt(epsilon2) > u) else 1j*np.sqrt(u**2 - epsilon2)   
        
    z_dip_barra = k1*(2*zp - z)
    
    exp_electron_f = lambda u: np.exp(1j*alpha_z1(u)*z_dip_barra)
    J0 =  lambda u : special.jv(0,u*k1*R)  
    J2 =  lambda u : special.jv(2,u*k1*R)


    
    term_Ex_py = lambda u: J0(u) - J2(u)*np.cos(2*phi)
    term_Ex_px = lambda u: J2(u)*np.sin(2*phi)
    

    list_all_modes = ['TE1A','TE1B','TE2A','TE2B']
    list_poles = []
    for type_of_mode in list_all_modes:  
        rta = function_poles(type_of_mode,omega_omegaWG,d)
        if rta != None:
            list_poles.append(rta)


    [inf_limits,sup_limits] = integration_limits(omega_omegaWG,d)
            
    term0_re = 0
    for j in range(len(list_poles)):

        polo_j = list_poles[j]
#            print(polo_j)
        qj_polo = polo_j/k1      
        den =  (polo_j*d)**2

        chi_prima = d*np.sqrt(epsilon2 - qj_polo**2)*k1
        chi = d*np.sqrt(qj_polo**2 - 1)*k1        
        den =  (polo_j*d)**2
        
        Rj_aux = 2*chi_prima*chi/den
        
        term1 = (chi**2 + chi_prima**2)/(chi*chi_prima)
        term2 =  ((eta*chi)**2 + chi_prima**2)/(eta*chi_prima)
        
        Rj_term_final = 1/(2*term1 + term2)
        
        Rj = Rj_aux*Rj_term_final
                
        
        rs = lambda u: Rj*qj_polo/(u - qj_polo - 1j*epsilon_pole) 
    
        termEx_num_re_f = lambda u: np.real(u*rs(u)*(px*term_Ex_px(u) + py*term_Ex_py(u))*exp_electron_f(u)/alpha_z1(u))
        termEx_num_im_f = lambda u: np.imag(u*rs(u)*(px*term_Ex_px(u) + py*term_Ex_py(u))*exp_electron_f(u)/alpha_z1(u))
 
        inf_limits_int = inf_limits[j]/k1
        sup_limits_int = sup_limits[j]/k1
        
    
        termEy_pole_aprox_int_re,err = integrate.quad(termEx_num_re_f, inf_limits_int, sup_limits_int,  limit=limitt) 
        termEy_pole_aprox_int_im,err = integrate.quad(termEx_num_im_f, inf_limits_int, sup_limits_int,  limit=limitt) 

        termEy_pole_aprox_int = termEy_pole_aprox_int_re + 1j*termEy_pole_aprox_int_im
        
        term0_re = term0_re + termEy_pole_aprox_int
    
    return 1j*0.5*term0_re*(k1**3)  

#
##%%
#    

