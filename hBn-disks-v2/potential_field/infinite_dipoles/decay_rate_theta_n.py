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
path_constants =  path_basic.replace('/potential_field/infinite_dipoles','')
#print('Importar modulos necesarios para este codigo')

try:
    sys.path.insert(1, path_constants)
    from hBn_PP import epsilon_x,hBn_lambda_p,hBn_Rp
except ModuleNotFoundError:
    print('hBn_PP.py no se encuentra en ' + path_constants)


try:
    sys.path.insert(1, path_basic)
    from dipole_moment_dif_sign import dipole_moment_pole_aprox_resonance_v1,dipole_moment_ana_resonance_v1,dipole_moment_anav1_for_decay_rate_resonance,dipole_moment_pole_aprox_resonance_v1_for_decay_rate,dipole_moment_ana_perfect_dipole_resonance
except ModuleNotFoundError:
    print('potential.py no se encuentra en ' + path_basic)



try:
    sys.path.insert(1, path_basic)
    from potential_induced import potential_ana_resonance_v1,potential_num_resonance
except ModuleNotFoundError:
    print('potential_induced.py no se encuentra en ' + path_basic)
    
try:
    sys.path.insert(1, path_constants)
    from constants import constantes
except ModuleNotFoundError:
    print('constants.py no se encuentra en ' + path_constants)

pi,hb,c,alfac,mu1,mu2 = constantes()
aux = c*hb

#%%


def decay_rate_theta_inf_dipoles_ana_res(omegac,epsi_silica,d_nano,int_v,zp,a,b,n):     
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

    x, y, z = 0,0,0
    E = omegac*aux
    n1 = epsi1*mu1
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1
   
#    d_micros = d_nano*1e-3
#    Rp = hBn_Rp(E,epsi1,epsi3)
#    lambda_p_v = hBn_lambda_p(E,epsi1,epsi3)*d_micros
#    kp = 2*np.pi/lambda_p_v
#    alfa_p = kp/omegac 
    
    #    list_dipoles = np.linspace(-Nmax,Nmax,2*Nmax + 1)
    Rp = 1
#    theta_degree = 45
#    theta = theta_degree*np.pi/180
    px,py,pz  = dipole_moment_pole_aprox_resonance_v1(omegac,epsi_silica,d_nano,int_v,b,zp)  
#    list_dipoles = np.linspace(-Nmax,Nmax,2*Nmax + 1)
#            
    kx = omegac*int_v + 2*np.pi*n/a     
#    expo_kx = np.exp(1j*kx*x)
#
    ### cte que falta en dip moment : e/(2*np.pi*v)
    charge_e_cgs = 1.602176634*1e-20
    
    hbar_cgs = 1.054571817*1e-27
    
    cte_dip = int_v/(2*np.pi)
    
    px,py,pz = px*cte_dip, py*cte_dip, pz*cte_dip
    
    cte_final = (charge_e_cgs/(hbar_cgs*c))**2
#    print(cte_final)
#    
#    ky = kp*np.sin(theta)

#    cte = 1/((2*np.pi)**2*a)
    
#    cte2 = alfac*int_v*1e15/(np.pi) ## cambio de unidades + agregar lo que faltaba en el momento dipolar
    x, y, z = 0,0,0
    phi_n = potential_num_resonance(omegac,epsi_silica,d_nano,int_v,b,zp,x,y,z,a,n)
    
    rta = cte_final*a*np.abs(phi_n)**2/(2*np.pi*np.real(Rp))
    
#    try: 
#        
#        rta>0
#        
#    except ValueError: 
#        print(zp/lambda_p_v)
    
    return rta



def decay_rate_theta_inf_dipoles_ana_res_div_gamma0(omegac,epsi_silica,d_nano,int_v,zp,a,b,n):     
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

#    x, y, z = 0,0,0
    E = omegac*aux

    k1 = omegac
    
    epsi_x = epsilon_x(E)
    
#    if np.imag(np.sqrt(epsi_x*epsi_z))>0:
#        epsi_m = np.sqrt(epsi_x*epsi_z)
#    else:
#        epsi_m = -np.sqrt(epsi_x*epsi_z)
#    

    epsi_HBN_par = epsi_x
   
    d_micro = d_nano*1e-3
    alfa_p = epsi_silica(E)*2/(omegac*d_micro*(epsi_HBN_par-1))
    kp = alfa_p*omegac
    Rp = 1


#        
    
    px,py,pz  = dipole_moment_pole_aprox_resonance_v1(omegac,epsi_silica,d_nano,int_v,b,zp)  
#    list_dipoles = np.linspace(-Nmax,Nmax,2*Nmax + 1)
#            
    kx = omegac*int_v + 2*np.pi*n/a     
#    expo_kx = np.exp(1j*kx*x)
#
    ### cte que falta en dip moment : e/(2*np.pi*v)
#    charge_e_cgs = 1.602176634*1e-20
#    
#    hbar_cgs = 1.054571817*1e-27
    
    v = c/int_v
    cte_dip = 1/(2*np.pi*v)
    
    px,py,pz = px*cte_dip, py*cte_dip, pz*cte_dip
    
#    cte_final = (charge_e_cgs/(hbar_cgs*c))**2
#    print(cte_final)
#    
#    ky = kp*np.sin(theta)

#    cte = 1/((2*np.pi)**2*a)
    
#    cte2 = alfac*int_v*1e15/(np.pi) ## cambio de unidades + agregar lo que faltaba en el momento dipolar
#    den = np.sqrt(kp**2 - kx**2)
   # return np.imag(final_2*cte*kp*np.cos(theta))
#   
#    x, y, z = 0,0,0
#    phi_n = potential_ana_resonance_v1(omegac,epsi_silica,d_nano,int_v,b,zp,x,y,z,a,n)
    
    den = np.sqrt(kp**2 - kx**2)
   # return np.imag(final_2*cte*kp*np.cos(theta))
    phi_n = np.exp(-2*kp*zp)*Rp*kp*(px*kx/den + py + 1j*pz*kp/den )/(2*np.pi*a)    
    
    
    
    rta = a*np.abs(phi_n)**2/(2*np.pi*np.abs(Rp))
    
#    cte_aux = cte_aux*1e9 ### cambiar unidades

    gamma = np.sqrt(1 - (int_v)**(-2))**(-1)
    gamma = 1
    alpha = -3*epsi_silica(E)/(4*1j*k1**3)

    arg = np.abs(b)*omegac*int_v/gamma
    K1 = special.kn(1,arg)
    K0 = special.kn(0,arg)
    
#    print(K1)
   

    factor_gamma0 = (2*omegac*int_v/(v*gamma))**2
    gamma0 = factor_gamma0*(K0**2/gamma**2 + K1**2)*np.imag(alpha)/np.pi  ## decay rate de 1 dipolo # pero sin el "e/hbar" se cancela con el momento dipolar^2
    
    

    return np.real(rta/(gamma0*hb))






# normalizado con el paper 149 
def decay_rate_theta_inf_dipoles_ana_res_div_gamma0_v3(omegac,epsi_silica,d_nano_film,d_thickness_disk_nano,D_disk_nano,int_v,zp,a,b,n):     
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

#    x, y, z = 0,0,0
    E = omegac*aux
    
#    d_micros = d_nano*1e-3
#    Rp = hBn_Rp(E,epsi_silica(E),epsi_silica(E))
#    lambda_p_v = hBn_lambda_p(E,epsi_silica(E),epsi_silica(E))*d_micros
#    kp = 2*np.pi/lambda_p_v
#    alfa_p = kp/omegac 
    
    
#    Rp = 1
#    epsi_x = epsilon_x(E)
#    epsi_HBN_par = epsi_x
#    d_micro = d_nano_film*1e-3
#    alfa_p = epsi_silica(E)*2/(omegac*d_micro*(epsi_HBN_par-1))
#    kp = alfa_p*omegac

    d_micros = d_nano_film*1e-3
    Rp = hBn_Rp(E,epsi_silica(E),epsi_silica(E))
    lambda_p_v = hBn_lambda_p(E,epsi_silica(E),epsi_silica(E))*d_micros
    kp = 2*np.pi/lambda_p_v
    alfa_p = kp/omegac 
#    
    px,py,pz  = dipole_moment_ana_perfect_dipole_resonance(omegac,epsi_silica,d_nano_film,d_thickness_disk_nano,D_disk_nano,int_v,b,zp)  
#    list_dipoles = np.linspace(-Nmax,Nmax,2*Nmax + 1)
#            
    kx = omegac*int_v + 2*np.pi*n/a     
    den = np.sqrt(kp**2 - kx**2)
   # return np.imag(final_2*cte*kp*np.cos(theta))

    kp_2 = np.sqrt(kp**2)
#    term_kp = 1 + kp/kp_2
    term_kp_2 = kp_2 + kp
    
#    term_kp_2 = 2*kp
#    kp_2 = kp
    
    phi_n = np.exp(-2*kp_2*zp)*np.abs(Rp)*(px*kx*term_kp_2/den + py*term_kp_2 + 1j*pz*kp_2*term_kp_2/den )/(4*np.pi*a)

    
    cte_formula = a/(12*np.pi*np.abs(Rp)) ## hay un extra 1/(2pi) en la formula de phi. necesario para silver 
    cte_formula = 12*2*(np.pi**3)*a/(Rp) ## hay un extra 1/(2pi) en la formula de phi. necesario para silver 
    cte_formula = 128*(np.pi**7)*a/(Rp)    
#      
#    cte_formula = a/(12*Rp) ## hay un extra 1/(2pi) en la formula de phi
#    
#    cte_formula = a*192*2*np.pi**4/np.abs(Rp)  ## hay un extra (1/(2pi))^N en la formula de phi. necesario para grafeno  
##
#

    gamma = (1 - (int_v**(-2)) )**(-1/2)
    arg = np.abs(b)*omegac*int_v/gamma
    K1 = special.kn(1,arg)
    K0 = special.kn(0,arg)

    
    factor_K = K0**2 + K1**2  ## decay rate de 1 dipolo # pero sin el "e/hbar" se cancela con el momento dipolar^2
    extra_cte_adimensional = a*omegac ## para poder comparar con diferentes materiales y que no dependa del periodo "a" 
    seno_theta_n = den/kp  
    
#    px_dir,py_dir,pz_dir = dipole_moment_anav2_for_decay_rate_resonance_dir(omegac,int_v,b,zp)        
#    denominador = np.abs(px_dir)**2 +  np.abs(py_dir)**2 +  np.abs(pz_dir)**2

    k_prima = omegac*np.sqrt(epsi_silica(E))
#    epsi1 = 1
#    k_prima = omegac*np.sqrt(epsi1)
#        
    rta = (np.abs(phi_n)**2)*cte_formula*k_prima*(int_v**(-2))/(factor_K*seno_theta_n)        
        
    return rta



