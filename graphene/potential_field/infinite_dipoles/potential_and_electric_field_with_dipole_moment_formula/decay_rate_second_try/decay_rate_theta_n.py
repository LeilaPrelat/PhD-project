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
path_constants =  path_basic.replace('/potential_field/infinite_dipoles/potential_and_electric_field_with_dipole_moment_formula/decay_rate_second_try','')
#print('Importar modulos necesarios para este codigo')

try:
    sys.path.insert(1, path_constants)
    from graphene_sigma import sigma_DL
except ModuleNotFoundError:
    print('graphene_sigma.py no se encuentra en ' + path_basic)

try:
    sys.path.insert(1, path_basic)
    from dipole_moment import dipole_moment_sin_integrar_en_y, dipole_moment_sin_integrar_en_y_resonance, dipole_moment_anav2_res,dipole_moment_anav2_for_decay_rate_res, dipole_moment_anav2_for_decay_rate_resonance_dir
except ModuleNotFoundError:
    print('dipole_moment.py no se encuentra en ' + path_basic)
    
try:
    sys.path.insert(1, path_constants)
    from constants import constantes
except ModuleNotFoundError:
    print('constants.py no se encuentra en ' + path_constants)

pi,hb,c,alfac,mu1,mu2 = constantes()
aux = c*hb

#%%


def decay_rate_theta_inf_dipoles_ana_old(omegac,epsi1,epsi2,hbmu,hbgama,int_v,zp,a,b,n,omega0,kappa_factor_omega0,kappa_r_factor,theta):     
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

    x, y, z = 0,0,zp 
    E = omegac*aux
    n1 = epsi1*mu1
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1
   
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    Rp = 2*epsi1/(epsi1 + epsi2)
    alfa_p = 1j*(epsi1 + epsi2)/(cond*cte1)
    kp = alfa_p*k1
    
    
    px,py,pz  = dipole_moment_sin_integrar_en_y(omegac,epsi1,epsi2,hbmu,hbgama,int_v,b,zp,theta,omega0,kappa_factor_omega0,kappa_r_factor)  
#    list_dipoles = np.linspace(-Nmax,Nmax,2*Nmax + 1)
#            
    kx = omegac*int_v + 2*np.pi*n/a     
    expo_kx = np.exp(1j*kx*x)
    
    
    ky = kp*np.sin(theta)
    term_den = np.sqrt(ky**2 + kx**2)
#    den_dir = np.sqrt(np.abs(z)**2 + np.abs(y)**2)
#    K0 = special.kn(0,kx*den_dir)
#    K1 = special.kn(1,kx*den_dir)
    
    exp_electron = np.exp(-term_den*(2*zp - z))*np.exp(1j*ky*np.abs(y))*expo_kx

#    term1 =  -2*1j*px*kx*K0 + 2*py*kx*np.abs(y)*K1/den_dir  + pz*np.sign(z)*2*kx*np.abs(z)*K1/den_dir


    rp = Rp*kp/(term_den - kp)
    term2 = 1j*px*kx*rp/term_den
    
    term3 = 1j*py*ky*rp/term_den
    
    term4 = -pz*kp*rp
    
 

    final =  (term2 + term3 + term4)*exp_electron

    Ex = 1j*kx*final
    Ey = 1j*ky*final
    Ez = term_den*final

    final_2 = np.conjugate(px)*Ex + np.conjugate(py)*Ey + np.conjugate(pz)*Ez


    cte = 1/((2*np.pi)**2*a)
    
   # return np.imag(final_2*cte*kp*np.cos(theta))

#    cte2 = alfac*int_v*1e15/(np.pi) ## cambio de unidades + agregar lo que faltaba en el momento dipolar
    
    return np.imag(final_2*cte)


#%%



def decay_rate_theta_inf_dipoles_ana_res_old(omegac,epsi1,epsi2,hbmu,hbgama,int_v,zp,a,b,n,theta):     
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
   
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    Rp = 2*epsi1/(epsi1 + epsi2)
    alfa_p = 1j*(epsi1 + epsi2)/(cond*cte1)
    kp = alfa_p*k1
    
    
    px,py,pz  = dipole_moment_sin_integrar_en_y_resonance(omegac,epsi1,epsi2,hbmu,hbgama,int_v,b,zp,theta)  
#    list_dipoles = np.linspace(-Nmax,Nmax,2*Nmax + 1)
#            
    kx = omegac*int_v + 2*np.pi*n/a     
    expo_kx = np.exp(1j*kx*x)
    
    
    ky = kp*np.sin(theta)
    term_den = np.sqrt(ky**2 + kx**2)
#    den_dir = np.sqrt(np.abs(z)**2 + np.abs(y)**2)
#    K0 = special.kn(0,kx*den_dir)
#    K1 = special.kn(1,kx*den_dir)
    
    exp_electron = np.exp(-term_den*(2*zp - z))*np.exp(1j*ky*np.abs(y))*expo_kx

#    term1 =  -2*1j*px*kx*K0 + 2*py*kx*np.abs(y)*K1/den_dir  + pz*np.sign(z)*2*kx*np.abs(z)*K1/den_dir


    rp = Rp*kp/(term_den - kp)
    term2 = 1j*px*kx*rp/term_den
    
    term3 = 1j*py*ky*rp/term_den
    
    term4 = -pz*kp*rp
    
 

    final =  (term2 + term3 + term4)*exp_electron

    Ex = 1j*kx*final
    Ey = 1j*ky*final
    Ez = term_den*final

    final_2 = np.conjugate(px)*Ex + np.conjugate(py)*Ey + np.conjugate(pz)*Ez


    cte = 1/((2*np.pi)**2*a)
    
#    cte2 = alfac*int_v*1e15/(np.pi) ## cambio de unidades + agregar lo que faltaba en el momento dipolar
    
   # return np.imag(final_2*cte*kp*np.cos(theta))

    
    return np.imag(final_2*cte)


#%%


def decay_rate_theta_inf_dipoles_ana_res(omegac,epsi1,epsi2,hbmu,hbgama,int_v,zp,a,b,n):     
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
    n1 = epsi1*mu1
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1
   
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    Rp = 2*epsi1/(epsi1 + epsi2)
    alfa_p = 1j*(epsi1 + epsi2)/(cond*cte1)
    kp = alfa_p*k1
    
    
    px,py,pz  = dipole_moment_anav2_res(omegac,epsi1,epsi2,hbmu,hbgama,int_v,b,zp)  
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
    den = np.sqrt(kp**2 - kx**2)
   # return np.imag(final_2*cte*kp*np.cos(theta))
    phi_n = np.exp(-2*kp*zp)*Rp*kp*(px*kx/den + py + 1j*pz*kp/den )/(2*np.pi*a)
    
    rta = cte_final*a*np.abs(phi_n)**2/(2*np.pi*Rp)
    
    return rta





def decay_rate_theta_inf_dipoles_ana_res_div_gamma0(omegac,epsi1,epsi2,hbmu,hbgama,int_v,zp,a,b,n):     
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
    n1 = epsi1*mu1
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1
   
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    Rp = 2*epsi1/(epsi1 + epsi2)
    alfa_p = 1j*(epsi1 + epsi2)/(cond*cte1)
    kp = alfa_p*k1
    
    
    px,py,pz  = dipole_moment_anav2_for_decay_rate_res(omegac,epsi1,epsi2,hbmu,hbgama,int_v,b,zp)  
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
    den = np.sqrt(kp**2 - kx**2)
   # return np.imag(final_2*cte*kp*np.cos(theta))
    phi_n = np.exp(-2*kp*zp)*Rp*kp*(px*kx/den + py + 1j*pz*kp/den )/(2*np.pi*a)
    
    rta = a*np.abs(phi_n)**2/(2*np.pi*Rp)
    
#    cte_aux = cte_aux*1e9 ### cambiar unidades

    gamma = np.sqrt(1 - (int_v)**(-2))**(-1)
    alpha = -3*epsi1/(4*1j*k1**3)

    arg = np.abs(b)*omegac*int_v/gamma
    K1 = special.kn(1,arg)
    K0 = special.kn(0,arg)
    
#    print(K1)
   

    factor_gamma0 = (2*omegac*int_v/(v*gamma))**2
    gamma0 = factor_gamma0*(K0**2/gamma**2 + K1**2)*np.imag(alpha)/np.pi  ## decay rate de 1 dipolo # pero sin el "e/hbar" se cancela con el momento dipolar^2
    
    return rta/(gamma0*hb)



def decay_rate_theta_inf_dipoles_ana_res_div_gamma0_v2(omegac,epsi1,epsi2,hbmu,hbgama,int_v,zp,a,b,n):     
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
    n1 = epsi1*mu1
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1
   
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    Rp = 2*epsi1/(epsi1 + epsi2)
    alfa_p = 1j*(epsi1 + epsi2)/(cond*cte1)
    kp = alfa_p*k1
    
    
    px,py,pz  = dipole_moment_anav2_for_decay_rate_res(omegac,epsi1,epsi2,hbmu,hbgama,int_v,b,zp)  
#    list_dipoles = np.linspace(-Nmax,Nmax,2*Nmax + 1)
#            
    kx = omegac*int_v + 2*np.pi*n/a     
#    expo_kx = np.exp(1j*kx*x)
#
    ### cte que falta en dip moment : e/(2*np.pi*v)
#    charge_e_cgs = 1.602176634*1e-20
#    
#    hbar_cgs = 1.054571817*1e-27

    
#    px,py,pz = px*cte_dip, py*cte_dip, pz*cte_dip
    
#    cte_final = (charge_e_cgs/(hbar_cgs*c))**2
#    print(cte_final)
#    
#    ky = kp*np.sin(theta)

#    cte = 1/((2*np.pi)**2*a)
    
#    cte2 = alfac*int_v*1e15/(np.pi) ## cambio de unidades + agregar lo que faltaba en el momento dipolar
    den = np.sqrt(kp**2 - kx**2)
   # return np.imag(final_2*cte*kp*np.cos(theta))
    phi_n = np.exp(-2*kp*zp)*Rp*kp*(px*kx/den + py + 1j*pz*kp/den )/(2*np.pi*a)

    cte_formula = a/(2*np.pi*Rp)
    
#    usar_dif_p = 1
#    usar_mismo_p = 0
#    if usar_dif_p == 1:  
        
    px_dir,py_dir,pz_dir = dipole_moment_anav2_for_decay_rate_resonance_dir(omegac,int_v,b,zp)
    
    denominador = np.abs(px_dir)**2 +  np.abs(py_dir)**2 +  np.abs(pz_dir)**2

    rta = np.abs(phi_n)**2*(omegac**3)*cte_formula/denominador    
    
#    if usar_mismo_p == 1: ## aca los momentos p se cancelan  
#        
#        px_dir,py_dir,pz_dir = dipole_moment_anav2_for_decay_rate_res(omegac,epsi1,epsi2,hbmu,hbgama,int_v,b,zp) # multiplicar por e/(2*pi*v)
#
#        denominador = np.abs(px_dir)**2 +  np.abs(py_dir)**2 +  np.abs(pz_dir)**2
#
#        alffa_eff_x = 1j*(2*omegac**3/(3*epsi1) +  np.imag(rtaself_x))**(-1)
#        alffa_eff_y = 1j*(2*omegac**3/(3*epsi1) +  np.imag(rtaself_y))**(-1)
#        alffa_eff_z = 1j*(2*omegac**3/(3*epsi1) +  np.imag(rtaself_z))**(-1)
#    
#        alfa_eff = np.abs(alffa_eff_x)**2 + np.abs(alffa_eff_y)**2 + np.abs(alffa_eff_z)**2 
#
#    
#        rta = np.imag(Green_self*alfa_eff*2*(omegac**3)/denominador )

    return rta 





def decay_rate_theta_inf_dipoles_ana_res_div_gamma0_angle(omegac,epsi1,epsi2,hbmu,hbgama,int_v,zp,a,b,n):     
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
    n1 = epsi1*mu1
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1
   
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    Rp = 2*epsi1/(epsi1 + epsi2)
    alfa_p = 1j*(epsi1 + epsi2)/(cond*cte1)
    kp = alfa_p*k1
    
    
    px,py,pz  = dipole_moment_anav2_res(omegac,epsi1,epsi2,hbmu,hbgama,int_v,b,zp)  
#    list_dipoles = np.linspace(-Nmax,Nmax,2*Nmax + 1)
#            
  
#    expo_kx = np.exp(1j*kx*x)
#
    ### cte que falta en dip moment : e/(2*np.pi*v)
#    charge_e_cgs = 1.602176634*1e-20
#    
#    hbar_cgs = 1.054571817*1e-27
    
    v = c/int_v
    cte_dip = 1/(2*np.pi*v)
    
    px,py,pz = px*cte_dip, py*cte_dip, pz*cte_dip
    

    def theta(E_meV,a):
    
        E = E_meV*1e-3  
        cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
        alfa_p = 1j*(epsi1 + epsi2)/(cond)
        omegac = E/(hb*c)
        kp = alfa_p*omegac    
        
        theta0 = np.arccos((omegac*int_v + 2*np.pi*n/a)/np.real(kp))
    #    print(theta0)
        return theta0*180/np.pi

    kx = omegac*int_v + 2*np.pi*n/a   
    den = np.sqrt(kp**2 - kx**2)  ### ky 
    
    Emev = E*1e3
    den = kp*np.sin(theta(Emev,a))
    kx = kp*np.cos(theta(Emev,a))
    
   # return np.imag(final_2*cte*kp*np.cos(theta))
    phi_n = np.exp(-2*kp*zp)*Rp*kp*(px*kx/den + py + 1j*pz*kp/den )/(2*np.pi*a)
    
    rta = a*np.abs(phi_n)**2/(2*np.pi*Rp)
    
#    cte_aux = cte_aux*1e9 ### cambiar unidades

    gamma = np.sqrt(1 - (int_v)**(-2))**(-1)
    alpha = -3*epsi1/(4*1j*k1**3)

    arg = np.abs(b)*omegac*int_v/gamma
    K1 = special.kn(1,arg)
    K0 = special.kn(0,arg)
    
#    print(K1)
   

    factor_gamma0 = (2*omegac*int_v/(v*gamma))**2
    gamma0 = factor_gamma0*(K0**2/gamma**2 + K1**2)*np.imag(alpha)/np.pi  ## decay rate de 1 dipolo # pero sin el "e/hbar" se cancela con el momento dipolar^2
    
    return rta/(gamma0*hb)

