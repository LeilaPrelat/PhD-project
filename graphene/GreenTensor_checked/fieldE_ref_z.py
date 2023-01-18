
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila

campo externo directo numerico con la convencion de z hacia abajo
en z = 0
No hay solucion analitica porque es para cualquier x,y (solo hay sol analitica en x=y=0)
Dipolo en el 0,0
"""
from scipy import integrate
from scipy import special
import numpy as np
import sys
import os 

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/GreenTensor_checked','')
#print('Importar modulos necesarios para este codigo')

try:
    sys.path.insert(1, path_basic)
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

def Efield_ref_num(omegac,epsi1,epsi2,hbmu,hbgama,x,z,xe,zp,int_v):
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
    numerico habiendo aplicado QE al green tensor
    al momento de calcular el campo E. 
    """
    E = omegac*aux    
    # x_y = ky/k0
#    n1 = epsi1*mu1
#    cte1 = np.sqrt(n1)
#    k1 = omegac*cte1
#    k1_2 = k1**2
    
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    
    alpha_p = 1j*(epsi1 + epsi2)/(cond)
    Rp = 2*epsi1/(epsi1 + epsi2)
    kp = alpha_p*omegac
    term_pole = lambda alpha_parallel : Rp*alpha_p/(alpha_parallel-alpha_p)    

    J0 = lambda alpha_parallel: special.jn(0,alpha_parallel*omegac*np.abs(x-xe))

    term_den = lambda alpha_parallel: alpha_parallel**2 + int_v**2
    
    exp_electron = lambda alpha_parallel: np.exp(-alpha_parallel*omegac*(np.abs(z) + 2*zp))



    cota_inf1 = 0.01
    cota_sup2 = 1000
    
############ I0 5 ########################################################
    Int05_B_function_re = lambda u: u**3*np.real(term_pole(u)*J0(u))*exp_electron(u)/term_den(u)
    Int05_B_function_im = lambda u: u**3*np.imag(term_pole(u)*J0(u))*exp_electron(u)/term_den(u)

    int05B_re,err = integrate.quad(Int05_B_function_re, cota_inf1, cota_sup2)
    int05B_im,err = integrate.quad(Int05_B_function_im, cota_inf1, cota_sup2)
#########################################################################

    v = c/int_v
    sign_v = np.sign(v)
    cte_aux  = 1j*sign_v
    rta05 = cte_aux*(int05B_re + 1j*int05B_im) 
    
    return -rta05*(omegac**2)


#%%

def Efield_ref_fresnel(omegac,epsi1,epsi2,hbmu,hbgama,x,z,xe,zp,int_v):
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
    numerico habiendo aplicado QE al green tensor
    al momento de calcular el campo E. 
    """
    E = omegac*aux    
    # x_y = ky/k0
#    n1 = epsi1*mu1
#    cte1 = np.sqrt(n1)
#    k1 = omegac*cte1
#    k1_2 = k1**2
    
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    
    rp_num = lambda u: epsi2*1j*u - epsi1*1j*u - cond*u**2
    rp_den = lambda u: epsi2*1j*u + epsi1*1j*u - cond*u**2
    rp = lambda u: rp_num(u)/rp_den(u)

    J0 = lambda alpha_parallel: special.jn(0,alpha_parallel*omegac*np.abs(x-xe))

    term_den = lambda alpha_parallel: alpha_parallel**2 + int_v**2
    
    exp_electron = lambda alpha_parallel: np.exp(-alpha_parallel*omegac*(np.abs(z) + 2*zp))



    cota_inf1 = 0.01
    cota_sup2 = 1000
    
############ I0 5 ########################################################
    Int05_B_function_re = lambda u: u**3*np.real(rp(u)*J0(u))*exp_electron(u)/term_den(u)
    Int05_B_function_im = lambda u: u**3*np.imag(rp(u)*J0(u))*exp_electron(u)/term_den(u)

    int05B_re,err = integrate.quad(Int05_B_function_re, cota_inf1, cota_sup2)
    int05B_im,err = integrate.quad(Int05_B_function_im, cota_inf1, cota_sup2)
#########################################################################

    v = c/int_v
    sign_v = np.sign(v)
    cte_aux  = 1j*sign_v
    rta05 = cte_aux*(int05B_re + 1j*int05B_im) 
    
    return -rta05*(omegac**2)

#%%

def Efield_ref_ana(omegac,epsi1,epsi2,hbmu,hbgama,x,z,xe,zp,int_v):
    """
    Parameters
    ----------
    omegac : omega/c = k0 en 1/micrometros    
    epsi1 : epsilon del medio de arriba del plano
    int_v : fraccion de la vel de la c que vale la velocidad uniforme del electron en x, si v = c/int
    b : electron position en z, b < 0 (dipolo en 0 y eje de z apunta hacia abajo)  
    Returns
    -------
    External field ref / (e/omega)
    """
    E = omegac*aux    
    # x_y = ky/k0
#    n1 = epsi1*mu1
#    cte1 = np.sqrt(n1)
#    k1 = omegac*cte1
#    k1_2 = k1**2
    

    
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    
    alpha_p = 1j*(epsi1 + epsi2)/(cond)
    Rp = 2*epsi1/(epsi1 + epsi2)
    kp = alpha_p*omegac


    H0 = special.hankel1(0,kp*np.abs(x-xe))


    v = c/int_v
    sign_v = np.sign(v)
    kp_4 = kp**4
        
    den_term = kp**2 + (omegac*int_v)**2
    
    expo_term = np.exp(-kp*(np.abs(z) + 2*zp))
    
    rta_final = np.pi*H0*expo_term*sign_v*Rp*kp_4/den_term
    
    return rta_final


#%%


