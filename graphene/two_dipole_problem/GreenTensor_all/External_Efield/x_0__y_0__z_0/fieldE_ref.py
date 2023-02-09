
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila

campo externo reflejado en el plano con la convencion de z hacia abajo
en el punto x  = 0, y = 0,z = 0
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
path_constants =  path_basic.replace('/External_Efield/x_0__y_0__z_0','')
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

#%%

def Efield_NUM_QE(omegac,epsi1,epsi2,hbmu,hbgama,int_v,b,zp):
    """
    Parameters
    ----------
    omegac : omega/c = k0 en 1/micrometros    
    epsi1 : epsilon del medio de arriba del plano
    epsi2 : epsilon del medio de abajo del plano
    hbmu : chemical potential in eV  
    hbgama : collision frequency in eV
    int_v : fraccion de la vel de la c que vale la velocidad uniforme del electron en x, si v = c/int
    b: electron position en z, b < 0 (dipolo en 0 y eje de z apunta hacia abajo) 
    zp : coordenada zp del plano, zp > 0  (dipolo en 0 y eje de z apunta hacia abajo) 
    Returns
    -------
    External field reflected (green tensor reflected)
    numerico habiendo aplicado QE al green tensor. 
    Dipolo en el 0,0,0. PARA EL PUNTO 0,0,0
    """
    E = omegac*aux    
    omega = omegac*c #=omega/c    
    n1 = epsi1*mu1
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1
    # k1_2 = k1**2
    k1_2 = k1**2

# green tensor evaluated in : DIPOLO EN EL 0,0,0 
#    x,y,z = 0, 0, 0
#    yD,zD = 0, b

    # Rbarra = lambda w: k1*w # green tensor evaluated in x,y,z = 0, 0, 0 and yD,zD = 0, b
    z_dip_barra = k1*(np.abs(b) + 2*zp)  #NO ESTOY SEGURA DE ESTA EXPONENCIAL EH 
    atan = lambda w : np.cos(np.pi - 2*np.arctan(np.abs(w))) 

    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)

    rp_num = lambda u: epsi2*1j*u - epsi1*1j*u - cte1*cond*(u**2)
    rp_den = lambda u: epsi2*1j*u + epsi1*1j*u - cte1*cond*(u**2)
    rp = lambda u: rp_num(u)/rp_den(u)
    
    rs_num = lambda u: 1j*u - 1j*u - cond/cte1
    rs_den = lambda u: 1j*u + 1j*u + cond/cte1
    rs = lambda u: rs_num(u)/rs_den(u)

    expB = lambda u: np.exp(-u*z_dip_barra) 
    exp_electron = lambda w: np.exp(1j*w*int_v/n1)
    J0 = lambda w,u: special.jn(0,u*w) 
    J2 = lambda w,u: special.jn(2,u*w) 

    cota_sup1 = 1000
    cota_sup2 = 80*k1
############ I0 5 + I2 5 #################################################################################################
#    Int5_B_function_re = lambda w,u: np.real((J0(w,u) + J2(w,u))*rs(u)*expB(u)*exp_electron(w))
#    Int5_B_function_im = lambda w,u: np.imag((J0(w,u) + J2(w,u))*rs(u)*expB(u)*exp_electron(w))

#    int5B_re,err = integrate.dblquad(Int5_B_function_re, 0, cota_sup1, -cota_sup2, cota_sup2)
#    int5B_im,err = integrate.dblquad(Int5_B_function_im, 0, cota_sup1, -cota_sup2, cota_sup2)
############ I0 5 + I2 5 + I0 6 + I2 6 #################################################################################################
    IntB_function_re = lambda w,u: np.real((J0(w,u) + atan(w)*J2(w,u))*((u**2)*rp(u) + rs(u))*expB(u)*exp_electron(w))
    IntB_function_im = lambda w,u: np.imag((J0(w,u) + atan(w)*J2(w,u))*((u**2)*rp(u) + rs(u))*expB(u)*exp_electron(w))

    intB_re,err = integrate.dblquad(IntB_function_re, 1, cota_sup1, -cota_sup2, cota_sup2)
    intB_im,err = integrate.dblquad(IntB_function_im, 1, cota_sup1, -cota_sup2, cota_sup2)
##########################################################################################################################


    charge_e = 1.602*1e-19*c
    cte_aux = -1j*(charge_e*omega)/(2*np.pi*int_v)

    cte = k1_2*0.5*cte_aux #signo menos
    
#    rta5 = (int5B_re + 1j*int5B_im)*cte 
    rta = (intB_re + 1j*intB_im)*cte 
    
    return rta

#%%

def Efield_NUM(omegac,epsi1,epsi2,hbmu,hbgama,int_v,b,zp):
    """
    Parameters
    ----------
    omegac : omega/c = k0 en 1/micrometros    
    epsi1 : epsilon del medio de arriba del plano
    epsi2 : epsilon del medio de abajo del plano
    hbmu : chemical potential in eV  
    hbgama : collision frequency in eV
    int_v : fraccion de la vel de la c que vale la velocidad uniforme del electron en x, si v = c/int
    b: electron position en z, b < 0 (dipolo en 0 y eje de z apunta hacia abajo)
    zp : coordenada zp del plano, zp > 0  (dipolo en 0 y eje de z apunta hacia abajo) 
    Returns
    -------
    External field reflected (green tensor reflected)
    dividido por 1j*e/(2*pi)
    numerico sin haber aplicado QE al green tensor
    Dipolo en el 0,0,0
    """
    E = omegac*aux  
    omega = omegac*c #=omega/c    
    n1 = epsi1*mu1
    n2 = epsi2*mu2
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1
    k1_2 = k1**2

# green tensor evaluated in : DIPOLO EN EL 0,0,0 
#    x,y,z = 0, 0, 0
#    yD,zD = 0, b

#    Rbarra = lambda w: k1*w # green tensor evaluated in x,y,z = 0, 0, 0 and yD,zD = 0, b
    z_dip_barra = k1*(np.abs(b) + 2*zp)  #NO ESTOY SEGURA DE ESTA EXPONENCIAL EH , b
    atan = lambda w : np.cos(np.pi - 2*np.arctan(np.abs(w)))

#    alpha_z1 = lambda u: np.sqrt(1-u**2)
#    alpha_z2 = lambda u: np.sqrt(u**2-1)
    aux2 = n2/n1
    alpha_z1 = lambda u: np.sqrt(1-u**2) if u<1 else 1j*np.sqrt(u**2-1)
    alpha_z2 = lambda u: np.sqrt(aux2-u**2) if u<aux2 else 1j*np.sqrt(u**2-aux2)
        
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)

    rp_num = lambda u: epsi2*alpha_z1(u) - epsi1*alpha_z2(u) + cte1*cond*alpha_z1(u)*alpha_z2(u)
    rp_den = lambda u: epsi2*alpha_z1(u) + epsi1*alpha_z2(u) + cte1*cond*alpha_z1(u)*alpha_z2(u)
    rp = lambda u: rp_num(u)/rp_den(u)
    
    rs_num = lambda u: alpha_z1(u) - alpha_z2(u) - cond/cte1
    rs_den = lambda u: alpha_z1(u) + alpha_z2(u) + cond/cte1
    rs = lambda u: rs_num(u)/rs_den(u)

    exp = lambda u: np.exp(1j*alpha_z1(u)*z_dip_barra) 
    exp_electron = lambda w: np.exp(1j*w*int_v/n1)
    J0 = lambda w,u: special.jn(0,u*w) 
    J2 = lambda w,u: special.jn(2,u*w) 
    
    cota_sup1 = 0.95
    cota_inf1 = 1.05
    
    cota_sup1A = 1000
    cota_sup2A = 80*k1

############ I0 5 + I2 5  + I0 6 + I2 6 #################################################################################################
    Int_A_function_re = lambda w,u: np.real((J0(w,u) + atan(w)*J2(w,u))*u*exp(u)*exp_electron(w)*(-rp(u)*alpha_z1(u) + rs(u)/alpha_z1(u)))
    Int_A_function_im = lambda w,u: np.imag((J0(w,u) + atan(w)*J2(w,u))*u*exp(u)*exp_electron(w)*(-rp(u)*alpha_z1(u) + rs(u)/alpha_z1(u)))

    intA_re,err = integrate.dblquad(Int_A_function_re, 0, cota_sup1, -cota_sup2A, cota_sup2A)
    intA_im,err = integrate.dblquad(Int_A_function_im, 0, cota_sup1, -cota_sup2A, cota_sup2A)

    Int_B_function_re = lambda w,u: np.real((J0(w,u) + atan(w)*J2(w,u))*u*exp(u)*exp_electron(w)*(-rp(u)*alpha_z1(u) + rs(u)/alpha_z1(u)))
    Int_B_function_im = lambda w,u: np.imag((J0(w,u) + atan(w)*J2(w,u))*u*exp(u)*exp_electron(w)*(-rp(u)*alpha_z1(u) + rs(u)/alpha_z1(u)))

    intB_re,err = integrate.dblquad(Int_B_function_re, cota_inf1, cota_sup1A, -cota_sup2A, cota_sup2A)
    intB_im,err = integrate.dblquad(Int_B_function_im, cota_inf1, cota_sup1A, -cota_sup2A, cota_sup2A)
        
#########################################################################################################################################

    charge_e = 1.602*1e-19*c
    cte_aux = charge_e*omega/(2*np.pi*int_v) ## se cancela el -1j con el 1j de las integrales I0 5, I2 5, etc (*)

    cte = k1_2*0.5*cte_aux #(*)
    rta = (intB_re + intA_re + 1j*(intA_im + intB_im))*cte
    
    return rta

#%%

def Efield_ana1b(omegac,epsi1,epsi2,hbmu,hbgama,int_v,b,zp,z):
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
    dividido por 1j*e/(2*pi)
    analitico opcion 1.b
    """
    E = omegac*aux  
    omega = omegac*c
#    k0 = omegac
    n1 = epsi1*mu1
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1
    k1_2 = k1**2
    n_v1 = int_v/cte1
    
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    Rp = 2*epsi1/(epsi1 + epsi2)
    alfa_p = 1j*(epsi1 + epsi2)/(cond*cte1)

    exp_electron = np.exp(-alfa_p*k1*(np.abs(z) + np.abs(b) + 2*zp))

    charge_e = 1.602*1e-19*c
    cte_aux = charge_e*omega/int_v 
    
    cte = np.pi*k1_2*Rp*(alfa_p**(5/2))/2
    cte2 = 1/np.sqrt(alfa_p + n_v1) + 1/np.sqrt(alfa_p - n_v1) 
    cte3 = (1 - 1j)*1j*np.sqrt(-1j)
    
    final_expression = cte*cte2*cte3*cte_aux*exp_electron
    
    return final_expression

#%%