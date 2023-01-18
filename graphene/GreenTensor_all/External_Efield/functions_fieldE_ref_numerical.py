
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila

funciones dentro de
el campo externo reflejado 
 
"""
from scipy import special
import numpy as np
import sys
import os 

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/External_Efield','')
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

def Efield_NUM_QE(omegac,epsi1,epsi2,hbmu,hbgama,int_v,b,zp,x,y,z,alpha_parallel,xD_tilde):
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
    x : punto x donde miramos el campo en micrones
    y : punto y donde miramos el campo en micrones
    z : punto z donde miramos el campo en micrones    
    Returns
    -------
    External field reflected (green tensor reflected)
    numerico habiendo aplicado QE al green tensor. 
    Dipolo en el 0,0,0.
    """
    E = omegac*aux    
    omega = omegac*c #=omega/c    
    n1 = epsi1*mu1
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1
    # k1_2 = k1**2
    k1_2 = k1**2
    x_bar = x*k1
    y_bar = y*k1

# green tensor evaluated in : DIPOLO EN EL 0,0,0 
#    x,y,z = 0, 0, 0
#    yD,zD = 0, b

    u,w = alpha_parallel, xD_tilde  

    # Rbarra = lambda w: k1*w # green tensor evaluated in x,y,z = 0, 0, 0 and yD,zD = 0, b
    z_dip_barra = k1*(z + np.abs(b) + 2*zp)  #NO ESTOY SEGURA DE ESTA EXPONENCIAL EH
    atan = np.cos(2*np.arctan2(np.abs(y_bar),np.abs(x_bar + w)))

    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)

    rp_num = epsi2*1j*u - epsi1*1j*u - cte1*cond*(u**2)
    rp_den = epsi2*1j*u + epsi1*1j*u - cte1*cond*(u**2)
    rp = rp_num/rp_den
    
    rs_num = 1j*u - 1j*u - cond/cte1
    rs_den = 1j*u + 1j*u + cond/cte1
    rs = rs_num/rs_den

    expB = np.exp(-u*z_dip_barra) 
    exp_electron = np.cos(w*int_v/n1) + 1j*np.sin(w*int_v/n1)
    J0 = special.jn(0,u*np.sqrt((x_bar+w)**2 + y_bar**2)) 
    J2 = special.jn(2,u*np.sqrt((x_bar+w)**2 + y_bar**2)) 

############ I0 5 + I2 5 #################################################################################################
#    Int5_B_function_re = lambda w,u: np.real((J0(w,u) + J2(w,u))*rs(u)*expB(u)*exp_electron(w))
#    Int5_B_function_im = lambda w,u: np.imag((J0(w,u) + J2(w,u))*rs(u)*expB(u)*exp_electron(w))

#    int5B_re,err = integrate.dblquad(Int5_B_function_re, 0, cota_sup1, -cota_sup2, cota_sup2)
#    int5B_im,err = integrate.dblquad(Int5_B_function_im, 0, cota_sup1, -cota_sup2, cota_sup2)
    
    ff = (J0 + atan*J2)*((u**2)*rp + rs)*expB*exp_electron
############ I0 5 + I2 5 + I0 6 + I2 6 #################################################################################################
    IntB_function_re = np.real(ff)
    IntB_function_im = np.imag(ff)

##########################################################################################################################

    charge_e = 1.602*1e-19*c
    cte_aux = 1j*(charge_e*omega)/(2*int_v)

    cte = k1_2*0.5*cte_aux #signo menos
    
#    rta5 = (int5B_re + 1j*int5B_im)*cte 
    rta = (IntB_function_re + 1j*IntB_function_im)*cte 
    
    return -rta

#%%

def Efield_NUM(omegac,epsi1,epsi2,hbmu,hbgama,int_v,b,zp,x,y,z,alpha_parallel,xD_tilde):
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
    x : punto x donde miramos el campo en micrones
    y : punto y donde miramos el campo en micrones
    z : punto z donde miramos el campo en micrones 
    Returns
    -------
    External field reflected (green tensor reflected)
    dividido por 1j*e/(2*pi)
    numerico sin haber aplicado QE al green tensor
    Dipolo en el 0,0,0. para un punto no necesariamente
    0,0,0
    """
    E = omegac*aux  
    omega = omegac*c #=omega/c    
    n1 = epsi1*mu1
    n2 = epsi2*mu2
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1
    k1_2 = k1**2
    x_bar = x*k1
    y_bar = y*k1    

# green tensor evaluated in : DIPOLO EN EL 0,0,0 
#    x,y,z = 0, 0, 0
#    yD,zD = 0, b

    u,w = alpha_parallel, xD_tilde  
#    Rbarra = lambda w: k1*w # green tensor evaluated in x,y,z = 0, 0, 0 and yD,zD = 0, b
    z_dip_barra = k1*(z + np.abs(b) + 2*zp)  #NO ESTOY SEGURA DE ESTA EXPONENCIAL EH
    atan = np.cos(2*np.arctan2(np.abs(y_bar),np.abs(x_bar + w)))
    
#    alpha_z1 = lambda u: np.sqrt(1-u**2)
#    alpha_z2 = lambda u: np.sqrt(u**2-1)
    aux2 = n2/n1
    alpha_z1 = np.sqrt(1-u**2) if u<1 else 1j*np.sqrt(u**2-1)
    alpha_z2 = np.sqrt(aux2-u**2) if u<aux2 else 1j*np.sqrt(u**2-aux2)
        
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)

    rp_num = epsi2*alpha_z1 - epsi1*alpha_z2 + cte1*cond*alpha_z1*alpha_z2
    rp_den = epsi2*alpha_z1 + epsi1*alpha_z2 + cte1*cond*alpha_z1*alpha_z2
    rp = rp_num/rp_den
    
    rs_num = alpha_z1 - alpha_z2 - cond/cte1
    rs_den = alpha_z1 + alpha_z2 + cond/cte1
    rs = rs_num/rs_den

    exp = np.exp(1j*alpha_z1*z_dip_barra) 
    exp_electron = np.cos(w*int_v/n1) + 1j*np.sin(w*int_v/n1)
    J0 = special.jn(0,u*np.sqrt((x_bar+w)**2 + y_bar**2)) 
    J2 = special.jn(2,u*np.sqrt((x_bar+w)**2 + y_bar**2)) 
    
    
    ff = (J0 + atan*J2)*u*exp*exp_electron*(rs/alpha_z1 - rp*alpha_z1)
    
############ I0 5 + I2 5  + I0 6 + I2 6 #################################################################################################
    Int_function_re = np.real(ff)
    Int_function_im = np.imag(ff)

#########################################################################################################################################

    charge_e = 1.602*1e-19*c
    cte_aux = charge_e*omega/(2*int_v) ## se cancela el -1j con el 1j de las integrales I0 5, I2 5, etc (*)

    cte = k1_2*0.5*cte_aux #(*)
    rta = (Int_function_re + 1j*Int_function_im)*cte
    
    return rta

#%%
