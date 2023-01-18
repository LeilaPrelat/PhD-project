
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
path_constants =  path_basic.replace('/External_Efield','')
#print('Importar modulos necesarios para este codigo')

try:
    sys.path.insert(1, path_constants)
    from constants import constantes
except ModuleNotFoundError:
    print('constants.py no se encuentra en ' + path_constants)

pi,hb,c,alfac,mu1,mu2 = constantes()
aux = c*hb

#%%

def Efield_NUM_QE(omegac,epsi1,int_v,b,x,y,z):
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

    # green tensor evaluated in x,y,z = 0, 0, 0 and yD,zD = 0, b
    z_dip_barra = k1*np.abs(z-b)               # green tensor evaluated in x,y,z = 0, 0, 0 and yD,zD = 0, b
    atan = lambda w : np.cos(2*np.arctan2(np.abs(y_bar),np.abs(x_bar - w))) 
    
    expB = lambda u: np.exp(-u*z_dip_barra) 
    
    sin_electron = lambda w: np.sin(w*int_v/n1)
    cos_electron = lambda w: np.cos(w*int_v/n1)
    
    
    J0 = lambda w,u: special.jn(0,u*np.sqrt((x_bar-w)**2 + y_bar**2)) 
    J2 = lambda w,u: special.jn(2,u*np.sqrt((x_bar-w)**2 + y_bar**2)) 

    cota_sup1 = 3001
    cota_sup2 = 80*k1
############ I0 5 ########################################################
    Int05_B_function_re = lambda w,u: expB(u)*cos_electron(w)*(J0(w,u) + atan(w)*J2(w,u))*(u**2 - 1)
    Int05_B_function_im = lambda w,u: expB(u)*sin_electron(w)*(J0(w,u) + atan(w)*J2(w,u))*(u**2 - 1)

    int05B_re,err = integrate.dblquad(Int05_B_function_re, 1, cota_sup1, -cota_sup2, cota_sup2)
    int05B_im,err = integrate.dblquad(Int05_B_function_im, 1, cota_sup1, -cota_sup2, cota_sup2)
#########################################################################

    charge_e = 1.602*1e-19*c
    cte_aux = 1j*(charge_e*omega)/(2*np.pi*int_v)

    cte = k1_2*0.5*cte_aux #signo menos
    rta05 = (int05B_re + 1j*int05B_im)*cte 
    
    return rta05

#%%

def Efield_NUM(omegac,epsi1,int_v,b,x,y,z):
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
    numerico sin haber aplicado QE al green tensor
    al momento de calcular el campo E
    """

#    omega = omegac*c #=omega/c    
    n1 = epsi1*mu1
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1
    k1_2 = k1**2
    x_bar = x*k1
    y_bar = y*k1

# green tensor evaluated in : DIPOLO EN EL 0,0,0 
#    x,y,z = 0, 0, 0
#    yD,zD = 0, b

     # green tensor evaluated in x,y,z = 0, 0, 0 and yD,zD = 0, b
    z_dip_barra = k1*np.abs(z - b)               # green tensor evaluated in x,y,z = 0, 0, 0 and yD,zD = 0, b
    atan = lambda w : np.cos(2*np.arctan2(np.abs(y_bar),np.abs(x_bar - w)))

    alpha_z1 = lambda u: np.sqrt(1 - u**2)
    alpha_z2 = lambda u: np.sqrt(u**2 - 1)
    
    expA = lambda u: np.exp(1j*alpha_z1(u)*z_dip_barra) 
    expB = lambda u: np.exp(-alpha_z2(u)*z_dip_barra)
     
    exp_electron = lambda w: np.exp(1j*w*int_v/n1)
    J0 = lambda w,u: special.jn(0,u*np.sqrt((x_bar-w)**2 + y_bar**2)) 
    J2 = lambda w,u: special.jn(2,u*np.sqrt((x_bar-w)**2 + y_bar**2)) 
    
    cota_sup1 = 0.95
    cota_inf1 = 1.05
#    cota_sup2 = 1e5
    cota_sup1A = 3000
    cota_sup2A = 80*k1

############ I0 5 ########################################################
    Int05_A_function_re = lambda w,u: np.real((J0(w,u) + atan(w)*J2(w,u))*u*expA(u)*exp_electron(w)*(alpha_z1(u) + 1/alpha_z1(u)))
    Int05_A_function_im = lambda w,u: np.imag((J0(w,u) + atan(w)*J2(w,u))*u*expA(u)*exp_electron(w)*(alpha_z1(u) + 1/alpha_z1(u)))

    int05A_re,err = integrate.dblquad(Int05_A_function_re, 0, cota_sup1, -cota_sup2A, cota_sup2A)
    int05A_im,err = integrate.dblquad(Int05_A_function_im, 0, cota_sup1, -cota_sup2A, cota_sup2A)

    Int05_B_function_re = lambda w,u: np.real(1j*(J0(w,u) + J2(w,u))*u*expB(u)*exp_electron(w)*(alpha_z2(u) - 1/alpha_z2(u)))
    Int05_B_function_im = lambda w,u: np.imag(1j*(J0(w,u) + J2(w,u))*u*expB(u)*exp_electron(w)*(alpha_z2(u) - 1/alpha_z2(u)))

    int05B_re,err = integrate.dblquad(Int05_B_function_re, cota_inf1, cota_sup1A, -cota_sup2A, cota_sup2A)
    int05B_im,err = integrate.dblquad(Int05_B_function_im, cota_inf1, cota_sup1A, -cota_sup2A, cota_sup2A)
    
##########################################################################

    charge_e = 1.602*1e-19
    cte_aux = (charge_e*c*omegac*c)/(2*np.pi*int_v) ## se cancela el i con el de la cte y se cancela el menos

    cte = k1_2*0.5*cte_aux
    rta05 = (int05B_re + int05A_re + 1j*(int05B_im + int05A_im))*cte 
    
    return rta05

#%%