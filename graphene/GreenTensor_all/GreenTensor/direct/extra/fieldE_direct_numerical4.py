
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila

green tensor : integrales resueltas numericamente
luego de aplicar la aprox QE + sin aplicar la aprox QE

diferencia con la version2 : version simplificada del caso
sin QE (ver overleaf formula 169)
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
path_constants =  path_basic.replace('GreenTensor/' + 'direct','')
#print('Importar modulos necesarios para este codigo')

try:
    sys.path.insert(1, path_constants)
    from constants import constantes
except ModuleNotFoundError:
    print('constants.py no se encuentra en ' + path_constants)

pi,hb,c,alfac,mu1,mu2 = constantes()
aux = c*hb

#%%

def Efield_NUM_QE_2terms(omegac,epsi1,int_v,b):
    """
    Parameters
    ----------
    omegac : omega/c = k0 en 1/micrometros    
    epsi1 : epsilon del medio de arriba del plano
    int_v : fraccion de la vel de la c que vale la velocidad uniforme del electron en x, si v = c/int
    b: electron position en z, b < 0 (dipolo en 0 y eje de z apunta hacia abajo) 
    Returns
    -------
    External field direct (green tensor direct)
    numerico habiendo aplicado QE al green tensor
    y despreciando
    2 terminos (rosa y violeta)
    al momento de calcular el campo E. 
    Dipolo en el 0,0,0.
    """
    
    omega = omegac*c #=omega/c    
    n1 = epsi1*mu1
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1
    # k1_2 = k1**2
    k1_2 = k1**2

# green tensor evaluated in : DIPOLO EN EL 0,0,0 
#    x,y,z = 0, 0, 0
#    yD,zD = 0, b

    # green tensor evaluated in x,y,z = 0, 0, 0 and yD,zD = 0, b
    z_dip_barra = k1*np.abs(b)               # green tensor evaluated in x,y,z = 0, 0, 0 and yD,zD = 0, b
#    phi = 0  # green tensor evaluated in 
    
    expB = lambda u: np.exp(-u*z_dip_barra) 
    sin_electron = lambda w: np.sin(w*int_v/n1)
    cos_electron = lambda w: np.cos(w*int_v/n1)
    J0 = lambda w,u: special.jn(0,u*w) 


    cota_sup1 = 100
    cota_sup2 = 60*k1
############ I0 5 ########################################################
    Int05_B_function_re = lambda w,u: J0(w,u)*expB(u)*cos_electron(w)*(u**2 - 1)
    Int05_B_function_im = lambda w,u: J0(w,u)*expB(u)*sin_electron(w)*(u**2 - 1)

    int05B_re,err = integrate.dblquad(Int05_B_function_re, 0, cota_sup1, -cota_sup2, cota_sup2)
    int05B_im,err = integrate.dblquad(Int05_B_function_im, 0, cota_sup1, -cota_sup2, cota_sup2)
#########################################################################

    charge_e = 1.602*1e-19*c
    cte_aux = 1j*(charge_e*omega)/(2*np.pi*int_v)

    cte = k1_2*0.5*cte_aux #signo menos
    rta05 = (int05B_re + 1j*int05B_im)*cte 
    
    return rta05

#%%

def Efield_NUM_2terms(omegac,epsi1,int_v,b):
    """
    Parameters
    ----------
    omegac : omega/c = k0 en 1/micrometros    
    epsi1 : epsilon del medio de arriba del plano
    int_v : fraccion de la vel de la c que vale la velocidad uniforme del electron en x, si v = c/int
    b: electron position en z, b < 0 (dipolo en 0 y eje de z apunta hacia abajo) 
    Returns
    -------
    External field direct (green tensor direct)
    dividido por 1j*e/(2*pi)
    numerico sin haber aplicado QE al green tensor
    y despreciando
    2 terminos (rosa y violeta)
    al momento de calcular el campo E
    """

    omega = omegac*c #=omega/c    
    n1 = epsi1*mu1
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1
    k1_2 = k1**2

# green tensor evaluated in : DIPOLO EN EL 0,0,0 
#    x,y,z = 0, 0, 0
#    yD,zD = 0, b

     # green tensor evaluated in x,y,z = 0, 0, 0 and yD,zD = 0, b
    z_dip_barra = k1*np.abs(b)               # green tensor evaluated in x,y,z = 0, 0, 0 and yD,zD = 0, b
#    phi = 0  # green tensor evaluated in 

    alpha_z1 = lambda u: np.sqrt(1 - u**2)
    alpha_z2 = lambda u: np.sqrt(u**2 - 1)
    
    sinA = lambda u: np.sin(alpha_z1(u)*z_dip_barra) 
    cosA = lambda u: np.cos(alpha_z1(u)*z_dip_barra) 
    
    expB = lambda u: np.exp(-alpha_z2(u)*z_dip_barra) # no esta el 1j en alfa2
    
    sin_electron = lambda w: np.sin(w*int_v/n1)
    cos_electron = lambda w: np.cos(w*int_v/n1)
    
 #   exp_electron = lambda w: np.exp(1j*omegac*w*int_v)
    
    J0 = lambda w,u: special.jn(0,u*w) 
    
    cota_sup1 = 0.8
    cota_inf1 = 1.2
    
    cota_sup1A = 100
    cota_sup2A = 100*k1
############ parte A ########################################################
    Int05_A_function_re = lambda w,u: J0(w,u)*u*(cosA(u)*cos_electron(w) - sinA(u)*sin_electron(w))*(alpha_z1(u) + 1/alpha_z1(u))
#    Int05_A_function_im = lambda w,u: np.imag(J0(w,u)*u*expA(u)*exp_electron(w)*(alpha_z1(u) + 1/alpha_z1(u)))
    # esta parte imaginaria es cero, ver overleaf

    int05A_re,err = integrate.dblquad(Int05_A_function_re, 0, cota_sup1, 0, cota_sup2A)
#    int05A_im,err = integrate.dblquad(Int05_A_function_im, 0, cota_sup1, lambda u: 0, lambda u:cota_sup2)

###################### parte B ###############################################

    Int05_B_function_re = lambda w,u: -J0(w,u)*u*expB(u)*sin_electron(w)*(alpha_z2(u) - 1/alpha_z2(u)) #hay un 1j afuera, por eso la parte real es con sin y con un menos

    Int05_B_function_im = lambda w,u: J0(w,u)*u*expB(u)*cos_electron(w)*(alpha_z2(u) - 1/alpha_z2(u))
                                                                    # saque el 1j afuera asi que hay
                                                                    # un menos en 1/alfa2
                                                                                                            
    int05B_re,err = integrate.dblquad(Int05_B_function_re, cota_inf1, cota_sup1A, -cota_sup2A, cota_sup2A)
    int05B_im,err = integrate.dblquad(Int05_B_function_im, cota_inf1, cota_sup1A, -cota_sup2A, cota_sup2A)
    
##########################################################################

    charge_e = 1.602*1e-19*c
    cte_aux = (charge_e*omega)/(2*np.pi*int_v) ## se cancela el i con el de la cte y se cancela el menos

    cte = k1_2*cte_aux
    rta05 = (int05A_re + (int05B_re + 1j*int05B_im)*0.5)*cte 
    
    return rta05

#%%