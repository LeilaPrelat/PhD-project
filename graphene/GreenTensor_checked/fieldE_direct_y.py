
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
    sys.path.insert(1, path_constants)
    from constants import constantes
except ModuleNotFoundError:
    print('constants.py no se encuentra en ' + path_constants)

pi,hb,c,alfac,mu1,mu2 = constantes()
aux = c*hb

#%%

def Efield_NUM_QE(omegac,epsi1,x,z,ze,int_v):
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

# green tensor evaluated in : DIPOLO EN EL 0,0,0 
#    x,y,z = 0, 0, 0
#    yD,zD = 0, b


#    phi = np.arctan2(np.abs(y-ye),np.abs(x-xe))
    
    
    cte_all = k1*np.abs(z-ze)**2
    
    function = lambda xe_monio: ( cte_all + xe_monio**2 )**(-3/2)


    cos_electron = lambda xe_monio: np.cos(xe_monio*int_v)
    
#
#    cota_sup1 = -400

    signo_vy = np.sign(int_v)
    
    cota_inf1 = 0.01
    cota_sup2 = 1000
    
############ I0 5 ########################################################
    Int05_B_function_re = lambda xe_monio: np.real(function(xe_monio)*cos_electron(xe_monio))
    Int05_B_function_im = lambda xe_monio: np.imag(function(xe_monio)*cos_electron(xe_monio))

    int05B_re,err = integrate.quad(Int05_B_function_re, cota_inf1, cota_sup2)
    int05B_im,err = integrate.quad(Int05_B_function_im, cota_inf1, cota_sup2)
#########################################################################

#    cte = k1_2*cte_aux #signo menos
    rta05 = (int05B_re + 1j*int05B_im) 
    
    return -rta05*1j*2*k1_2*signo_vy*np.exp(1j*omegac*x*int_v)

#%%

def Efield_ANA(omegac,epsi1,x,z,ze,int_v):
    """
    Parameters
    ----------
    omegac : omega/c = k0 en 1/micrometros    
    epsi1 : epsilon del medio de arriba del plano
    int_v : fraccion de la vel de la c que vale la velocidad uniforme del electron en x, si v = c/int
    b : electron position en z, b < 0 (dipolo en 0 y eje de z apunta hacia abajo)  
    Returns
    -------
    External field direct componente sin el e/omega
    """
    
    # x_y = ky/k0
#    n1 = epsi1*mu1
#    cte1 = np.sqrt(n1)
#    k1 = omegac*cte1
#    k1_2 = k1**2
    
    arg = np.abs(ze - z)*omegac*int_v
    K1 = special.kn(1,arg)
#    K1 = special.kn(1,arg)
#    K2 = special.kn(2,arg)
    # cte = omega/(v*np.abs(b))
    
 #   v = c/int_v
    
#    cte_v = 1 - (1/int_v)**2
#
#    cte_aux = 1j*2*(omegac*int_v)**2

    rta_final = K1*omegac*int_v/np.abs(ze - z)
    
    signo_vy = np.sign(int_v)
    
    return -1j*2*signo_vy*rta_final*np.exp(1j*omegac*x*int_v)





