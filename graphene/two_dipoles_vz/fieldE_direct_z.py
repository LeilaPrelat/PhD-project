
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
path_constants =  path_basic.replace('/two_dipoles_vx','')
#print('Importar modulos necesarios para este codigo')

try:
    sys.path.insert(1, path_constants)
    from constants import constantes
except ModuleNotFoundError:
    print('constants.py no se encuentra en ' + path_constants)

pi,hb,c,alfac,mu1,mu2 = constantes()
aux = c*hb

#%%

def Efield_NUM_QE(omegac,epsi1,R,z,int_v):
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


#    phi = np.arctan2(np.abs(y-ye),np.abs(x-xe)


    den = lambda z_monio: z_monio**2 + (omegac*R)**2

    
    function = lambda z_monio: 3*(z_monio**2)*(den(z_monio)**(-5/2)) - den(z_monio)**(-3/2) 


    cos_electron = lambda z_monio: np.cos(z_monio*int_v*omegac)
    
#
    cota_inf1 =  0.01*k1
    cota_sup2 = 1000*k1
    
############ I0 5 ########################################################
    Int05_B_function_re = lambda z_monio: np.real(function(z_monio)*cos_electron(z_monio))
    Int05_B_function_im = lambda z_monio: np.imag(function(z_monio)*cos_electron(z_monio))

    int05B_re,err = integrate.quad(Int05_B_function_re, cota_inf1, cota_sup2)
    int05B_im,err = integrate.quad(Int05_B_function_im, cota_inf1, cota_sup2)
#########################################################################

#    cte = k1_2*cte_aux #signo menos

    signo_vz = np.sign(int_v)
    rta05 = (int05B_re + 1j*int05B_im) 
    
    return -1j*signo_vz*rta05*np.exp(1j*omegac*z*int_v)*(omegac**2)

#%%

def Efield_ANA(omegac,epsi1,R,z,int_v):
    """
    Parameters
    ----------
    omegac : omega/c = k0 en 1/micrometros    
    epsi1 : epsilon del medio de arriba del plano
    int_v : fraccion de la vel de la c que vale la velocidad uniforme del electron en x, si v = c/int
    b : electron position en z, b < 0 (dipolo en 0 y eje de z apunta hacia abajo)  
    Returns
    -------
    External field direct del paper 149
    """
    arg = np.abs(R*omegac*int_v)
    K0 = special.kn(0,arg)
    
    signo_vz = np.sign(int_v)
    
    cte_aux = 1j*(omegac*int_v)**2

 #   rta_final = K0*cte_v*cte_aux 
    
    return signo_vz*K0*cte_aux*np.exp(1j*omegac*z*int_v)


#%%


