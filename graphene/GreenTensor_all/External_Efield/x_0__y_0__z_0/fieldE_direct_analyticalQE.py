#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila

campo externo directo analitico con la convencion de z hacia abajo
en el punto x  = 0, y = 0,z = 0
"""
import numpy as np
import sys
import os 
from scipy import special

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/External_Efield/x_0__y_0__z_0','')
#print('Importar modulos necesarios para este codigo')

try:
    sys.path.insert(1, path_constants)
    from constants import constantes
except ModuleNotFoundError:
    print('constants.py no se encuentra en ' + path_constants)

pi,hb,c,alfac,mu1,mu2 = constantes()
aux = c*hb

#%%

def Efield_ANA_QE_2terms(omegac,epsi1,int_v,b):
    """
    Parameters
    ----------
    omegac : omega/c = k0 en 1/micrometros    
    epsi1 : epsilon del medio de arriba del plano
    int_v : fraccion de la vel de la c que vale la velocidad uniforme del electron en x, si v = c/int
    b : electron position en z, b < 0 (dipolo en 0 y eje de z apunta hacia abajo)  
    Returns
    -------
    External field direct (green tensor direct)
    analitico luego de aplicar QE y de despreciar
    2 terminos (rosa y violeta)
    """
    
    # x_y = ky/k0
    n1 = epsi1*mu1
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1
    k1_2 = k1**2
    
    arg = np.abs(b)*omegac*int_v
    K0 = special.kn(0,arg)
    K1 = special.kn(1,arg)
    K2 = special.kn(2,arg)
    # cte = omega/(v*np.abs(b))
    
    omega = omegac*c
    v = c/int_v

    charge_e = 1.602*1e-19
    cte_aux = 1j*(omega*charge_e)/(2*np.pi)

    rta05 = - k1_2*v*K0
    rta06 = omega*0.5*(omegac*K2*int_v - K1/np.abs(b)) 
    
    return (rta05 + rta06)*cte_aux

#%%


def Efield_ANA_paper149(omegac,epsi1,int_v,b):
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
    
    # x_y = ky/k0
    n1 = epsi1*mu1
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1
    k1_2 = k1**2
    
    arg = np.abs(b)*omegac*int_v
    K0 = special.kn(0,arg)
    K1 = special.kn(1,arg)
    K2 = special.kn(2,arg)
    # cte = omega/(v*np.abs(b))
    
    omega = omegac*c
    v = c/int_v

    charge_e = 1.602*1e-19
    cte_aux = 1j*(omega*charge_e)/(2*np.pi)

    rta05 = - k1_2*v*K0
    rta06 = omega*0.5*(omegac*K2*int_v - K1/np.abs(b)) 
    
    return (rta05 + rta06)*cte_aux

#%%


