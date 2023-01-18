#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila
"""
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
    from graphene_sigma import sigma_DL
except ModuleNotFoundError:
    print('graphene_sigma.py no se encuentra en ' + path_basic)
    

try:
    sys.path.insert(1, path_basic)
    from dipole_moment import dipole_moment_x_num,dipole_moment_x_ana
except ModuleNotFoundError:
    print('dipole_moment.py no se encuentra en ' + path_basic)




try:
    sys.path.insert(1, path_basic)
    from fieldE_ref import Efield_ref_ana
except ModuleNotFoundError:
    print('fieldE_ref.py no se encuentra en ' + path_basic)
    


try:
    sys.path.insert(1, path_basic)
    from fieldE_direct_x import Efield_ANA
except ModuleNotFoundError:
    print('fieldEdir_direct_numerical.py no se encuentra en ' + path_basic)
      


try:
    sys.path.insert(1, path_constants)
    from constants import constantes
except ModuleNotFoundError:
    print('constants.py no se encuentra en ' + path_constants)

pi,hb,c,alfac,mu1,mu2 = constantes()
aux = c*hb


#%%

def EELS_parallel_f_dir(omegac,epsi1,epsi2,hbmu,hbgama,int_v,b,zp,L):     
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
    E = omegac*aux    

    n1 = epsi1*mu1
    cte1 = np.sqrt(n1)
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    Rp = 2*epsi1/(epsi1 + epsi2)
    alfa_p = 1j*(epsi1 + epsi2)/(cond*cte1)
    kp = alfa_p*omegac
    kp_2 = kp**2
    
    expo = np.exp(-2*kp*(zp-b))
    
#    omega = omegac*c    
    rta =  Rp*kp_2*expo/np.sqrt(kp_2 - (omegac*int_v)**2)
    
    
    v = c/int_v
    cte_final = alfac*c*L/(v**2) ## adimensional

    
    return 4*np.real(rta)*cte_final

#%%

#%%

def EELS_1_ana(omegac,epsi1,epsi2,hbmu,hbgama,int_v,b,zp,x1,z1,x2,z2,L,omega01,kappa_factor_omega01,kappa_r_factor1,omega02,kappa_factor_omega02,kappa_r_factor2):     
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
    E = omegac*aux

    
    p1 = dipole_moment_x_ana(omegac,epsi1,epsi2,hbmu,hbgama,int_v,b,zp,x1,z1,x2,z2,omega01,kappa_factor_omega01,kappa_r_factor1,omega02,kappa_factor_omega02,kappa_r_factor2)
    
    Eref1 = Efield_ref_ana(omegac,epsi1,epsi2,hbmu,hbgama,x1,z1,b,zp,int_v)
    Edir1 = Efield_ANA(omegac,epsi1,x1,z1,b,int_v)

    rta = -1j*p1*(Eref1 + Edir1)
    
    
    rta2 = EELS_parallel_f_dir(omegac,epsi1,epsi2,hbmu,hbgama,int_v,b,zp,L)
    
    cte_final = alfac*c ## aparece e^2/hbar  (aparece e dentro de los campos)
    
    
    cte_que_no_puse_en_los_campos = (omegac*c)**2 # e/omega 
#    print('1:', 1j*p1)
#    print('2:', (Eref1 + Edir1))
#    
    expression_on_overleaf = rta*cte_final/(np.pi*cte_que_no_puse_en_los_campos) #queda en unidad de segundos
    
    #queremos en unidad de energia 
    
    final = expression_on_overleaf*omegac*c
    
    
    return rta*alfac/(np.pi*omegac) + rta2

#%%

def EELS_1_num(omegac,epsi1,epsi2,hbmu,hbgama,int_v,b,zp,x1,z1,x2,z2,L,omega01,kappa_factor_omega01,kappa_r_factor1,omega02,kappa_factor_omega02,kappa_r_factor2):     
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

    E = omegac*aux    

    p1 = dipole_moment_x_num(omegac,epsi1,epsi2,hbmu,hbgama,int_v,b,zp,x1,z1,x2,z2,omega01,kappa_factor_omega01,kappa_r_factor1,omega02,kappa_factor_omega02,kappa_r_factor2)
    
    Eref1 = Efield_ref_ana(omegac,epsi1,epsi2,hbmu,hbgama,x1,z1,b,zp,int_v)
    Edir1 = Efield_ANA(omegac,epsi1,x1,z1,b,int_v)


    rta = -1j*p1*(Eref1 + Edir1)
 
    rta2 = EELS_parallel_f_dir(omegac,epsi1,epsi2,hbmu,hbgama,int_v,b,zp,L)    
    
    cte_final = alfac*c ## aparece e^2/hbar  (aparece e dentro de los campos)
    
    
    cte_que_no_puse_en_los_campos = (omegac*c)**2
#    print('1:', 1j*p1)
#    print('2:', (Eref1 + Edir1))
#
    expression_on_overleaf = rta*cte_final/(np.pi*cte_que_no_puse_en_los_campos)

    final = expression_on_overleaf*omegac*c
    
    return rta*alfac/(np.pi*omegac) + rta2


#%%
    

