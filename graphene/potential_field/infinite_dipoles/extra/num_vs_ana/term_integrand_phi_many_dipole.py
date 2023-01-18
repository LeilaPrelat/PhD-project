#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila

green tensor : integrales resueltas numericamente
luego de aplicar la aprox QE + sin aplicar la aprox QE
"""

import numpy as np
import sys
import os 

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/potential_field/infinite_dipoles/num_vs_ana','')
#print('Importar modulos necesarios para este codigo')

try:
    sys.path.insert(1, path_constants)
    from graphene_sigma import sigma_DL
except ModuleNotFoundError:
    print('graphene_sigma.py no se encuentra en ' + path_constants)

try:
    sys.path.insert(1, path_constants)
    from constants import constantes
except ModuleNotFoundError:
    print('constants.py no se encuentra en ' + path_basic)

pi,hb,c,alfac,mu1,mu2 = constantes()
aux = c*hb

#%%

def phi_many_dipoles_integrand_term1(omegac,epsi1,epsi2,hbmu,hbgama,x,y,a,N,zp,int_v,px,alpha_y):
    """
    Parameters
    ----------
    omegac : omega/c = k0 en 1/micrometros    
    epsi1 : epsilon del medio de arriba del plano
    epsi2 : epsilon del medio de abajo del plano
    hbmu : chemical potential in eV  
    hbgama : collision frequency in eV
    x : coordenada x 
    y : coordenada y 
    z : coordenada z
    a : distancia en x entre dipolos
    N : cantidad de dipolos
    zp : coordenada zp del plano
    Returns
    -------
    potential creado por muchos dipolos 
    
    """
    E = omegac*aux

    # x_y = ky/k0
    n1 = epsi1*mu1
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1
    
    cond = 4*pi*alfac*sigma_DL(E,hbmu,hbgama)
    # Rp = 2*epsi1/(epsi1 + epsi2)
    # alfa_p = 1j*(epsi1 + epsi2)/(cond*cte1)
    # kp = alfa_p*k1     

    int_tot0 = 0 
#    print('hola')
    list_dipolos = np.linspace(-N,N,2*N+1)
    
    factor = a*(2*np.pi)**2
    
#    print(2*pi/(omegac*a))
    for n in list_dipolos:
#        print(n)
        alpha_x = int_v + 2*pi*n/(omegac*a)
        
        
        
        alpha_parallel = np.sqrt(alpha_x**2 + alpha_y**2)
        # exp_z1 = lambda ky: np.exp(-k_parallel(ky)*np.abs(z))
        exp_z = np.exp(-alpha_parallel*2*k1*np.abs(zp))
        # rp = lambda ky: Rp*kp/(k_parallel(ky) - kp)        
#        exp_xy = np.exp(1j*k1*(alpha_y*y + alpha_x*x))

        exp_xy = np.exp(1j*k1*(alpha_y*y + alpha_x*x))
        
        rp_num = epsi2*1j*alpha_parallel - epsi1*1j*alpha_parallel - cond*(alpha_parallel**2)
        rp_den = epsi2*1j*alpha_parallel + epsi1*1j*alpha_parallel - cond*(alpha_parallel**2)
        rp = rp_num/rp_den
        
    ################### I1 ################################################
        
    ################### I1 ################################################
#        Int1 = 1j*px*k1*alpha_x*exp_z*exp_xy*rp/alpha_parallel
        Int1 = alpha_x*exp_z*exp_xy*rp/alpha_parallel    
    ################### I3 ###############################################

    ################### I6 ###############################################    
        
        int_tot0 = int_tot0 + k1*Int1
#        print(n,int_tot0)
#        print(int_tot0)    
    return int_tot0/factor

#%%

def phi_many_dipoles_integrand_term2(omegac,epsi1,epsi2,hbmu,hbgama,x,y,a,N,zp,int_v,py,alpha_y):
    """
    Parameters
    ----------
    omegac : omega/c = k0 en 1/micrometros    
    epsi1 : epsilon del medio de arriba del plano
    epsi2 : epsilon del medio de abajo del plano
    hbmu : chemical potential in eV  
    hbgama : collision frequency in eV
    x : coordenada x 
    y : coordenada y 
    z : coordenada z
    a : distancia en x entre dipolos
    N : cantidad de dipolos
    zp : coordenada zp del plano
    Returns
    -------
    potential creado por muchos dipolos 
    
    """
    E = omegac*aux

    # x_y = ky/k0
    n1 = epsi1*mu1
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1
    
    cond = 4*pi*alfac*sigma_DL(E,hbmu,hbgama)
    # Rp = 2*epsi1/(epsi1 + epsi2)
    # alfa_p = 1j*(epsi1 + epsi2)/(cond*cte1)
    # kp = alfa_p*k1    
    
    int_tot0 = 0 
#    print('hola')
    list_dipolos = np.linspace(-N,N,2*N+1)
    
    factor = a*(2*np.pi)**2
    for n in list_dipolos:
#        print(n)
        alpha_x = int_v + 2*pi*n/(omegac*a)
        
        alpha_parallel =  np.sqrt(alpha_x**2 + alpha_y**2)
        # exp_z1 = lambda ky: np.exp(-k_parallel(ky)*np.abs(z))
        exp_z = np.exp(-alpha_parallel*2*k1*np.abs(zp))
        
#        exp_xy = np.exp(1j*k1*(alpha_y*y + alpha_x*x))

        exp_xy = np.exp(1j*k1*(alpha_y*y + alpha_x*x))

        rp_num =  epsi2*1j*alpha_parallel - epsi1*1j*alpha_parallel - cond*(alpha_parallel**2)
        rp_den =  epsi2*1j*alpha_parallel + epsi1*1j*alpha_parallel - cond*(alpha_parallel**2)
        rp =  rp_num/rp_den
        
    ################### I1 ################################################
        
    ################### I1 ################################################
#        Int1 =  1j*py*k1*alpha_y*exp_z*exp_xy*rp/alpha_parallel
        Int1 =  alpha_y*exp_z*exp_xy*rp/alpha_parallel    
    
    ################### I3 ###############################################

    ################### I6 ###############################################    
        
        int_tot0 = int_tot0 + Int1*k1
#        print(int_tot0)    
    return int_tot0/factor

#%%


def phi_many_dipoles_integrand_term3(omegac,epsi1,epsi2,hbmu,hbgama,x,y,a,N,zp,int_v,pz,alpha_y):
    """
    Parameters
    ----------
    omegac : omega/c = k0 en 1/micrometros    
    epsi1 : epsilon del medio de arriba del plano
    epsi2 : epsilon del medio de abajo del plano
    hbmu : chemical potential in eV  
    hbgama : collision frequency in eV
    x : coordenada x 
    y : coordenada y 
    z : coordenada z
    a : distancia en x entre dipolos
    N : cantidad de dipolos
    zp : coordenada zp del plano
    Returns
    -------
    potential creado por muchos dipolos 
    
    """
    E = omegac*aux

    # x_y = ky/k0
    n1 = epsi1*mu1
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1
    
    cond = 4*pi*alfac*sigma_DL(E,hbmu,hbgama)
    # Rp = 2*epsi1/(epsi1 + epsi2)
    # alfa_p = 1j*(epsi1 + epsi2)/(cond*cte1)
    # kp = alfa_p*k1    

    int_tot0 = 0 
#    print('hola')
    list_dipolos = np.linspace(-N,N,2*N+1)
    
    factor = a*(2*np.pi)**2
    for n in list_dipolos:
#        print(n)
        alpha_x = int_v + 2*pi*n/(omegac*a)
        
        alpha_parallel =  np.sqrt(alpha_x**2 + alpha_y**2)
        # exp_z1 = lambda ky: np.exp(-k_parallel(ky)*np.abs(z))
        exp_z = np.exp(-alpha_parallel*2*k1*np.abs(zp))
        
#        exp_xy = np.exp(1j*k1*(alpha_y*y + alpha_x*x))

        exp_xy = np.exp(1j*k1*(alpha_y*y + alpha_x*x))

        rp_num = epsi2*1j*alpha_parallel - epsi1*1j*alpha_parallel - cond*(alpha_parallel**2)
        rp_den =  epsi2*1j*alpha_parallel + epsi1*1j*alpha_parallel - cond*(alpha_parallel**2)
        rp =  rp_num/rp_den
        
    ################### I1 ################################################
        
    ################### I1 ################################################
        Int1 =  -exp_z*exp_xy*rp
    
    
    ################### I3 ###############################################

    ################### I6 ###############################################    
        
        int_tot0 = int_tot0 + Int1*k1
#        print(int_tot0)    
    return int_tot0/factor



#%%
