#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila

green tensor : integrales resueltas numericamente
luego de aplicar la aprox QE + sin aplicar la aprox QE
"""
from scipy import integrate
import numpy as np
from scipy import special
import sys
import os 

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/potential_field/many_potential_Javier_formula','')
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

def phi_many_dipoles_num3(omegac,epsi1,epsi2,hbmu,hbgama,x,y,z,a,N,zp,int_v,px,py,pz):
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
    Rp = 2*epsi1/(epsi1 + epsi2)
    alfa_p = 1j*(epsi1 + epsi2)/(cond*cte1)
    kp = alfa_p*k1    
    
    v = c/int_v
    cota_inf = 1
    cota_sup = 1000    

#    phi = np.arctan2(np.abs(y),np.abs(x))    
    R = np.sqrt(x**2 + y**2)

    int_tot0 = 0 
#    print('hola')
    list_dipolos = np.linspace(-N,N,2*N+1)
    cte_final = (2*np.pi)**2

    for n in list_dipolos:
#        print(n)
        kx = omegac*c/v - 2*pi*n/a
        
        k_parallel = lambda ky: np.sqrt(kx**2 + ky**2)

        exp_zp =  lambda ky: np.exp(-2*zp*k_parallel(ky))
#        exp_z2 = lambda ky: np.exp(-k_parallel(ky)*2*np.abs(zp))
        
        exp_xy = lambda ky: np.exp(1j*ky*y + 1j*kx*x)
        rp = lambda ky: Rp*kp/(k_parallel(ky) - kp)
    ################### I1 ################################################


        K0 = lambda ky : special.kv(0,ky*R)
        K1 = lambda ky : special.kv(1,ky*R)
        
        term1 = kx*(-1j*px*K0 + 2*py*K1) + pz*np.sign(z)*1j/y
        term1_re = np.real(term1)
        term1_im = np.imag(term1)
        
        
        Int1_B_function_re = lambda ky: np.real(1j*(px*kx + ky*py)*rp*exp_xy(ky)*exp_zp(ky)/k_parallel(ky))
        Int1_B_function_im = lambda ky: np.imag(1j*(px*kx + ky*py)*rp*exp_xy(ky)*exp_zp(ky)/k_parallel(ky))
        

        int1_re,err = integrate.quad(Int1_B_function_re, cota_inf, cota_sup)
        int1_im,err = integrate.quad(Int1_B_function_im, cota_inf, cota_sup)
    
    ################### I2 ###############################################
        Int2_B_function_re = lambda ky: np.real(-pz*rp*exp_xy(ky)*exp_zp(ky))
        Int2_B_function_im = lambda ky: np.imag(-pz*rp*exp_xy(ky)*exp_zp(ky))
    
        int2_re,err = integrate.quad(Int2_B_function_re, cota_inf, cota_sup)
        int2_im,err = integrate.quad(Int2_B_function_im, cota_inf, cota_sup)

    ################### I6 ###############################################    
        
        int_re_tot = int1_re + int2_re + term1_re
        int_im_tot = int1_im + int2_im + term1_im
        
        int_tot0 = int_tot0 + (int_re_tot + 1j*int_im_tot)/(a*cte_final)
#        print(int_tot0)    
    return int_tot0

#%%
    
def phi_many_dipoles_ana3(omegac,epsi1,epsi2,hbmu,hbgama,x,y,z,a,N,zp,int_v,px,py,pz):
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
    +
    """
    E = omegac*aux

    # x_y = ky/k0
    n1 = epsi1*mu1
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1
    
    cond = 4*pi*alfac*sigma_DL(E,hbmu,hbgama)
    Rp = 2*epsi1/(epsi1 + epsi2)
    alfa_p = 1j*(epsi1 + epsi2)/(cond*cte1)
    kp = alfa_p*k1    

    R = np.sqrt(x**2 + y**2)
    exp_zp = np.exp(-2*zp*kp)
        
    v = c/int_v
    term_final = 0
    list_dipolos = np.linspace(-N,N,2*N+1)
    cte_final = (2*np.pi)**2
    for n in list_dipolos:
#        print(n)
        kx = omegac*c/v - 2*pi*n/a
        
        if kp <  kx:
            ky = np.sqrt(kp**2 - kx**2)
        else:
            ky = 1j*np.sqrt(kx**2 - kp**2)


#        exp_z2 = lambda ky: np.exp(-k_parallel(ky)*2*np.abs(zp))
        
        exp_xy = np.exp(1j*ky*y + 1j*kx*x)
#        rp = Rp*kp/(k_parallel(ky) - kp)
    ################### I1 ################################################


        K0 = special.kv(0,ky*R)
        K1 = special.kv(1,ky*R)
        

#        term1_re = np.real(term1)
#        term1_im = np.imag(term1)
        
        term3 = -2*np.pi*(px*kx + py*ky)*Rp*exp_zp*exp_xy
        term4 = pz*np.pi*2*1j*Rp*kp*exp_zp*exp_xy
        
        print(kx*(-1j*px*K0 + 2*py*K1) + pz*np.sign(z)*1j/y)
        term1 = kx*(-1j*px*K0 + 2*py*K1) + pz*np.sign(z)*1j/y 


    ################### I6 ###############################################    
        
        term_final = term_final + (term1 + term3 + term4)/(a*cte_final)
        
        
    return term_final

#%%

