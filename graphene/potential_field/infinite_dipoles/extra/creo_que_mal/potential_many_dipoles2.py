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

def phi_many_dipoles_num2(omegac,epsi1,epsi2,hbmu,hbgama,x,y,z,a,N,zp,int_v,px,py,pz):
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

    int_tot0 = 0 
#    print('hola')
    list_dipolos = np.linspace(-N,N,2*N+1)
    for n in list_dipolos:
#        print(n)
        kx = omegac*c/v - 2*pi*n/a
        
        k_parallel = lambda ky: np.sqrt(kx**2 + ky**2)
        exp_z1 = lambda ky: np.exp(-k_parallel(ky)*np.abs(z))
        exp_z2 = lambda ky: np.exp(-k_parallel(ky)*(2*np.abs(zp) + np.abs(z)))
        
        exp_xy = lambda ky: np.exp(1j*ky*y + 1j*kx*x)
        rp = lambda ky: Rp*kp/(k_parallel(ky) - kp)
    ################### I1 ################################################
        Int1_B_function_re = lambda ky: np.real(-1j*(px*kx + ky*py)*exp_z1(ky)*exp_xy(ky)/k_parallel(ky))
        Int1_B_function_im = lambda ky: np.imag(-1j*(px*kx + ky*py)*exp_z1(ky)*exp_xy(ky)/k_parallel(ky))

        int1_re,err = integrate.quad(Int1_B_function_re, cota_inf, cota_sup)
        int1_im,err = integrate.quad(Int1_B_function_im, cota_inf, cota_sup)
    
    ################### I2 ###############################################
        Int2_B_function_re = lambda ky: np.real(pz*np.sign(z)*exp_z1(ky)*exp_xy(ky))
        Int2_B_function_im = lambda ky: np.imag(pz*np.sign(z)*exp_z1(ky)*exp_xy(ky))
    
        int2_re,err = integrate.quad(Int2_B_function_re, cota_inf, cota_sup)
        int2_im,err = integrate.quad(Int2_B_function_im, cota_inf, cota_sup)
    ################### I3 ###############################################
        Int3_B_function_re = lambda ky: np.real(px*kx*1j*rp(ky)*exp_xy(ky)*exp_z2(ky)/k_parallel(ky))
        Int3_B_function_im = lambda ky: np.imag(px*kx*1j*rp(ky)*exp_xy(ky)*exp_z2(ky)/k_parallel(ky))
    
        int3_re,err = integrate.quad(Int3_B_function_re, cota_inf, cota_sup)
        int3_im,err = integrate.quad(Int3_B_function_im, cota_inf, cota_sup)
    ################### I4 ###############################################
        Int4_B_function_re = lambda ky: np.real(py*ky*1j*rp(ky)*exp_xy(ky)*exp_z2(ky)/k_parallel(ky))
        Int4_B_function_im = lambda ky: np.imag(py*ky*1j*rp(ky)*exp_xy(ky)*exp_z2(ky)/k_parallel(ky))
#    
        int4_re,err = integrate.quad(Int4_B_function_re, cota_inf, cota_sup)
        int4_im,err = integrate.quad(Int4_B_function_im, cota_inf, cota_sup)
    ################### I5 ###############################################
        Int5_B_function_re = lambda ky: np.real(-pz*np.sign(z)*rp(ky)*exp_xy(ky)*exp_z2(ky))
        Int5_B_function_im = lambda ky: np.imag(-pz*np.sign(z)*rp(ky)*exp_xy(ky)*exp_z2(ky))
    
        int5_re,err = integrate.quad(Int5_B_function_re, cota_inf, cota_sup)
        int5_im,err = integrate.quad(Int5_B_function_im, cota_inf, cota_sup)

    ################### I6 ###############################################    
        
        int_re_tot = int1_re + int2_re + int3_re + int4_re + int5_re
        int_im_tot = int1_im + int2_im + int3_im + int4_im + int5_im 
        
        int_tot0 = int_tot0 + a*(int_re_tot + 1j*int_im_tot)
#        print(int_tot0)    
    return int_tot0

#%%
    
def phi_many_dipoles_ana2(omegac,epsi1,epsi2,hbmu,hbgama,x,y,z,a,N,zp,int_v,px,py,pz):
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
#    kp = alfa_p*k1    
    
    v = c/int_v
    term_final = 0
    list_dipolos = np.linspace(-N,N,2*N+1)
    for n in list_dipolos:
#        print(n)
        kx = omegac*c/v - 2*pi*n/a
        alfa_x = kx/k1
        
        alfa_y = np.sqrt(alfa_p**2 - alfa_x**2)
        
        exp_z1 = np.exp(-2*alfa_p*k1*zp)
        exp_z2 = np.exp(-alfa_p*k1*np.abs(z))
        
        exp_xy = np.exp(1j*alfa_y*k1*y + 1j*alfa_x*k1*x)
        
        term1 = -1j*(px*alfa_x + py*alfa_y)/alfa_p
        term2 = pz*np.sign(z)
        term3 = 1j*Rp*(px*alfa_x + py*alfa_y)*exp_z1
        term4 = -pz*alfa_p*np.sign(z)*Rp*exp_z1
        
        
        factor_final = exp_xy*exp_z2*a*k1
        
        ################### I1 ################################################
        
        term_final = (term1 + term2 + term3 + term4)*factor_final + term_final 
        
        
    return term_final

#%%

