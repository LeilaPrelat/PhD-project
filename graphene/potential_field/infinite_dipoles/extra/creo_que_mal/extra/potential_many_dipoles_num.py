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
path_constants =  path_basic.replace('/potential_field/many_potential','')
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
    print('constants.py no se encuentra en ' + path_basic)

pi,hb,c,alfac,mu1,mu2 = constantes()
aux = c*hb

#%%

def phi_many_dipoles_num(omegac,epsi1,epsi2,hbmu,hbgama,x,y,z,a,N,zp,int_v,px,py,pz):
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
    
    cond = 4*pi*alfac*sigma_DL(E,hbmu,hbgama) #no hace falta dividir por c
    Rp = 2*epsi1/(epsi1 + epsi2)
    alfa_p = 1j*(epsi1 + epsi2)/(cond*cte1)
    kp = alfa_p*k1    
    
    v = c/int_v
    cota_inf = 1
    cota_sup = 800    

    int_tot0 = lambda ky: 0 
    for n in range(N):
        print(n)
        kx = omegac*c/v - 2*pi*n/a
        
        k_parallel = lambda ky: np.sqrt(kx**2 + ky**2)
        exp_z = lambda ky: np.exp(-k_parallel(ky)*np.abs(z))
        exp_xy = lambda ky: np.exp(1j*ky*y + 1j*kx*x)
        rp = lambda ky: np.exp(-2*k_parallel(ky)*zp)*Rp*kp/(k_parallel(ky) - kp)
    ############ I1 ###########################################################
        int1_re = lambda ky: np.real(-1j*(px*kx + ky*py)*exp_z(ky)*exp_xy(ky)/k_parallel(ky))
        int1_im = lambda ky: np.imag(-1j*(px*kx + ky*py)*exp_z(ky)*exp_xy(ky)/k_parallel(ky))
    
    ################### I2 ###############################################
        int2_re = lambda ky: np.real(pz*exp_z(ky)*ky*exp_xy(ky))
        int2_im = lambda ky: np.imag(pz*exp_z(ky)*ky*exp_xy(ky))
    
#        int3_re,err = integrate.quad(Int3_B_function_re, cota_inf, cota_sup)
#        int3_im,err = integrate.quad(Int3_B_function_im, cota_inf, cota_sup)
    ################### I4 ###############################################
        int3_re = lambda ky: np.real(px*1j*exp_z(ky)*kx*rp(ky)*exp_xy(ky)/k_parallel(ky))
        int3_im = lambda ky: np.imag(px*1j*exp_z(ky)*kx*rp(ky)*exp_xy(ky)/k_parallel(ky))
    
#        int4_re,err = integrate.quad(Int4_B_function_re, cota_inf, cota_sup)
#        int4_im,err = integrate.quad(Int4_B_function_im, cota_inf, cota_sup)
    ################### I5 ###############################################
        int4_re = lambda ky: np.real(py*1j*exp_z(ky)*ky*rp(ky)*exp_xy(ky)/k_parallel(ky))
        int4_im = lambda ky: np.imag(py*1j*exp_z(ky)*ky*rp(ky)*exp_xy(ky)/k_parallel(ky))
#    
#        int5_re,err = integrate.quad(Int5_B_function_re, cota_inf, cota_sup)
#        int5_im,err = integrate.quad(Int5_B_function_im, cota_inf, cota_sup)
    ################### I5 ###############################################
        int5_re = lambda ky: np.real(-pz*exp_z(ky)*ky*rp(ky)*exp_xy(ky))
        int5_im = lambda ky: np.imag(-pz*exp_z(ky)*ky*rp(ky)*exp_xy(ky))
    
#        int6_re,err = integrate.quad(Int6_B_function_re, cota_inf, cota_sup)
#        int6_im,err = integrate.quad(Int6_B_function_im, cota_inf, cota_sup)
    ################### I6 ###############################################    
        
        int_re_tot = lambda ky: int1_re(ky) + int2_re(ky) + int3_re(ky) + int4_re(ky) + int5_re(ky) 
        int_im_tot = lambda ky: int1_im(ky) + int2_im(ky) + int3_im(ky) + int4_im(ky) + int5_im(ky) 
        
        int_tot0 = lambda ky: int_tot0(ky) + a*(int_re_tot(ky) + 1j*int_im_tot(ky))
    
    
    int_final, err = integrate.quad(int_tot0, cota_inf, cota_sup)
    
    return int_final

#%%
