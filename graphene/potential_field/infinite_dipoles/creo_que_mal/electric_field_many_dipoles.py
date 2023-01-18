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
    
    cond = 4*pi*alfac*sigma_DL(E,hbmu,hbgama)
    Rp = 2*epsi1/(epsi1 + epsi2)
    alfa_p = 1j*(epsi1 + epsi2)/(cond*cte1)
    kp = alfa_p*k1    
    
    v = c/int_v
    cota_inf = 1
    cota_sup = 1000    

    int_tot0_x = 0 
    int_tot0_y = 0
    int_tot0_z = 0
    for n in range(N):
#        print(n)
        kx = omegac*c/v - 2*pi*n/a
        
        k_parallel = lambda ky: np.sqrt(kx**2 + ky**2)
        exp_z = lambda ky: np.e**(-k_parallel(ky)*np.abs(z))
        exp_xy = lambda ky: np.e**(1j*ky*y + 1j*kx*x)
        rp = lambda ky: Rp*kp/(k_parallel(ky) - kp)
    ################### I1 ################################################
        Int1_B_function_re_x = lambda ky: np.real(1j*rp(ky)*(px*kx + ky*py)*exp_z(ky)*exp_xy(ky)/k_parallel(ky))
        Int1_B_function_im_x = lambda ky: np.imag(1j*rp(ky)*(px*kx + ky*py)*exp_z(ky)*exp_xy(ky)/k_parallel(ky))

        int1_re_x,err = integrate.quad(Int1_B_function_re_x, cota_inf, cota_sup)
        int1_im_x,err = integrate.quad(Int1_B_function_im_x, cota_inf, cota_sup)

        Int1_B_function_re_y = lambda ky: np.real(1j*ky*rp(ky)*(px*kx + ky*py)*exp_z(ky)*exp_xy(ky)/k_parallel(ky))
        Int1_B_function_im_y = lambda ky: np.imag(1j*ky*rp(ky)*(px*kx + ky*py)*exp_z(ky)*exp_xy(ky)/k_parallel(ky))

        int1_re_y,err = integrate.quad(Int1_B_function_re_y, cota_inf, cota_sup)
        int1_im_y,err = integrate.quad(Int1_B_function_im_y, cota_inf, cota_sup)
 
        Int1_B_function_re_z = lambda ky: np.real(1j*rp(ky)*(px*kx + ky*py)*exp_z(ky)*exp_xy(ky))
        Int1_B_function_im_z = lambda ky: np.imag(1j*rp(ky)*(px*kx + ky*py)*exp_z(ky)*exp_xy(ky))

        int1_re_z,err = integrate.quad(Int1_B_function_re_z, cota_inf, cota_sup)
        int1_im_z,err = integrate.quad(Int1_B_function_im_z, cota_inf, cota_sup)    

    ################### I2 ###############################################
        Int2_B_function_re_x = lambda ky: np.real(-rp(ky)*pz*np.sign(z)*exp_z(ky)*exp_xy(ky))
        Int2_B_function_im_x = lambda ky: np.imag(-rp(ky)*pz*np.sign(z)*exp_z(ky)*exp_xy(ky))
    
        int2_re_x,err = integrate.quad(Int2_B_function_re_x, cota_inf, cota_sup)
        int2_im_x,err = integrate.quad(Int2_B_function_im_x, cota_inf, cota_sup)
        
        Int2_B_function_re_y = lambda ky: np.real(-ky*rp(ky)*pz*np.sign(z)*exp_z(ky)*ky*exp_xy(ky))
        Int2_B_function_im_y = lambda ky: np.imag(-ky*rp(ky)*pz*np.sign(z)*exp_z(ky)*ky*exp_xy(ky))
    
        int2_re_y,err = integrate.quad(Int2_B_function_re_y, cota_inf, cota_sup)
        int2_im_y,err = integrate.quad(Int2_B_function_im_y, cota_inf, cota_sup)


        Int2_B_function_re_z = lambda ky: np.real(-rp(ky)*pz*np.sign(z)*exp_z(ky)*exp_xy(ky)*k_parallel(ky))
        Int2_B_function_im_z = lambda ky: np.imag(-rp(ky)*pz*np.sign(z)*exp_z(ky)*exp_xy(ky)*k_parallel(ky))
    
        int2_re_z,err = integrate.quad(Int2_B_function_re_z, cota_inf, cota_sup)
        int2_im_z,err = integrate.quad(Int2_B_function_im_z, cota_inf, cota_sup)        
    ################### I3 ###############################################
        Int3_B_function_re_x = lambda ky: np.real(-(px*kx + py*ky)*1j*exp_z(ky)*exp_xy(ky)/k_parallel(ky))
        Int3_B_function_im_x = lambda ky: np.imag(-(px*kx + py*ky)*1j*exp_z(ky)*exp_xy(ky)/k_parallel(ky))
    
        int3_re_x,err = integrate.quad(Int3_B_function_re_x, cota_inf, cota_sup)
        int3_im_x,err = integrate.quad(Int3_B_function_im_x, cota_inf, cota_sup)
        
        
        Int3_B_function_re_y = lambda ky: np.real(-ky*(px*kx + py*ky)*1j*exp_z(ky)*exp_xy(ky)/k_parallel(ky))
        Int3_B_function_im_y = lambda ky: np.imag(-ky*(px*kx + py*ky)*1j*exp_z(ky)*exp_xy(ky)/k_parallel(ky))
    
        int3_re_y,err = integrate.quad(Int3_B_function_re_y, cota_inf, cota_sup)
        int3_im_y,err = integrate.quad(Int3_B_function_im_y, cota_inf, cota_sup)


        Int3_B_function_re_z = lambda ky: np.real(-(px*kx + py*ky)*1j*exp_z(ky)*exp_xy(ky))
        Int3_B_function_im_z = lambda ky: np.imag(-(px*kx + py*ky)*1j*exp_z(ky)*exp_xy(ky))
    
        int3_re_z,err = integrate.quad(Int3_B_function_re_z, cota_inf, cota_sup)
        int3_im_z,err = integrate.quad(Int3_B_function_im_z, cota_inf, cota_sup)
        
        
    ################### I4 ###############################################
        Int4_B_function_re_x = lambda ky: np.real(pz*np.sign(z)*exp_z(ky)*rp(ky)*exp_xy(ky))
        Int4_B_function_im_x = lambda ky: np.imag(pz*np.sign(z)*exp_z(ky)*rp(ky)*exp_xy(ky))
#    
        int4_re_x,err = integrate.quad(Int4_B_function_re_x, cota_inf, cota_sup)
        int4_im_x,err = integrate.quad(Int4_B_function_im_x, cota_inf, cota_sup)


        Int4_B_function_re_y = lambda ky: np.real(ky*pz*np.sign(z)*exp_z(ky)*ky*rp(ky)*exp_xy(ky))
        Int4_B_function_im_y = lambda ky: np.imag(ky*pz*np.sign(z)*exp_z(ky)*ky*rp(ky)*exp_xy(ky))
#    
        int4_re_y,err = integrate.quad(Int4_B_function_re_y, cota_inf, cota_sup)
        int4_im_y,err = integrate.quad(Int4_B_function_im_y, cota_inf, cota_sup)


        Int4_B_function_re_z = lambda ky: np.real(pz*np.sign(z)*exp_z(ky)*rp(ky)*exp_xy(ky)*k_parallel(ky))
        Int4_B_function_im_z = lambda ky: np.imag(pz*np.sign(z)*exp_z(ky)*rp(ky)*exp_xy(ky)*k_parallel(ky))
#    
        int4_re_z,err = integrate.quad(Int4_B_function_re_z, cota_inf, cota_sup)
        int4_im_z,err = integrate.quad(Int4_B_function_im_z, cota_inf, cota_sup)


    ################### I5 ##############################################
        
        int_re_tot_x = int1_re_x + int2_re_x + int3_re_x + int4_re_x
        int_im_tot_x = int1_im_x + int2_im_x + int3_im_x + int4_im_x  
        
        int_tot0_x = int_tot0_x - a*1j*kx*(int_re_tot_x + 1j*int_im_tot_x)


        int_re_tot_y = int1_re_y + int2_re_y + int3_re_y + int4_re_y
        int_im_tot_y = int1_im_y + int2_im_y + int3_im_y + int4_im_y 

        int_tot0_y = int_tot0_y - a*1j*(int_re_tot_y + 1j*int_im_tot_y)
        
        
        int_re_tot_z = int1_re_z + int2_re_z + int3_re_z + int4_re_z
        int_im_tot_z = int1_im_z + int2_im_z + int3_im_z + int4_im_z  
        
        int_tot0_z = int_tot0_z + np.sign(z)*a*(int_re_tot_z + 1j*int_im_tot_z)
    
    return int_tot0_x, int_tot0_y, int_tot0_z

#%%
    


