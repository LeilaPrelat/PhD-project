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
path_constants =  path_basic.replace('/potential_field/many_potential_Javier_formula/creo_que_mal','')
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

    int_tot0_x = 0 
    int_tot0_y = 0
    int_tot0_z = 0
    list_dipolos = np.linspace(-N,N,2*N+1)
    for n in list_dipolos:
#        print(n)
        kx = omegac*c/v - 2*pi*n/a
        
        k_parallel = lambda ky: np.sqrt(kx**2 + ky**2)
        exp_z1 = lambda ky: np.exp(-k_parallel(ky)*np.abs(z))
        exp_z2 = lambda ky: np.exp(-k_parallel(ky)*(2*np.abs(zp) + np.abs(z)))
        
        
        exp_xy = lambda ky: np.e**(1j*ky*y + 1j*kx*x)
        rp = lambda ky: Rp*kp/(k_parallel(ky) - kp)
        
        #multiplicar todo por -1j*k componente
    ################### I1 ################################################
        Int1_B_function_re_x = lambda ky: np.real(-kx*(px*kx + ky*py)*exp_z1(ky)*exp_xy(ky)/k_parallel(ky))
        Int1_B_function_im_x = lambda ky: np.imag(-kx*(px*kx + ky*py)*exp_z1(ky)*exp_xy(ky)/k_parallel(ky))

        int1_re_x,err = integrate.quad(Int1_B_function_re_x, cota_inf, cota_sup)
        int1_im_x,err = integrate.quad(Int1_B_function_im_x, cota_inf, cota_sup)
    
    ################### I2 ###############################################
        Int2_B_function_re_x = lambda ky: np.real(-1j*kx*pz*np.sign(z)*exp_z1(ky)*exp_xy(ky))
        Int2_B_function_im_x = lambda ky: np.imag(-1j*kx*pz*np.sign(z)*exp_z1(ky)*exp_xy(ky))
    
        int2_re_x,err = integrate.quad(Int2_B_function_re_x, cota_inf, cota_sup)
        int2_im_x,err = integrate.quad(Int2_B_function_im_x, cota_inf, cota_sup)
    ################### I3 ###############################################
        Int3_B_function_re_x = lambda ky: np.real(px*(kx**2)*rp(ky)*exp_xy(ky)*exp_z2(ky)/k_parallel(ky))
        Int3_B_function_im_x = lambda ky: np.imag(px*(kx**2)*rp(ky)*exp_xy(ky)*exp_z2(ky)/k_parallel(ky))
    
        int3_re_x,err = integrate.quad(Int3_B_function_re_x, cota_inf, cota_sup)
        int3_im_x,err = integrate.quad(Int3_B_function_im_x, cota_inf, cota_sup)
    ################### I4 ###############################################
        Int4_B_function_re_x = lambda ky: np.real(py*ky*kx*rp(ky)*exp_xy(ky)*exp_z2(ky)/k_parallel(ky))
        Int4_B_function_im_x = lambda ky: np.imag(py*ky*kx*rp(ky)*exp_xy(ky)*exp_z2(ky)/k_parallel(ky))
#    
        int4_re_x,err = integrate.quad(Int4_B_function_re_x, cota_inf, cota_sup)
        int4_im_x,err = integrate.quad(Int4_B_function_im_x, cota_inf, cota_sup)
    ################### parte con x ###############################################
        Int5_B_function_re_x = lambda ky: np.real(1j*kx*pz*np.sign(z)*rp(ky)*exp_xy(ky)*exp_z2(ky))
        Int5_B_function_im_x = lambda ky: np.imag(1j*kx*pz*np.sign(z)*rp(ky)*exp_xy(ky)*exp_z2(ky))
    
        int5_re_x,err = integrate.quad(Int5_B_function_re_x, cota_inf, cota_sup)
        int5_im_x,err = integrate.quad(Int5_B_function_im_x, cota_inf, cota_sup)
    ################### parte con y ###############################################    
        
        Int1_B_function_re_y = lambda ky: np.real(-ky*(px*kx + ky*py)*exp_xy(ky)/k_parallel(ky))
        Int1_B_function_im_y = lambda ky: np.imag(-ky*(px*kx + ky*py)*exp_xy(ky)/k_parallel(ky))

        int1_re_y,err = integrate.quad(Int1_B_function_re_y, cota_inf, cota_sup)
        int1_im_y,err = integrate.quad(Int1_B_function_im_y, cota_inf, cota_sup)
    
    ################### I2 ###############################################
        Int2_B_function_re_y = lambda ky: np.real(-1j*pz*np.sign(z)*exp_z1(ky)*exp_xy(ky))
        Int2_B_function_im_y = lambda ky: np.imag(-1j*pz*np.sign(z)*exp_z1(ky)*exp_xy(ky))
    
        int2_re_y,err = integrate.quad(Int2_B_function_re_y, cota_inf, cota_sup)
        int2_im_y,err = integrate.quad(Int2_B_function_im_y, cota_inf, cota_sup)
    ################### I3 ###############################################
        Int3_B_function_re_y = lambda ky: np.real(px*kx*ky*rp(ky)*exp_xy(ky)*exp_z2(ky)/k_parallel(ky))
        Int3_B_function_im_y = lambda ky: np.imag(px*kx*ky*rp(ky)*exp_xy(ky)*exp_z2(ky)/k_parallel(ky))
    
        int3_re_y,err = integrate.quad(Int3_B_function_re_y, cota_inf, cota_sup)
        int3_im_y,err = integrate.quad(Int3_B_function_im_y, cota_inf, cota_sup)
    ################### I4 ###############################################
        Int4_B_function_re_y = lambda ky: np.real(py*(ky**2)*rp(ky)*exp_xy(ky)*exp_z2(ky)/k_parallel(ky))
        Int4_B_function_im_y = lambda ky: np.imag(py*(ky**2)*rp(ky)*exp_xy(ky)*exp_z2(ky)/k_parallel(ky))
#    
        int4_re_y,err = integrate.quad(Int4_B_function_re_y, cota_inf, cota_sup)
        int4_im_y,err = integrate.quad(Int4_B_function_im_y, cota_inf, cota_sup)
    ################### I5 ###############################################
        Int5_B_function_re_y = lambda ky: np.real(1j*ky*pz*np.sign(z)*rp(ky)*exp_xy(ky)*exp_z2(ky))
        Int5_B_function_im_y = lambda ky: np.imag(1j*ky*pz*np.sign(z)*rp(ky)*exp_xy(ky)*exp_z2(ky))
    
        int5_re_y,err = integrate.quad(Int5_B_function_re_y, cota_inf, cota_sup)
        int5_im_y,err = integrate.quad(Int5_B_function_im_y, cota_inf, cota_sup)
        
    ################### parte con z : multiplicar por k parallel*np.sign(z) ###############################################   
        
    
        Int1_B_function_re_z = lambda ky: np.real(-1j*(px*kx + ky*py)*exp_z1(ky)*exp_xy(ky)*np.sign(z))
        Int1_B_function_im_z = lambda ky: np.imag(-1j*(px*kx + ky*py)*exp_z1(ky)*exp_xy(ky)*np.sign(z))

        int1_re_z,err = integrate.quad(Int1_B_function_re_z, cota_inf, cota_sup)
        int1_im_z,err = integrate.quad(Int1_B_function_im_z, cota_inf, cota_sup)
    
    ################### I2 ###############################################
        Int2_B_function_re_z = lambda ky: np.real(pz*np.sign(z)*np.sign(z)*exp_z1(ky)*exp_xy(ky)*k_parallel(ky))
        Int2_B_function_im_z = lambda ky: np.imag(pz*np.sign(z)*np.sign(z)*exp_z1(ky)*exp_xy(ky)*k_parallel(ky))
    
        int2_re_z,err = integrate.quad(Int2_B_function_re_z, cota_inf, cota_sup)
        int2_im_z,err = integrate.quad(Int2_B_function_im_z, cota_inf, cota_sup)
    ################### I3 ###############################################
        Int3_B_function_re_z = lambda ky: np.real(px*kx*1j*rp(ky)*exp_xy(ky)*exp_z2(ky)*np.sign(z))
        Int3_B_function_im_z = lambda ky: np.imag(px*kx*1j*rp(ky)*exp_xy(ky)*exp_z2(ky)*np.sign(z))
    
        int3_re_z,err = integrate.quad(Int3_B_function_re_z, cota_inf, cota_sup)
        int3_im_z,err = integrate.quad(Int3_B_function_im_z, cota_inf, cota_sup)
    ################### I4 ###############################################
        Int4_B_function_re_z = lambda ky: np.real(py*ky*1j*rp(ky)*exp_xy(ky)*exp_z2(ky)*np.sign(z))
        Int4_B_function_im_z = lambda ky: np.imag(py*ky*1j*rp(ky)*exp_xy(ky)*exp_z2(ky)*np.sign(z))
#    
        int4_re_z,err = integrate.quad(Int4_B_function_re_z, cota_inf, cota_sup)
        int4_im_z,err = integrate.quad(Int4_B_function_im_z, cota_inf, cota_sup)
    ################### I5 ###############################################
        Int5_B_function_re_z = lambda ky: np.real(-pz*np.sign(z)*rp(ky)*exp_xy(ky)*exp_z2(ky)*np.sign(z)*k_parallel(ky))
        Int5_B_function_im_z = lambda ky: np.imag(-pz*np.sign(z)*rp(ky)*exp_xy(ky)*exp_z2(ky)*np.sign(z)*k_parallel(ky))
    
        int5_re_z,err = integrate.quad(Int5_B_function_re_z, cota_inf, cota_sup)
        int5_im_z,err = integrate.quad(Int5_B_function_im_z, cota_inf, cota_sup)
    
    ################### parte con z ###############################################   
    
        int_re_tot_x = int1_re_x + int2_re_x + int3_re_x + int4_re_x
        int_im_tot_x = int1_im_x + int2_im_x + int3_im_x + int4_im_x  
        
        int_tot0_x = int_tot0_x - a*(int_re_tot_x + 1j*int_im_tot_x)


        int_re_tot_y = int1_re_y + int2_re_y + int3_re_y + int4_re_y
        int_im_tot_y = int1_im_y + int2_im_y + int3_im_y + int4_im_y 

        int_tot0_y = int_tot0_y - a*(int_re_tot_y + 1j*int_im_tot_y)
        
        
        int_re_tot_z = int1_re_z + int2_re_z + int3_re_z + int4_re_z
        int_im_tot_z = int1_im_z + int2_im_z + int3_im_z + int4_im_z  
        
        int_tot0_z = int_tot0_z + np.sign(z)*a*(int_re_tot_z + 1j*int_im_tot_z)
    
#    final = np.abs()**2 + np.abs(int_tot0_y)**2 + np.abs(int_tot0_z)**2
    return int_tot0_x, int_tot0_y, int_tot0_z

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
    analitico
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

    term_final_x = 0
    term_final_y = 0
    term_final_z = 0
    list_dipolos = np.linspace(-N,N,2*N+1)
    for n in list_dipolos:
#        print(n)
        kx = omegac*c/v - 2*pi*n/a
        alfa_x = kx/k1
        
        alfa_y = np.sqrt(alfa_p**2 - alfa_x**2)
        ky = alfa_y*k1
        
        exp_z1 = np.exp(-2*alfa_p*k1*zp)
        exp_z2 = np.exp(-alfa_p*k1*np.abs(z))
        
        exp_xy = np.exp(1j*alfa_y*k1*y + 1j*alfa_x*k1*x)
        
        term1 = 1j*(px*alfa_x + py*alfa_y)/alfa_p
        term2 = pz*np.sign(z)
        term3 = 1j*Rp*(px*alfa_x + py*alfa_y)*exp_z1
        term4 = pz*alfa_p*np.sign(z)*Rp*exp_z1
        
        
        factor_final = exp_xy*exp_z2*k1*a
        
        ################### I1 ################################################
        
        term_final_x = - 1j*(term1 + term2 + term3 + term4)*factor_final*kx + term_final_x 
        term_final_y = - 1j*(term1 + term2 + term3 + term4)*factor_final*ky + term_final_y
        term_final_z = kp*np.sign(z)*(term1 + term2 + term3 + term4)*factor_final + term_final_z
    
    ################### parte con z ###############################################   
#    
#    final = np.abs()**2 + np.abs(term_final_y)**2 + np.abs(term_final_z)**2
    return term_final_x, term_final_y, term_final_z

#%%
    
