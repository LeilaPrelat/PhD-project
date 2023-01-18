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

def Ex_many_dipoles_num2(omegac,epsi1,epsi2,hbmu,hbgama,x,y,z,a,zp,int_v,px,py,pz):
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
    omega = omegac*c

    # x_y = ky/k0

    k1 = omegac

    v = c/int_v
    cota_inf = 1
    cota_sup = 1000    

    int_tot0 = 0 
#    print('hola')
    cond = 4*pi*alfac*sigma_DL(E,hbmu,hbgama)
    alfa_p = 1j*(epsi1 + epsi2)/(cond)
    kp = alfa_p*k1    
    
    v = c/int_v
    int_tot0 = 0 
#    print('hola')

#    print('hola')
    
    nmax =  (np.real(kp) - omega/v)*a/(2*np.pi)    #maximum order 
    nmax = omegac*(np.real(alfa_p) - int_v)*a/(2*np.pi)
    nmax = int(nmax)
#    print('k: f1',omegac)
    
    print('alfa_p: f1', kp)
#    print('a/2pi: f1', a/(2*np.pi))
#    print(nmax)
    
    N = 20
    list_dipolos = np.linspace(0,2*N,2*N + 1)
    
    
    factor = a*(2*np.pi)**2
    for n in list_dipolos:
#        print(n)
        kx = omegac*c/v + 2*pi*n/a
        
        k_parallel = lambda ky: np.sqrt(kx**2 + ky**2)
        # exp_z1 = lambda ky: np.exp(-k_parallel(ky)*np.abs(z))
        exp_z = lambda ky: np.exp(-k_parallel(ky)*2*np.abs(zp))
        
        exp_xy1 = lambda ky: np.exp(1j*kx*x)*np.cos(ky*y)
        exp_xy2 = lambda ky: np.exp(1j*kx*x)*np.sin(ky*y)        
        # rp = lambda ky: Rp*kp/(k_parallel(ky) - kp)

        rp_num = lambda ky: epsi2*1j*k_parallel(ky) - epsi1*1j*k_parallel(ky) - cond*k_parallel(ky)**2/omegac
        rp_den = lambda ky: epsi2*1j*k_parallel(ky) + epsi1*1j*k_parallel(ky) - cond*k_parallel(ky)**2/omegac
        rp = lambda ky: rp_num(ky)/rp_den(ky)
        
    ################### I1 ################################################
        arg = kx*y
        K1 = special.kn(1,arg)
        K0 = special.kn(0,arg)
        term1 = kx*(-px*1j*K0 + 2*py*K1) + 1j*pz*np.sign(z)/y
        
    ################### I1 ################################################
        Int1_B_function_re = lambda ky: np.real(2*1j*(px*kx*exp_xy1(ky) + ky*py*exp_xy2(ky))*exp_z(ky)*rp(ky)/k_parallel(ky))
        Int1_B_function_im = lambda ky: np.imag(2*1j*(px*kx*exp_xy1(ky) + ky*py*exp_xy2(ky))*exp_z(ky)*rp(ky)/k_parallel(ky))

        int1_re,err = integrate.quad(Int1_B_function_re, cota_inf, cota_sup)
        int1_im,err = integrate.quad(Int1_B_function_im, cota_inf, cota_sup)
    
    ################### I2 ###############################################
        Int2_B_function_re = lambda ky: np.real(-2*pz*exp_z(ky)*exp_xy1(ky)*rp(ky)/k_parallel(ky))
        Int2_B_function_im = lambda ky: np.imag(-2*pz*exp_z(ky)*exp_xy1(ky)*rp(ky)/k_parallel(ky))

        int2_re,err = integrate.quad(Int2_B_function_re, cota_inf, cota_sup)
        int2_im,err = integrate.quad(Int2_B_function_im, cota_inf, cota_sup)
    
    ################### I3 ###############################################

    ################### I6 ###############################################    
        
        int_tot0 = int_tot0 + (int1_re + 1j*int1_im + int2_re + 1j*int2_im + term1)*1j*kx
#        print(int_tot0)    

        print(int_tot0)
    return int_tot0/factor


#%%

def Ex_many_dipoles_pole_approx(omegac,epsi1,epsi2,hbmu,hbgama,x,y,z,a,zp,int_v,px,py,pz):
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
    omega = omegac*c

    # x_y = ky/k0

    k1 = omegac

    cota_inf = 1
    cota_sup = 1000    

#    print('hola')
    cond = 4*pi*alfac*sigma_DL(E,hbmu,hbgama)
    Rp = 2*epsi1/(epsi1 + epsi2)
    alfa_p = 1j*(epsi1 + epsi2)/(cond)
    kp = alfa_p*k1    
    
    v = c/int_v
    int_tot0 = 0 
#    print('hola')

#    print('hola')
    
    nmax =  (np.real(kp) - omega/v)*a/(2*np.pi)    #maximum order 
    nmax = omegac*(np.real(alfa_p) - int_v)*a/(2*np.pi)
    nmax = int(nmax)
#    print('k: f2',omegac)
#    
#    print('alfa_p: f2', alfa_p - int_v)
#    print('a/2pi: f2', a/(2*np.pi))
#    print(nmax)
#    print(nmax)

    if nmax > 0:
        list_dipolos = np.linspace(0,nmax,nmax+1)
    else: 
        list_dipolos = np.linspace(nmax,0,-nmax+1)

    N = 20
    list_dipolos = np.linspace(0,2*N,2*N + 1)
    
    factor = a*(2*np.pi)**2
    for n in list_dipolos:
#        print(n)
        kx = omegac*c/v + 2*pi*n/a
        
        k_parallel = lambda ky: np.sqrt(kx**2 + ky**2)
        # exp_z1 = lambda ky: np.exp(-k_parallel(ky)*np.abs(z))
        exp_z = lambda ky: np.exp(-k_parallel(ky)*2*np.abs(zp))
        
        exp_xy1 = lambda ky: np.exp(1j*kx*x)*np.cos(ky*y)
        exp_xy2 = lambda ky: np.exp(1j*kx*x)*np.sin(ky*y)  
        # rp = lambda ky: Rp*kp/(k_parallel(ky) - kp)


        rp  = lambda ky : Rp*kp/(k_parallel(ky) - kp)

        
    ################### I1 ################################################
        arg = kx*y
        K1 = special.kn(1,arg)
        K0 = special.kn(0,arg)
        term1 = kx*(-px*1j*K0 + 2*py*K1) + 1j*pz*np.sign(z)/y
        
    ################### I1 ################################################
        Int1_B_function_re = lambda ky: np.real(2*1j*(px*kx*exp_xy1(ky) + ky*py*exp_xy2(ky))*exp_z(ky)*rp(ky)/k_parallel(ky))
        Int1_B_function_im = lambda ky: np.imag(2*1j*(px*kx*exp_xy1(ky) + ky*py*exp_xy2(ky))*exp_z(ky)*rp(ky)/k_parallel(ky))

        int1_re,err = integrate.quad(Int1_B_function_re, cota_inf, cota_sup)
        int1_im,err = integrate.quad(Int1_B_function_im, cota_inf, cota_sup)
    
    ################### I2 ###############################################
        Int2_B_function_re = lambda ky: np.real(-2*pz*exp_z(ky)*exp_xy1(ky)*rp(ky)/k_parallel(ky))
        Int2_B_function_im = lambda ky: np.imag(-2*pz*exp_z(ky)*exp_xy1(ky)*rp(ky)/k_parallel(ky))

        int2_re,err = integrate.quad(Int2_B_function_re, cota_inf, cota_sup)
        int2_im,err = integrate.quad(Int2_B_function_im, cota_inf, cota_sup)
    
    ################### I3 ###############################################

    ################### I6 ###############################################    
        
        int_tot0 = int_tot0 + (int1_re + 1j*int1_im + int2_re + 1j*int2_im + term1)*1j*kx
#        print(int_tot0)    
#
#        print(int_tot0)
    return int_tot0/factor

#%%

def Ex_many_dipoles_ana2(omegac,epsi1,epsi2,hbmu,hbgama,x,y,z,a,zp,int_v,px,py,pz):
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
    omega = omegac*c
    # x_y = ky/k0

    k1 = omegac
    
    cond = 4*pi*alfac*sigma_DL(E,hbmu,hbgama)
    Rp = 2*epsi1/(epsi1 + epsi2)
    alfa_p = 1j*(epsi1 + epsi2)/(cond)
    kp = alfa_p*k1    
    
    v = c/int_v
    int_tot0 = 0 
#    print('hola')
    
    nmax =  (np.real(kp) - omega/v)*a/(2*np.pi)    #maximum order 
    nmax = omegac*(np.real(alfa_p) - int_v)*a/(2*np.pi)
    nmax = int(nmax)
#    print('k: f3',omegac)
#    
#    print('alfa_p: f3', alfa_p - int_v)
#    print('a/2pi: f3', a/(2*np.pi))
#    print(nmax)
    
    if nmax > 0:
        list_dipolos = np.linspace(0,nmax,nmax+1)
    else: 
        list_dipolos = np.linspace(nmax,0,-nmax+1)
    
    factor = a*(2*np.pi)**2
    for n in list_dipolos:
#        print(n)
        kx = omegac*int_v + 2*pi*n/a
        ky = np.sqrt(kp**2 - kx**2)  # tiene que haber un menos si o si
        
#        if np.imag(ky) > 0 : #convergencia
#            ky = ky
#        else:
#            ky = - ky
        
        # k_parallel = lambda ky: np.sqrt(kx**2 + ky**2)
        # exp_z1 = lambda ky: np.exp(-k_parallel(ky)*np.abs(z))
        exp_z = np.exp(-kp*2*np.abs(zp))
        
        exp_xy =  np.exp(1j*kx*x)*np.cos(ky*y)
        # rp = lambda ky: Rp*kp/(k_parallel(ky) - kp)
        
    ################### I1 ################################################
        arg = kx*y
        K1 = special.kn(1,arg)
        K0 = special.kn(0,arg)
        
        # print(arg)
        
        term1 = kx*(-px*1j*K0 + 2*py*K1) + 1j*pz*np.sign(z)/y
        
        term2 = -2*np.pi*Rp*kp*px*kx*exp_xy*exp_z/ky
        
        term3 = -pz*2*np.pi*1j*Rp*(kp**2)*exp_xy*exp_z/ky
        
    ################### I1 ################################################
    
        int_tot0 = (term1 + term2 + term3)*1j*kx + int_tot0
        
#        print(int_tot0)
        
    return int_tot0/factor

#%%

