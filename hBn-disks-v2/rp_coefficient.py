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
#print('Importar modulos necesarios para este codigo')

#try:
#    sys.path.insert(1, path_basic)
#    from hBn_PP import hBn_lambda_p, hBn_Rp, epsilon_x, epsilon_z
#except ModuleNotFoundError:
#    print('graphene_sigma.py no se encuentra en ' + path_basic)


try:
    sys.path.insert(1, path_basic)
    from constants import constantes
except ModuleNotFoundError:
    print('constants.py no se encuentra en ' + path_basic)

pi,hb,c,alfac,mu1,mu2 = constantes()
aux = c*hb


def epsilon_x(hbw):
    
    epsi_inf = 4.87
    hbgamma = 0.87*1e-3
    f_x = 1.83
    
    num = (170.1*1e-3)**2
    
    den = hbw*(hbw + 1j*hbgamma ) - num
    
    return epsi_inf - f_x*(num/den)


def epsilon_z(hbw):
    
    epsi_inf = 2.95
    hbgamma = 0.25*1e-3
    f_x = 0.61
    
    num = (92.5*1e-3)**2
    
    den = hbw*(hbw + 1j*hbgamma ) - num
    
    return epsi_inf - f_x*(num/den)

#def hBn_lambda_p(hbw,epsi_silica):
#    """
#    Parameters
#    ----------
#    hbw : energia = hbar*omega en eV
#    carrier_density : carrier density
#    d_nano : thickness in nm
#    decay_rate : 
#    masa_eff : 
#
#    Returns
#    -------
#    devuelve la conductividad de un metal
#    """
#    epsi_x = epsilon_x(hbw)
#    epsi_z = epsilon_z(hbw)
#
#    epsi_m = np.sqrt(epsi_x*epsi_z)
#    
##    if np.imag(epsi_m) > 0:
##        epsi_m = epsi_m
##    else:
##        epsi_m = -epsi_m
#
#    epsi1 = epsi_silica(hbw)
#    epsi3 = epsi_silica(hbw)
#    
#    num = (epsi1 - epsi_m)*(epsi3 - epsi_m)
#    den = (epsi1 + epsi_m)*(epsi3 + epsi_m)
#    
#    tot = np.log(num/den)
#
#    cte_epsi = np.sqrt(epsi_x/epsi_z)
#
#########################################################################
#### el self image green tensor no funciona bien con este cambio de signo 
#### pero sin este cambio, el rp de la pole aprox en el primer polo es negativo 
#    
##    if np.imag(cte_epsi) > 0: 
##        cte_epsi = cte_epsi
##    else:
##        cte_epsi = -cte_epsi
#########################################################################    
#
#
#    return 4*np.pi*cte_epsi/tot
#
#def hBn_Rp(hbw,epsi_silica):
#    """
#    Parameters
#    ----------
#    hbw : energia = hbar*omega en eV
#    carrier_density : carrier density
#    d_nano : thickness in nm
#    decay_rate : 
#    masa_eff : 
#
#    Returns
#    -------
#    devuelve la conductividad de un metal
#    """
#    epsi_x = epsilon_x(hbw)
#    epsi_z = epsilon_z(hbw)
#    
#    lambda_p_v = hBn_lambda_p(hbw,epsi_silica) ## lambda/d
#
#    epsi1 = epsi_silica(hbw)
#    
##    c = 3*10**(14)                ### light velocity in micron/seg
##    alfac = 1/137.0359            ### fine structure
#    
#    num = -epsi_z*epsi1
#    den = epsi_x*epsi_z - epsi1**2
#    
#    tot = num/den
#
#    return tot*lambda_p_v/np.pi

#%%

#def rp_pole_aprox_v2(omegac,epsi_silica,d_nano,k_parallel_nano):     
#    """    
#    Parameters
#    ----------
#    omegac : omega/c = k0 en 1/micrometros    
#    epsi1 : epsilon del medio de arriba del plano
#    epsi2 : epsilon del medio de abajo del plano
#    hbmu : chemical potential in eV  
#    hbgama : collision frequency in eV
#    z : coordenada z
#    xD : coordenada x del dipolo 
#    yD : coordenada y del dipolo
#    zD : coordenada z del dipolo 
#    zp : posicion del plano (>0)
#    px : coordenada x del dipolo 
#    py : coordenada y del dipolo
#    pz : coordenada z del dipolo
#    Returns
#    -------
#    px,py,pz en unidades de k*alfa_eff
#    """
#
#    E = omegac*aux
#    
##    d_micros = d_nano*1e-3
#    Rp = hBn_Rp(E,epsi_silica)
#    lambda_p_v = hBn_lambda_p(E,epsi_silica)*d_nano
#    kp = 2*np.pi/lambda_p_v
#   
#    rp = Rp*kp/(k_parallel_nano - kp)
#      
#    return rp

def rp_pole_aprox(omegac,epsi_silica,d_nano,zp_nano,k_parallel_nano):     
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
    px,py,pz en unidades de k*alfa_eff
    """

    E = omegac*aux
    
#    d_micros = d_nano*1e-3
    epsi_x = epsilon_x(E) ## paralelo 
#    epsi_z = epsilon_z(E)
    
    epsi_HBN_par = epsi_x
    
    r = (1 - epsi_silica(E))/(1 + epsi_silica(E))

    sigma_2D = d_nano*(epsi_HBN_par - 1 )/(4*pi)  ## se cancela el i*omega del sigma 
 
    
    exp_fresnel =  np.exp(-2*k_parallel_nano*zp_nano)
    exp_fresnel = 1 - 2*k_parallel_nano*zp_nano + (-2*k_parallel_nano*zp_nano)**2/2 + (-2*k_parallel_nano*zp_nano)**3/6
    exp_fresnel = 0
    
    rp = k_parallel_nano/(k_parallel_nano*(1-r*exp_fresnel) - epsi_silica(E)/(2*np.pi*sigma_2D) )
    rp = k_parallel_nano/(k_parallel_nano- epsi_silica(E)/(2*np.pi*sigma_2D) )
#
#    if np.imag(rp) < 0: 
#        rp_im = -np.imag(rp)
#    else:
#        rp_im = np.imag(rp)
#    
#    
#    rp_re = np.real(rp)
    
    return rp


def rp_pole_aprox_v2(omegac,epsi_silica,d_nano,zp_nano,k_parallel_nano):     ## aproximar el numerador k_parallel a kp 
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
    px,py,pz en unidades de k*alfa_eff
    """

    E = omegac*aux
    
#    d_micros = d_nano*1e-3
    epsi_x = epsilon_x(E) ## paralelo 
#    epsi_z = epsilon_z(E)
    
    epsi_HBN_par = epsi_x
    
    r = (1 - epsi_silica(E))/(1 + epsi_silica(E))

    sigma_2D = d_nano*(epsi_HBN_par - 1 )/(4*pi)  ## se cancela el i*omega del sigma 
 
    
    exp_fresnel =  np.exp(-2*k_parallel_nano*zp_nano)
    exp_fresnel = 1 - 2*k_parallel_nano*zp_nano + (-2*k_parallel_nano*zp_nano)**2/2 + (-2*k_parallel_nano*zp_nano)**3/6
    exp_fresnel = 0
    
    kp = epsi_silica(E)/(2*np.pi*sigma_2D)
    
    rp = k_parallel_nano/(k_parallel_nano*(1-r*exp_fresnel) - epsi_silica(E)/(2*np.pi*sigma_2D) )
    rp = kp/(k_parallel_nano -  kp)
#
#    if np.imag(rp) < 0: 
#        rp_im = -np.imag(rp)
#    else:
#        rp_im = np.imag(rp)
#    
#    
#    rp_re = np.real(rp)
    
    return rp

#%%

def rp_fresnel_num(omegac,epsi_silica,d_nano,zp_nano,k_parallel_nano):     
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
#    omegac_nano = omegac*1e-3
    
    epsi_x = epsilon_x(E) ## paralelo 
#    epsi_z = epsilon_z(E)
    
    epsi_HBN_par = epsi_x
#    epsi_HBN_perp = epsi_z

#    d_micros = d_nano*1e-3



    sigma_2D = d_nano*(epsi_HBN_par - 1 )/(4*pi)  ## se cancela el i*omega del sigma 


    r_prima = (1 - epsi_silica(E)/(2*pi*k_parallel_nano*sigma_2D) )**(-1)

    r = (1 - epsi_silica(E))/(1 + epsi_silica(E))


    exp_fresnel =  np.exp(-2*k_parallel_nano*zp_nano)

    r_fresnel = r_prima/(1 - r*r_prima*exp_fresnel )


    return r_fresnel

#%%
