#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 19 20:56:58 2020

@author: leila

Ex, Ey, Ez del campo inducido por un dipolo
el dipolo esta en x,y,z = 0,0,0

"""
import numpy as np
import os 
import sys
from scipy import integrate

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')

try:
    sys.path.insert(1, path_basic)
    from coef_Fresnel_ICFO import coefTM_TE
except ModuleNotFoundError:
    print('coef_Fresnel_ICFO.py no se encuentra en ' + path_basic)

try:
    sys.path.insert(1, path_basic)
    from polarizability import alpha_function
except ModuleNotFoundError:
    print('polarizability.py no se encuentra en ' + path_basic)

try:
    sys.path.insert(1, path_basic)
    from green_tensor_ICFO import green_tensor
except ModuleNotFoundError:
    print('green_tensor_ICFO.py no se encuentra en ' + path_basic)

try:
    sys.path.insert(1, path_basic)
    from constants_plane import constantes
except ModuleNotFoundError:
    print('constants_plane.py no se encuentra en ' + path_basic)

pi,hb,c,alfac,hbargama,mu1,mu2,epsi1,epsi2 = constantes()
aux = c*hb

#%%
    
#print('Definir la funcion Hz para el medio 1 y el 2')

def fieldE_ind_dip(omegac,hbmu,theta_inc,z_plane,z_electron,x,y,epsi_h,px,py,pz): 
    """
    Parameters
    ----------
    omegac : omega/c = k0 en 1/micrometros    
    hbmu: potencial quimico del grafeno en eV
    theta_inc: theta incidente en radianes : kz incidente y kx incidente
    z_plane: posicion del plano en micrometros (z_plane<0)
    z_electron: posicion del electron (z_electron>0) en micrometros
    x: x en micrometros
    y: y en micrometros
    epsi_h = epsilon del host medium (material del plano)
    px : coordenada x del dipolo 
    py : coordenada y del dipolo 
    pz : coordenada z del dipolo 
    Returns
        Ex, Ey, Ez del campo inducido por un dipolo
        el dipolo esta en x,y,z = 0,0,0
    -------
    """
    k0 = omegac #=omega/c
    xz1 = np.sqrt(epsi1*mu1)*np.cos(theta_inc)
    kz1 = xz1*k0
    # x_x = kx/k0
    
    alfa = alpha_function(omegac,hbmu,theta_inc,z_plane)
    green = green_tensor(omegac,hbmu,theta_inc,z_plane,epsi_h,px,py,pz)
    factor = 1/alfa - green 
    k1_2 = epsi1*mu1*(k0**2)
    k1 = np.sqrt(k1_2)
    factor_2 = factor*(k1_2**2)
    factor_final = factor_2*(np.e**(1j*kz1*z_electron)/(kz1**2))
    
    def coeficientes():
        
        coeff = coefTM_TE(omegac,hbmu,theta_inc,z_plane)
        coef_rs, coef_ts, coef_rp, coef_tp = coeff  
        

        coef_Bref = complex(coef_rs)
        coef_Bt = complex(coef_ts)
        coef_Aref = complex(coef_rp)
        coef_At = complex(coef_tp)
        
        return coef_Bref,coef_Bt,coef_Aref,coef_At

    def epsi(medio):
        if medio ==1:
            epsirta = epsi1
        elif medio ==2:
            epsirta = epsi2
        return epsirta 
        
    def mu(medio):
        if medio ==1:
            murta = mu1
        elif medio ==2:
            murta = mu2
        return murta 
            
    # def xt(medio):  # kz/k0
        
    #     def xt2(medio): # kz^2/k^0
    #         xt2rta = -x_x**2 + mu(medio)*epsi(medio)
    #         return xt2rta
        
    #     inside = xt2(medio)+0j
    #     xtmedio = (inside)**(1/2)   

    #     if (xtmedio*z_plane_barra).imag >= 0:
    # 	    xtmedio = xtmedio
    #     else:
    # 	    xtmedio = -xtmedio
    #     return xtmedio
    
    f1 = lambda qx, qy : qx*(np.e**(1j*(qx*x + qy*y)))/(np.sqrt(qx**2 + qy**2))    
    f2 = lambda qx, qy : np.sqrt(qx**2 + qy**2)*np.e**(1j*(qx*x + qy*y))
    Int1 = integrate.dblquad(f1, -np.inf, np.inf, lambda qx: -np.inf, lambda x: np.inf)
    Int2 = integrate.dblquad(f2, -np.inf, np.inf, lambda qx: -np.inf, lambda x: np.inf)
    
    
    coef_rs, coef_ts, coef_rp, coef_tp = coeficientes()


    term1 = Int1*(coef_rs + 1)
    term2 = Int1*(coef_rp + 1)*(kz1/k1)
    Ex = -term1 + term2
    
    Ey = term1 + term2
    
    Ez = Int2*(coef_rp + 1)/k1
    

    return np.array([Ex,Ey,Ez])*factor_final        
    
#%%        

def fieldE_ind_dip_x(omegac,hbmu,theta_inc,z_plane,z_electron,x,y,epsi_h,px,py,pz): 
    """
    Parameters
    ----------
    omegac : omega/c = k0 en 1/micrometros    
    hbmu: potencial quimico del grafeno en eV
    theta_inc: theta incidente en radianes : kz incidente y kx incidente
    z_plane: posicion del plano en micrometros (z_plane<0)
    z_electron: posicion del electron (z_electron>0) en micrometros
    x: x en micrometros
    y: y en micrometros
    epsi_h = epsilon del host medium (material del plano)
    px : coordenada x del dipolo 
    py : coordenada y del dipolo 
    pz : coordenada z del dipolo 
    Returns
        Ex del campo inducido por un dipolo
        el dipolo esta en x,y,z = 0,0,0
    -------
    """
    k0 = omegac #=omega/c
    xz1 = np.sqrt(epsi1*mu1)*np.cos(theta_inc)
    kz1 = xz1*k0
    # x_x = kx/k0
    
    alfa = alpha_function(omegac,hbmu,theta_inc,z_plane)
    
    # problema de convergencia con el tensor de green 
#    green = green_tensor(omegac,hbmu,theta_inc,z_plane,epsi_h,px,py,pz) 
#    factor = 1/alfa - green 
    factor = 1/alfa 
    k1_2 = epsi1*mu1*(k0**2)
    k1 = np.sqrt(k1_2)
    factor_2 = factor*(k1_2**2)
    factor_final = factor_2*(np.e**(1j*kz1*z_electron)/(kz1**2))
    
    def coeficientes():
        
        coeff = coefTM_TE(omegac,hbmu,theta_inc,z_plane)
        coef_rs, coef_ts, coef_rp, coef_tp = coeff  
        

        coef_Bref = complex(coef_rs)
        coef_Bt = complex(coef_ts)
        coef_Aref = complex(coef_rp)
        coef_At = complex(coef_tp)
        
        return coef_Bref,coef_Bt,coef_Aref,coef_At

    def epsi(medio):
        if medio ==1:
            epsirta = epsi1
        elif medio ==2:
            epsirta = epsi2
        return epsirta 
        
    def mu(medio):
        if medio ==1:
            murta = mu1
        elif medio ==2:
            murta = mu2
        return murta 
            
    # def xt(medio):  # kz/k0
        
    #     def xt2(medio): # kz^2/k^0
    #         xt2rta = -x_x**2 + mu(medio)*epsi(medio)
    #         return xt2rta
        
    #     inside = xt2(medio)+0j
    #     xtmedio = (inside)**(1/2)   

    #     if (xtmedio*z_plane_barra).imag >= 0:
    # 	    xtmedio = xtmedio
    #     else:
    # 	    xtmedio = -xtmedio
    #     return xtmedio
    
    f1 = lambda qx, qy : qx*(np.e**(1j*(qx*x + qy*y)))/(np.sqrt(qx**2 + qy**2))    
#    f2 = lambda qx, qy : np.sqrt(qx**2 + qy**2)*np.e**(1j*(qx*x + qy*y))
    Int1 = integrate.dblquad(f1, -np.inf, np.inf, lambda qx: -np.inf, lambda x: np.inf)
#    Int2 = integrate.dblquad(f2, -np.inf, np.inf, lambda qx: -np.inf, lambda x: np.inf)
        
    coef_rs, coef_ts, coef_rp, coef_tp = coeficientes()

    term1 = complex(Int1)*(coef_rs + 1)
    term2 = complex(Int1)*(coef_rp + 1)*(kz1/k1)
    Ex = -term1 + term2
    
    return Ex*factor_final        
    