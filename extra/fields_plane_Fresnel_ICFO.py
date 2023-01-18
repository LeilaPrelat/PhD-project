#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 19 20:56:58 2020

@author: leila

campos de un plano

"""
import numpy as np
import os 
import sys

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
    from constants_plane import constantes
except ModuleNotFoundError:
    print('constants_plane.py no se encuentra en ' + path_basic)

pi,hb,c,alfac,hbargama,mu1,mu2,epsi1,epsi2 = constantes()
aux = c*hb

#%%
    
#print('Definir la funcion Hz para el medio 1 y el 2')

def fieldE_plane(omegac,hbmu,theta_inc,z_plane,x,z,Ainc,Binc): 
    """
    Parameters
    ----------
    omegac : omega/c = k0 en 1/micrometros    
    hbmu: potencial quimico del grafeno en eV
    theta_inc: theta incidente en radianes
    z_plane: posicion del plano en micrometros (z_plane<0)
    x: x en micrometros
    z: z en micrometros
    Ainc: amplitud del campo incidente (medio 1 campo electrico)
    Binc: amplitud del campo incidente (medio 1 campo magnetico)
    Returns
        Ex, Ey, Ez de un plano
    -------
    """
    k0 = omegac #=omega/c
    z_barra = z*k0
    z_plane_barra = z_plane*k0
    # x_x = kx/k0
    
    def coeficientes():
        
        coeff = coefTM_TE(omegac,hbmu,theta_inc,z_plane)
        Bref,Bt,Aref,At = coeff  
        
        coef_Aref = complex(Aref)
        coef_At = complex(At)
        coef_Bref = complex(Bref)
        coef_Bt = complex(Bt)
        
        return coef_Bref*Binc,coef_Bt*Binc,coef_Aref*Ainc,coef_At*Ainc

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
            
    def xt(medio):  # kz/k0
        
        def xt2(medio): # kz^2/k^0
            xt2rta = -x_x**2 + mu(medio)*epsi(medio)
            return xt2rta
        
        inside = xt2(medio)+0j
        xtmedio = (inside)**(1/2)   

        if (xtmedio*z_plane_barra).imag >= 0:
    	    xtmedio = xtmedio
        else:
    	    xtmedio = -xtmedio
        return xtmedio
    
    exp_x = np.e**(1j*kx*x)
    exp1_mas = np.e**(1j*xt(1)*z_barra)
    exp1_menos = np.e**(-1j*xt(1)*z_barra)
    # exp2_mas = np.e**(1j*xt(2)*z_barra)
    exp2_mas = np.e**(1j*xt(2)*z_barra)        
    
    coef_Bref,coef_Bt,coef_Aref,coef_At = coeficientes()
    if z >= z_plane: #medio 1

        Ey = (mu(1)/xt(1))*(Binc*exp1_mas - coef_Bref*exp1_menos)
        Ex = Ainc*exp1_mas + coef_Aref*exp1_menos
        Ez = (-Ainc*exp1_mas + coef_Aref*exp1_menos)*x_x/xt(1)
    
        return np.array([Ex,Ey,Ez])*exp_x
    else:
                    
        Ey = coef_Bt*exp2_mas*mu(2)/xt(2)
        Ex = coef_At*exp2_mas
        Ez = -coef_At*exp2_mas*x_x/xt(2)
    
        return np.array([Ex,Ey,Ez])*exp_x        
    
#%%        

#print('Definir la funcion Ez para el medio 1 y el 2')

def fieldH_plane(omegac,hbmu,kx,z_plane,x,z,Ainc,Binc): 
    """
    Parameters
    ----------
    omegac : omega/c = k0 en 1/micrometros    
    hbmu: potencial quimico del grafeno en eV
    ky: ky en 1/micrometros
    z_plane: posicion del plano en micrometros (z_plane<0)
    x: x en micrometros
    y: y en micrometros
    z: z en micrometros
    Ainc: amplitud del campo incidente (medio 1 campo electrico)
    Binc: amplitud del campo incidente (medio 1 campo magnetico)
    Returns
        Hx, Hy, Hz de un plano
    -------
    """
    k0 = omegac #=omega/c
    z_barra = z*k0
    z_plane_barra = z_plane*k0
    x_x = kx/k0
    
    def coeficientes():
        
        coeff = coefTM_TE(omegac,hbmu,kx,z_plane)
        Bref,Bt,Aref,At = coeff  
        
        coef_Aref = complex(Aref)
        coef_At = complex(At)
        coef_Bref = complex(Bref)
        coef_Bt = complex(Bt)
        
        return coef_Bref*Binc,coef_Bt*Binc,coef_Aref*Ainc,coef_At*Ainc

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
            
    def xt(medio):  # kz/k0
        
        def xt2(medio): # kz^2/k^0
            xt2rta = -x_x**2 + mu(medio)*epsi(medio)
            return xt2rta
        
        inside = xt2(medio)+0j
        xtmedio = (inside)**(1/2)   

        if (xtmedio*z_plane_barra).imag >= 0:
    	    xtmedio = xtmedio
        else:
    	    xtmedio = -xtmedio
        return xtmedio
   
    exp_x = np.e**(1j*kx*x)
    exp1_mas = np.e**(1j*xt(1)*z_barra)
    exp1_menos = np.e**(-1j*xt(1)*z_barra)
    # exp2_mas = np.e**(1j*xt(2)*z_barra)
    exp2_mas = np.e**(1j*xt(2)*z_barra)        

    coef_Bref,coef_Bt,coef_Aref,coef_At = coeficientes()    
    
    if z>= z_plane: #medio 1
        
        Hy = (epsi(1)/xt(1))*(-Ainc*exp1_mas + coef_Aref*exp1_menos)
        Hx = Binc*exp1_mas + coef_Bref*exp1_menos
        Hz = (-Binc*exp1_mas + coef_Bref*exp1_menos)*x_x/xt(1)
    
        return np.array([Hx,Hy,Hz])*exp_x
    else:
            
        Hy = -coef_At*exp2_mas*epsi(2)/xt(2)
        Hx = coef_Bt*exp2_mas
        Hz = -coef_Bt*exp2_mas*x_x/xt(2)
    
        return np.array([Hx,Hy,Hz])*exp_x        
    
#%%

