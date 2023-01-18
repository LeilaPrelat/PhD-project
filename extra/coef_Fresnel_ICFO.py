#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila

coeficientes de la matriz de 4x4
caso inhomogeneo (con campos incidentes)  
"""

import numpy as np
import sys
import os 

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_graphene =  path_basic.replace('/' + 'plane','')
#print('Importar modulos necesarios para este codigo')

try:
    sys.path.insert(1, path_graphene + '/' + 'sigma_graphene')
    from graphene_sigma import sigma
except ModuleNotFoundError:
    print('graphene_sigma.py no se encuentra en ' + path_graphene + '/' + 'sigma_graphene')

try:
    sys.path.insert(1, path_basic)
    from constants_plane import constantes
except ModuleNotFoundError:
    print('constants_plane.py no se encuentra en ' + path_basic)

pi,hb,c,alfac,hbargama,mu1,mu2,epsi1,epsi2 = constantes()
aux = c*hb

#%%

def coefTM_TE(omegac,hbmu,theta_inc,z_plane):
    """
    Parameters
    ----------
    omegac : omega/c = k0 en 1/micrometros    
    hbmu: potencial quimico del grafeno en eV
    theta_inc: theta incidente en radianes : kz incidente y kx incidente
    z_plane: radio del cilindro en micrometros
    Returns
    -------
    4 coef: refleado y trasmitido del modo TM,
    refleado y trasmitido del modo TE
    """
    
    E = omegac*aux
    k0 = omegac #=omega/c
    a_barra = z_plane*k0 #adimensional
 
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
    
    xz1 = np.sqrt(epsi1*mu1)*np.cos(theta_inc)
    aux0 = np.sqrt(epsi2*mu2)/np.sqrt(epsi1*mu1)
    aux1 = np.arcsin(aux0*np.sin(theta_inc))
    xz2 = np.sqrt(epsi2*mu2)*np.cos(aux1)
            
    # def xt(medio):
        
    #     def xt2(medio):
    #         xt2rta = -x_y**2 + mu(medio)*epsi(medio)
    #         return xt2rta
        
    #     inside = xt2(medio)+0j
    #     xtmedio = (inside)**(1/2)   

    #     if (xtmedio*k0).imag >= 0:
    # 	    xtmedio = xtmedio
    #     else:
    # 	    xtmedio = -xtmedio
    #     return xtmedio
    
    sigmatot = sigma(E,hbmu,hbargama)
    cond3 = sigmatot*alfac*4*pi  #sigma devuelve la conductividad teniendo que multiplicar por alfac*c ---> no hay que dividir por c
     
    # expAref = np.e**(2*1j*xt(1)*a_barra) ## signo cambiado ---> ver paper 370 ######
    # expAtra = np.e**(-1j*(xt(2)-xt(1))*a_barra)

    expAref = np.e**(2*1j*xz1*a_barra) ## signo cambiado ---> ver paper 370 ######
    expAtra = np.e**(-1j*(xz2-xz1)*a_barra)

    ### recordar : hay que intencarmbiar los indices 1 y 2 cuando se compara con el paper 370, 
    ## aca: el medio 1 es el medio arriba del plano y el medio 2 es el medio abajo del plano
    ## overleaf: el medio 2 es el medio arriba del plano y el medio 1 es el medio abajo del plano

# FORMULAS DEL OVERLEAF : 70a (r_s) y 70b (r_p)

    # term1_s = mu(2)*xt(1)
    # term2_s = mu(1)*xt(2)
    term1_s = mu(2)*xz1
    term2_s = mu(1)*xz2
    num_s = term1_s - term2_s - cond3  ## rs
    den_s = term1_s + term2_s + cond3  ## rs
    
    coef_rs = expAref*num_s/den_s    ## rs
    coef_ts = expAtra*2*term2_s/den_s  ## ts
    
    # term1_p = epsi(2)*xt(1)
    # term2_p = epsi(1)*xt(2)
    term1_p = epsi(2)*xz1
    term2_p = epsi(1)*xz2
    num_p = term1_p - term2_p + cond3*xz2*xz1 ## rp
    den_p = term1_p + term2_p + cond3*xz2*xz1  ## rp
    
    coef_rp = expAref*num_p/den_p
    coef_tp = expAtra*2*term2_p/den_p
    
    return coef_rs, coef_ts, coef_rp, coef_tp

#%%
