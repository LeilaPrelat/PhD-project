#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila
"""
import numpy as np
import sys
import os 
from scipy import special

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/potential_field/potential','')
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
    print('constants.py no se encuentra en ' + path_constants)

pi,hb,c,alfac,mu1,mu2 = constantes()
aux = c*hb


#%%

def electric_potential(omegac,epsi1,epsi2,hbmu,hbgama,x,y,z,list_xD,yD,zD,zp,px,py,pz):     
    """    
    Parameters
    ----------
    omegac : omega/c = k0 en 1/micrometros    
    epsi1 : epsilon del medio de arriba del plano
    epsi2 : epsilon del medio de abajo del plano
    hbmu : chemical potential in eV  
    hbgama : collision frequency in eV
    z : coordenada z
    list_xD : coordenada x del dipolo 
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
#    k0 = omegac #=omega/c
    # x_y = ky/k0
    n1 = epsi1*mu1
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1
#    k1_2 = (k0*cte1)**2
 #   n_v1 = int_v/cte1
    
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    Rp = 2*epsi1/(epsi1 + epsi2)
    alfa_p = 1j*(epsi1 + epsi2)/(cond*cte1)
    kp = alfa_p*k1
    kp_2 = kp**2 
    
    z_dip_barra = k1*(np.abs(z) + 2*zp + np.abs(zD))        
    exp_electron = np.exp(-alfa_p*z_dip_barra)
    
    
    potential0 = 0
    
    for xD in list_xD:
        
        phi_1 = np.arctan2(np.abs(y - yD),np.abs(x - xD))
        phi_2 = np.arctan2(np.abs(y + yD),np.abs(x + xD))
        
        R1 = np.sqrt((x - xD)**2 + (y - yD)**2)
        R2 = np.sqrt((x + xD)**2 + (y + yD)**2)
        
        term_px_py1 = px*np.cos(phi_1) + py*np.sin(phi_1)
        term_px_py2 = px*np.cos(phi_2) + py*np.sin(phi_2)
        
        term_aux1 = np.abs(z-zD)**2 + R1**2
        term_aux2 = (np.abs(z))**2 + R2**2
        
        term1 = -term_px_py1*(np.abs(z-zD)**2/(term_aux1**(3/2)) - 1/(term_aux1**(1/2)))/R1
        
        term2 = -term_px_py2*((np.abs(z))**2/(term_aux2**(3/2)) - 1/(term_aux2**(1/2)))/R1
        
        J1 = special.jv(1,kp*R2)
        J0 = special.jv(0,kp*R2)
        
        term3 = -2*np.pi*1j*Rp*kp_2*term_px_py2*J1*exp_electron
        
        term4 = np.sign(z)*pz*np.abs(z-zD)/(term_aux1**(3/2))
        
        term5 = np.sign(z)*pz*(np.abs(z) + np.abs(zD))/(term_aux2**(3/2))
        
        term6 = 2*np.pi*1j*np.sign(z)*pz*Rp*kp_2*J0*exp_electron
        
        ffinal =  term1 + term2 + term3 + term4 + term5 + term6
        
        potential0 = ffinal + potential0
    return potential0
 
 #%%