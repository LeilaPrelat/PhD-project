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

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/potential_field/potential_and_electric_field/sum_dipoles/potential_and_electric_field_simple_version','')
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
limitt = 60

#%%

def potential_final_ana(omegac,epsi1,epsi2,hbmu,hbgama,x,y,z,a,zp,int_v,px,py,pz,Ndip):
    
    
    rta0 = 0
    
    for j in range(-Ndip,Ndip + 1):
    
        expo_j = np.exp(1j*omegac*int_v*a*j)
    
        phi_j = np.arctan2(y,x-a*j)
        
        Rj = np.sqrt((x-a*j)**2 + y**2)
        
        px_py_term = px*np.cos(phi_j) + py*np.sin(phi_j)
    
        E = omegac*aux
        cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
        Rp = 2*epsi1/(epsi1 + epsi2)
        alfa_p = 1j*(epsi1 + epsi2)/(cond)
        kp = alfa_p*omegac    
        kp_2 = kp**2
        
        
        expo = np.exp(-kp*(2*zp-z))
        
        H0 = special.hankel1(0,kp*Rj)
        H1 = special.hankel1(1,kp*Rj)
        
        rta = - (px_py_term*H1 + pz*np.sign(z)*H0 )*pi*1j*Rp*kp_2*expo*expo_j
    
        rta0 = rta0 + rta
    
    return rta

#%%
