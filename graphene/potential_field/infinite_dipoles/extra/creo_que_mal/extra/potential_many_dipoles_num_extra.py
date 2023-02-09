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

def phi_many_dipoles_num_extra(omegac,epsi1,epsi2,hbmu,hbgama,x,y,z,a,N,zp,int_v,ky,px,py,pz):
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
        
    def function(ky,n):
#        print(n)
        kx = omegac*c/v - 2*pi*n/a
        
        k_parallel = np.sqrt(kx**2 + ky**2)
        exp_z = np.exp(-k_parallel*np.abs(z))
        exp_xy = np.exp(1j*ky*y + 1j*kx*x)
        rp = Rp*kp/(k_parallel - kp)
    ################### I1 ################################################
        if px!= 0 or py!= 0:
            int1 = -1j*(px*kx + ky*py)*exp_z*exp_xy/k_parallel
        else:
            int1 = 0
    ################### I2 ###############################################
        if pz != 0:
            int2 = pz*exp_z*ky*exp_xy
        else:
            int2 = 0
    ################### I3 ###############################################
        if px != 0:
            int3 = px*1j*exp_z*rp*exp_xy/k_parallel
        else:
            int3 = 0
    
    ################### I4 ###############################################
        if py != 0:
            int4 = py*1j*exp_z*ky*rp*exp_xy/k_parallel
        else: 
            int4 = 0 
    ################### I5 ###############################################
        if pz!= 0:   
            int5 = -pz*exp_z*ky*rp*exp_xy
        else:
            int5 = 0
    ################### I6 ###############################################    
        
        int_tot = int1 + int2 + int3 + int4 + int5
        
        return a*int_tot
    
    
    int_tot0 = 0  
    N = int(N)
    for n in range(N):
        int_tot0 = function(ky,n) + int_tot0
    
    
    return int_tot0

#%%
