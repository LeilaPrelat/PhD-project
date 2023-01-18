#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 19 20:56:58 2020

@author: leila

relacion de dispersion
solucion analitica
para un plano de Ag

"""
import os 
import sys

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.optimize import fsolve 

sns.set()

save_plots = 1
save_txt = 0

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_ctes =  path_basic.replace('/' + 'new','')
#print('Importar modulos necesarios para este codigo')
path_save = path_basic + '/' + 'disp_relation_Silica'
try:
    sys.path.insert(1, path_ctes)
    from constants import constantes
except ModuleNotFoundError:
    print('constants.py no se encuentra en ' + path_ctes)

pi,hb,c,alfac,mu1,mu2,mu3 = constantes()
aux = c*hb

 
epsilon2 = 12
d = 20*1e-3

#%%

def k_parallel_air(omega):

    return omega/c

def k_parallel_medium(omega,epsilon2):

    return omega*np.sqrt(epsilon2)/c

#%%
    
#def k_parallel_WG_TM1(omega,epsilon2,d,k_parallel):
#    k0 = omega/c
#    
#    if k_parallel**2 <= epsilon2*(k0**2):    
#        chi_prima = d*np.sqrt(epsilon2*k0**2 - k_parallel**2)
#    else:
#        chi_prima = d*1j*np.sqrt(k_parallel**2 - epsilon2*k0**2)
#
#    if k_parallel**2 <= k0**2 :
#        chi = d*1j*np.sqrt(k0**2 - k_parallel**2)
#    else:
#        chi = d*np.sqrt(k_parallel**2 - k0**2)
#    
#    return np.tan(chi_prima/2) - epsilon2*chi/chi_prima



def k_parallel_WG_TM2A(omega,epsilon2,d,k_parallel):
    k0 = omega/c
    
    if k_parallel**2 <= epsilon2*(k0**2):    
        chi_prima = d*np.sqrt(epsilon2*k0**2 - k_parallel**2)
    else:
        chi_prima = d*1j*np.sqrt(k_parallel**2 - epsilon2*k0**2)

    if k_parallel**2 <= k0**2 :
        chi = d*1j*np.sqrt(k0**2 - k_parallel**2)
    else:
        chi = d*np.sqrt(k_parallel**2 - k0**2)
    
    return chi_prima/2 + np.arctan(chi_prima/(epsilon2*chi))

def k_parallel_WG_TM2B(omega,epsilon2,d,k_parallel):
    k0 = omega/c
    
    if k_parallel**2 <= epsilon2*(k0**2):    
        chi_prima = d*np.sqrt(epsilon2*k0**2 - k_parallel**2)
    else:
        chi_prima = d*1j*np.sqrt(k_parallel**2 - epsilon2*k0**2)

    if k_parallel**2 <= k0**2 :
        chi = d*1j*np.sqrt(k0**2 - k_parallel**2)
    else:
        chi = d*np.sqrt(k_parallel**2 - k0**2)
    
    return 1/np.tan(chi_prima/2) + epsilon2*chi/chi_prima


def k_parallel_WG_TM1A(omega,epsilon2,d,k_parallel):
    k0 = omega/c
    
    if k_parallel**2 <= epsilon2*(k0**2):    
        chi_prima = d*np.sqrt(epsilon2*k0**2 - k_parallel**2)
    else:
        chi_prima = d*1j*np.sqrt(k_parallel**2 - epsilon2*k0**2)

    if k_parallel**2 <= k0**2 :
        chi = d*1j*np.sqrt(k0**2 - k_parallel**2)
    else:
        chi = d*np.sqrt(k_parallel**2 - k0**2)
    
    return chi_prima/2 - np.arctan(epsilon2*chi/chi_prima)

def k_parallel_WG_TM1B(omega,epsilon2,d,k_parallel):
    k0 = omega/c
    
    if k_parallel**2 <= epsilon2*(k0**2):    
        chi_prima = d*np.sqrt(epsilon2*k0**2 - k_parallel**2)
    else:
        chi_prima = d*1j*np.sqrt(k_parallel**2 - epsilon2*k0**2)

    if k_parallel**2 <= k0**2 :
        chi = d*1j*np.sqrt(k0**2 - k_parallel**2)
    else:
        chi = d*np.sqrt(k_parallel**2 - k0**2)
    
    return np.tan(chi_prima/2) - epsilon2*chi/chi_prima

#%%
    

def k_parallel_WG_TE1A(omega,epsilon2,d,k_parallel):
    k0 = omega/c
    
    if k_parallel**2 <= epsilon2*(k0**2):    
        chi_prima = d*np.sqrt(epsilon2*k0**2 - k_parallel**2)
    else:
        chi_prima = d*1j*np.sqrt(k_parallel**2 - epsilon2*k0**2)

    if k_parallel**2 <= k0**2 :
        chi = d*1j*np.sqrt(k0**2 - k_parallel**2)
    else:
        chi = d*np.sqrt(k_parallel**2 - k0**2)
    
    
    return chi_prima/2 - np.arctan(chi/chi_prima)

def k_parallel_WG_TE1B(omega,epsilon2,d,k_parallel):
    k0 = omega/c
    
    if k_parallel**2 <= epsilon2*(k0**2):    
        chi_prima = d*np.sqrt(epsilon2*k0**2 - k_parallel**2)
    else:
        chi_prima = d*1j*np.sqrt(k_parallel**2 - epsilon2*k0**2)

    if k_parallel**2 <= k0**2 :
        chi = d*1j*np.sqrt(k0**2 - k_parallel**2)
    else:
        chi = d*np.sqrt(k_parallel**2 - k0**2)
    
    
    return np.tan(chi_prima/2) - chi/chi_prima


#def k_parallel_WG_TE1C(omega,epsilon2,d,k_parallel):
#    k0 = omega/c
#    
#    if k_parallel**2 <= epsilon2*(k0**2):    
#        chi_prima = d*np.sqrt(epsilon2*k0**2 - k_parallel**2)
#    else:
#        chi_prima = d*1j*np.sqrt(k_parallel**2 - epsilon2*k0**2)
#
#    if k_parallel**2 <= k0**2 :
#        chi = d*1j*np.sqrt(k0**2 - k_parallel**2)
#    else:
#        chi = d*np.sqrt(k_parallel**2 - k0**2)
#    
#    rta = np.arctan(chi/chi_prima)
#    rta[rta<0] +=  2*np.pi
#    
#    return chi_prima/2 - rta  ## OTRA RAMA


#%%

def k_parallel_WG_TE2A(omega,epsilon2,d,k_parallel):
    k0 = omega/c
    
    if k_parallel**2 <= epsilon2*(k0**2):    
        chi_prima = d*np.sqrt(epsilon2*k0**2 - k_parallel**2)
    else:
        chi_prima = d*1j*np.sqrt(k_parallel**2 - epsilon2*k0**2)

    if k_parallel**2 <= k0**2 :
        chi = d*1j*np.sqrt(k0**2 - k_parallel**2)
    else:
        chi = d*np.sqrt(k_parallel**2 - k0**2)
    
    
    rta = np.arctan(chi_prima/chi - 2*np.pi)
#    rta[rta<0] +=  2*np.pi
    return chi_prima/2 + rta  ## primera rama 

def k_parallel_WG_TE2C(omega,epsilon2,d,k_parallel):
    k0 = omega/c
    
    if k_parallel**2 <= epsilon2*(k0**2):    
        chi_prima = d*np.sqrt(epsilon2*k0**2 - k_parallel**2)
    else:
        chi_prima = d*1j*np.sqrt(k_parallel**2 - epsilon2*k0**2)

    if k_parallel**2 <= k0**2 :
        chi = d*1j*np.sqrt(k0**2 - k_parallel**2)
    else:
        chi = d*np.sqrt(k_parallel**2 - k0**2)
    
    return chi_prima/2 + np.arctan2(1,chi/chi_prima)  ## segunda rama 

0
#def k_parallel_WG_TE2C(omega,epsilon2,d,k_parallel):
#    k0 = omega/c
#    
#    if k_parallel**2 <= epsilon2*(k0**2):    
#        chi_prima = d*np.sqrt(epsilon2*k0**2 - k_parallel**2)
#    else:
#        chi_prima = d*1j*np.sqrt(k_parallel**2 - epsilon2*k0**2)
#
#    if k_parallel**2 <= k0**2 :
#        chi = d*1j*np.sqrt(k0**2 - k_parallel**2)
#    else:
#        chi = d*np.sqrt(k_parallel**2 - k0**2)
#    
#    return chi_prima/2 + np.arctan(chi_prima/chi) + 2*np.pi  ## 3era rama 

def k_parallel_WG_TE2B(omega,epsilon2,d,k_parallel):
    k0 = omega/c
    
    if k_parallel**2 <= epsilon2*(k0**2):    
        chi_prima = d*np.sqrt(epsilon2*k0**2 - k_parallel**2)
    else:
        chi_prima = d*1j*np.sqrt(k_parallel**2 - epsilon2*k0**2)

    if k_parallel**2 <= k0**2 :
        chi = d*1j*np.sqrt(k0**2 - k_parallel**2)
    else:
        chi = d*np.sqrt(k_parallel**2 - k0**2)
    
    return 1/np.tan(chi_prima/2) + chi/chi_prima

def k_parallel_WG_TE2C(omega,epsilon2,d,k_parallel):
    k0 = omega/c
    
    if k_parallel**2 <= epsilon2*(k0**2):    
        chi_prima = d*np.sqrt(epsilon2*k0**2 - k_parallel**2)
    else:
        chi_prima = d*1j*np.sqrt(k_parallel**2 - epsilon2*k0**2)

    if k_parallel**2 <= k0**2 :
        chi = d*1j*np.sqrt(k0**2 - k_parallel**2)
    else:
        chi = d*np.sqrt(k_parallel**2 - k0**2)
    
    return  1/np.tan(chi_prima/2 + np.pi/2) + chi/chi_prima  ## segunda rama 

#%%
    
tamfig = (4.5,3.5)
tamlegend = 10
tamletra = 12
tamtitle = 10
tamnum = 9
loc2 = [0,1]
pad = -2
lw = 1.5
hp = 0.3
mk = 2
labelpady = 0
labelpadx = 0
ms = 3

def graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad):
    plt.figure(figsize=tamfig)
    plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
    plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
    plt.tick_params(labelsize = tamnum, pad = pad)
    return   

#%%


title = '$\epsilon_2$ = %i, d = %inm' %(epsilon2,d*1e3)
    
cota_inf_omega_omegaWG = 0.00001
cota_sup_omega_omegaWG = 4
N = 800
list_omega_omegaWGA = np.linspace(cota_inf_omega_omegaWG,cota_sup_omega_omegaWG,N)
list_omega_omegaWG_TM2B_array = np.linspace(2.3,4,N)
list_omega_omegaWG_TM1B_array = np.linspace(2,4,N)
list_omega_omegaWG_TE1B_array = np.linspace(1.8,4,N)
list_omega_omegaWG_TE2A_array = np.linspace(0.99,4,N)
list_omega_omegaWG_TE2B_array = np.linspace(2.79,4,N)

#list_omega_omegaWG = np.linspace(4,0.005,N)
omegaWG = (np.pi*c/d)/(np.sqrt(epsilon2-1)) 
cte_new = c/omegaWG
        
list_k_parallel = np.linspace(0.001,10,N)

#%%

list_y_air = []
list_y_medium = []

zeros_eqTM1A = []
init_condit_TM1A = []
list_omega_omegaWG_TM1A = []

k_TM1 = 0
for omega_omegaWG in list_omega_omegaWGA:
    omega = omega_omegaWG*omegaWG
    k_air = k_parallel_air(omega)
    k_medium = k_parallel_medium(omega,epsilon2)   
    
#    print(k_air)
    
    list_y_air.append(k_air*cte_new)
    list_y_medium.append(k_medium*cte_new)

    init_condit = [k_air]
            
    def k_parallel_WG_TM1_1var(k_parallel):   
        return  np.abs(k_parallel_WG_TM1A(omega,epsilon2,d,k_parallel))
    
    if k_TM1 == 0 :
        resTM1 = fsolve(k_parallel_WG_TM1_1var, [k_air])
    else:
        resTM1 = fsolve(k_parallel_WG_TM1_1var, [init_condit_TM1A[k_TM1-1]],maxfev = 1000 )         
    resTM1_v = np.float(resTM1)
    if k_air < resTM1_v < k_medium:
        zeros_eqTM1A.append(resTM1_v)
        list_omega_omegaWG_TM1A.append(omega_omegaWG)
        k_TM1 = k_TM1 + 1
    #            if resTM1_v*cte_new < 7:
    #                init_condit_TM1A.append(resTM1_v)
    #            else:
    #                init_condit_TM1A.append(k_medium)
        
#        if omega_omegaWG < 1.6:
#            init_condit_TM1A.append(resTM1_v)
#        else:
#            init_condit_TM1A.append(k_medium)

        init_condit_TM1A.append(resTM1_v)



zeros_eqTM1B = []
init_condit_TM1B = []
list_omega_omegaWG_TM1B = []

k_TM1 = 0
for omega_omegaWG in list_omega_omegaWG_TM1B_array:
    omega = omega_omegaWG*omegaWG
    k_air = k_parallel_air(omega)
    k_medium = k_parallel_medium(omega,epsilon2)   
    
    init_condit = [k_air]
            
    def k_parallel_WG_TM1_1var(k_parallel):   
        return  np.abs(k_parallel_WG_TM1B(omega,epsilon2,d,k_parallel))
    
    if k_TM1 == 0 :
        resTM1 = fsolve(k_parallel_WG_TM1_1var, [k_air])
    else:
        resTM1 = fsolve(k_parallel_WG_TM1_1var, [init_condit_TM1B[k_TM1-1]],maxfev = 1000 )         
    resTM1_v = np.float(resTM1)
    if k_air < resTM1_v < k_medium:
        zeros_eqTM1B.append(resTM1_v)
        list_omega_omegaWG_TM1B.append(omega_omegaWG)
        k_TM1 = k_TM1 + 1
    #            if resTM1_v*cte_new < 7:
    #                init_condit_TM1A.append(resTM1_v)
    #            else:
    #                init_condit_TM1A.append(k_medium)
        
#        if omega_omegaWG < 1.6:
#            init_condit_TM1A.append(resTM1_v)
#        else:
#            init_condit_TM1A.append(k_medium)

        init_condit_TM1B.append(resTM1_v)
        
#%%
        
zeros_eqTM2A = []
init_condit_TM2A = []
list_omega_omegaWG_TM2A = []


#list_omega_omegaWGB = np.linspace(2,cota_sup_omega_omegaWG,N)

k_TM2 = 0
for omega_omegaWG in list_omega_omegaWGA:
    omega = omega_omegaWG*omegaWG
    k_air = k_parallel_air(omega)

    k_medium = k_parallel_medium(omega,epsilon2)   

    init_condit = [k_air]
            
    def k_parallel_WG_TM2_1var(k_parallel):   
        return  np.abs(k_parallel_WG_TM2B(omega,epsilon2,d,k_parallel))
    
    if k_TM2 == 0 :
        resTM2 = fsolve(k_parallel_WG_TM2_1var, [k_air])
    else:
        resTM2 = fsolve(k_parallel_WG_TM2_1var, [init_condit_TM2A[k_TM2-1]],maxfev = 1000 )         
    resTM2_v = np.float(resTM2)
    if k_air < resTM2_v < k_medium:
        zeros_eqTM2A.append(resTM2_v)
        list_omega_omegaWG_TM2A.append(omega_omegaWG)
        k_TM2 = k_TM2 + 1
    #            if resTM1_v*cte_new < 7:
    #                init_condit_TM1A.append(resTM1_v)
    #            else:
    #                init_condit_TM1A.append(k_medium)
        
#        if omega_omegaWG > 3:
#            init_condit_TM2A.append(resTM2_v)
#        else:
#            init_condit_TM2A.append(k_medium)

        init_condit_TM2A.append(resTM2_v)


zeros_eqTM2B = []
init_condit_TM2B = []
list_omega_omegaWG_TM2B = []
#list_omega_omegaWGB = np.linspace(2,cota_sup_omega_omegaWG,N)

k_TM2 = 0
for omega_omegaWG in list_omega_omegaWG_TM2B_array:
    omega = omega_omegaWG*omegaWG
    k_air = k_parallel_air(omega)

    k_medium = k_parallel_medium(omega,epsilon2)   

    init_condit = [k_air]
            
    def k_parallel_WG_TM2_1var(k_parallel):   
        return  np.abs(k_parallel_WG_TM2B(omega,epsilon2,d,k_parallel))
    
    if k_TM2 == 0 :
        resTM2 = fsolve(k_parallel_WG_TM2_1var, [k_air])
    else:
        resTM2 = fsolve(k_parallel_WG_TM2_1var, [init_condit_TM2B[k_TM2-1]],maxfev = 1000 )         
    resTM2_v = np.float(resTM2)
    if k_air < resTM2_v < k_medium:
        zeros_eqTM2B.append(resTM2_v)
        list_omega_omegaWG_TM2B.append(omega_omegaWG)
        k_TM2 = k_TM2 + 1
    #            if resTM1_v*cte_new < 7:
    #                init_condit_TM1A.append(resTM1_v)
    #            else:
    #                init_condit_TM1A.append(k_medium)
        
#        if omega_omegaWG > 3:
#            init_condit_TM2A.append(resTM2_v)
#        else:
#            init_condit_TM2A.append(k_medium)

        init_condit_TM2B.append(resTM2_v)

del init_condit_TM2B,init_condit_TM2A,init_condit_TM1A,init_condit_TM1B,resTM2,resTM1,k_TM2,k_TM1

del k_parallel_WG_TM2_1var,k_parallel_WG_TM1_1var

#%%


zeros_eqTE1A = []
init_condit_TE1A = []
list_omega_omegaWG_TE1A = []

k_TE1 = 0
for omega_omegaWG in list_omega_omegaWGA:
    omega = omega_omegaWG*omegaWG
    k_air = k_parallel_air(omega)
    k_medium = k_parallel_medium(omega,epsilon2)   
    
    init_condit = [k_air]
            
    def k_parallel_WG_TE1_1var(k_parallel):   
        return  np.abs(k_parallel_WG_TE1A(omega,epsilon2,d,k_parallel))
    
    if k_TE1 == 0 :
        resTE1 = fsolve(k_parallel_WG_TE1_1var, [k_air])
    else:
        resTE1 = fsolve(k_parallel_WG_TE1_1var, [init_condit_TE1A[k_TE1-1]],maxfev = 1000 )         
    resTE1_v = np.float(resTE1)
    if k_air < resTE1_v < k_medium:
        zeros_eqTE1A.append(resTE1_v)
        list_omega_omegaWG_TE1A.append(omega_omegaWG)
        k_TE1 = k_TE1 + 1
    #            if resTM1_v*cte_new < 7:
    #                init_condit_TM1A.append(resTM1_v)
    #            else:
    #                init_condit_TM1A.append(k_medium)
        
#        if omega_omegaWG < 1.6:
#            init_condit_TM1A.append(resTM1_v)
#        else:
#            init_condit_TM1A.append(k_medium)

        init_condit_TE1A.append(resTE1_v)
        


zeros_eqTE1B = []
init_condit_TE1B = []
list_omega_omegaWG_TE1B = []

k_TE1 = 0
#list_omega_omegaWG_TE1B_array = list_omega_omegaWG_TE1B_array[::-1]
for omega_omegaWG in list_omega_omegaWG_TE1B_array:
    omega = omega_omegaWG*omegaWG
    k_air = k_parallel_air(omega)
    k_medium = k_parallel_medium(omega,epsilon2)   
    
#    print(k_air)
    
    init_condit = [k_air]
            
    def k_parallel_WG_TE1_1var(k_parallel):   
        return  np.abs(k_parallel_WG_TE1B(omega,epsilon2,d,k_parallel))
    
    if k_TE1 == 0 :
        resTE1 = fsolve(k_parallel_WG_TE1_1var, [k_air])
    else:
        resTE1 = fsolve(k_parallel_WG_TE1_1var, [init_condit_TE1B[k_TE1-1]],maxfev = 1000 )         
    resTE1_v = np.float(resTE1)
    if k_air < resTE1_v < k_medium:
        zeros_eqTE1B.append(resTE1_v)
        list_omega_omegaWG_TE1B.append(omega_omegaWG)
        k_TE1 = k_TE1 + 1
    #            if resTM1_v*cte_new < 7:
    #                init_condit_TM1A.append(resTM1_v)
    #            else:
    #                init_condit_TM1A.append(k_medium)
        
#        if omega_omegaWG < 1.6:
#            init_condit_TM1A.append(resTM1_v)
#        else:
#            init_condit_TM1A.append(k_medium)

        init_condit_TE1B.append(resTE1_v)

#zeros_eqTE1B = zeros_eqTE1B[::-1]
#list_omega_omegaWG_TE1B = list_omega_omegaWG_TE1B[::-1]



#%%


zeros_eqTE2A = []
init_condit_TE2A = []
list_omega_omegaWG_TE2A = []

k_TE2 = 0
for omega_omegaWG in list_omega_omegaWG_TE2A_array:
    omega = omega_omegaWG*omegaWG
    k_air = k_parallel_air(omega)
    k_medium = k_parallel_medium(omega,epsilon2)   
    
    init_condit = [k_air]
            
    def k_parallel_WG_TE2_1var(k_parallel):   
        return  np.abs(k_parallel_WG_TE2B(omega,epsilon2,d,k_parallel))
    
    if k_TE2 == 0 :
        resTE2 = fsolve(k_parallel_WG_TE2_1var, [k_air])
    else:
        resTE2 = fsolve(k_parallel_WG_TE2_1var, [init_condit_TE2A[k_TE2-1]],maxfev = 1000 )         
    resTE2_v = np.float(resTE2)
    if k_air < resTE2_v < k_medium:
        zeros_eqTE2A.append(resTE2_v)
        list_omega_omegaWG_TE2A.append(omega_omegaWG)
        k_TE2 = k_TE2 + 1
    #            if resTM1_v*cte_new < 7:
    #                init_condit_TM1A.append(resTM1_v)
    #            else:
    #                init_condit_TM1A.append(k_medium)
        
#        if omega_omegaWG < 1.6:
#            init_condit_TM1A.append(resTM1_v)
#        else:
#            init_condit_TM1A.append(k_medium)

        init_condit_TE2A.append(resTE2_v)
        


zeros_eqTE2B = []
init_condit_TE2B = []
list_omega_omegaWG_TE2B = []

k_TE2 = 0
#list_omega_omegaWG_TE1B_array = list_omega_omegaWG_TE1B_array[::-1]
for omega_omegaWG in list_omega_omegaWG_TE2B_array:
    omega = omega_omegaWG*omegaWG
    k_air = k_parallel_air(omega)
    k_medium = k_parallel_medium(omega,epsilon2)   
    
#    print(k_air)
    
    init_condit = [k_air]
            
    def k_parallel_WG_TE2_1var(k_parallel):   
        return  np.abs(k_parallel_WG_TE2B(omega,epsilon2,d,k_parallel))
    
    if k_TE2 == 0 :
        resTE2 = fsolve(k_parallel_WG_TE2_1var, [k_air])
    else:
        resTE2 = fsolve(k_parallel_WG_TE2_1var, [init_condit_TE2B[k_TE2-1]],maxfev = 1000 )         
    resTE2_v = np.float(resTE2)
    if k_air < resTE2_v < k_medium:
        zeros_eqTE2B.append(resTE2_v)
        list_omega_omegaWG_TE2B.append(omega_omegaWG)
        k_TE2 = k_TE2 + 1
    #            if resTM1_v*cte_new < 7:
    #                init_condit_TM1A.append(resTM1_v)
    #            else:
    #                init_condit_TM1A.append(k_medium)
        
#        if omega_omegaWG < 1.6:
#            init_condit_TM1A.append(resTM1_v)
#        else:
#            init_condit_TM1A.append(k_medium)

        init_condit_TE2B.append(resTE2_v)


#%%

#
plt.figure(figsize=tamfig)    
plt.title(title,fontsize=tamtitle)
plt.plot(np.array(zeros_eqTM1A)*cte_new,list_omega_omegaWG_TM1A,'.',color = 'lightblue',ms = 4,label = 'TM1')
plt.plot(np.array(zeros_eqTM1B)*cte_new,list_omega_omegaWG_TM1B,'.',color = 'lightblue',ms = 4)
plt.plot(np.array(zeros_eqTM2A)*cte_new,list_omega_omegaWG_TM2A,'.',color = 'blue',ms = 4,label = 'TM2')
plt.plot(np.array(zeros_eqTM2B)*cte_new,list_omega_omegaWG_TM2B,'.',color = 'blue',ms = 4)

plt.plot(np.array(zeros_eqTE1A)*cte_new,list_omega_omegaWG_TE1A,'.',color = 'red',ms = 4,label = 'TE1')
plt.plot(np.array(zeros_eqTE1B)*cte_new,list_omega_omegaWG_TE1B,'--',color = 'red',ms = 4)

plt.plot(np.array(zeros_eqTE2A)*cte_new,list_omega_omegaWG_TE2A,'.',color = 'darkred',ms = 4,label = 'TE2')
plt.plot(np.array(zeros_eqTE2B)*cte_new,list_omega_omegaWG_TE2B,'--',color = 'darkred',ms = 4)


plt.xlabel(r'$k_\parallel$c/$\omega_{WG}$',fontsize=tamletra, labelpad = 0)
plt.ylabel(r'$\omega/\omega_{WG}$',fontsize=tamletra, labelpad = 0)
plt.plot(list_y_air,list_omega_omegaWGA,'--',color = 'grey',ms = lw, label = 'air')
plt.plot(list_y_medium,list_omega_omegaWGA,'--',color = 'black',ms = lw, label = 'medium')
plt.xlim([-0.5,8.5])
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
plt.grid(1)
plt.tight_layout()
if save_plots == 1 :
    plt.savefig('disp_relation_Silicon_totv2_d%inm.png' %(d*1e3))

#%%

plt.figure(figsize=tamfig)    
plt.title(title,fontsize=tamtitle)
plt.plot(np.array(zeros_eqTM1A)*cte_new,list_omega_omegaWG_TM1A,'.',color = 'blue',ms = ms,label = 'TM')
plt.plot(np.array(zeros_eqTM1B)*cte_new,list_omega_omegaWG_TM1B,'.',color = 'blue',ms = ms)
plt.plot(np.array(zeros_eqTM2A)*cte_new,list_omega_omegaWG_TM2A,'.',color = 'blue',ms = ms)
plt.plot(np.array(zeros_eqTM2B)*cte_new,list_omega_omegaWG_TM2B,'.',color = 'blue',ms = ms)

plt.plot(np.array(zeros_eqTE1A)*cte_new,list_omega_omegaWG_TE1A,'.',color = 'red',ms = ms,label = 'TE')
plt.plot(np.array(zeros_eqTE1B)*cte_new,list_omega_omegaWG_TE1B,'.',color = 'red',ms = ms)

plt.plot(np.array(zeros_eqTE2A)*cte_new,list_omega_omegaWG_TE2A,'.',color = 'red',ms = ms)
plt.plot(np.array(zeros_eqTE2B)*cte_new,list_omega_omegaWG_TE2B,'.',color = 'red',ms = ms)


plt.xlabel(r'$k_\parallel$c/$\omega_{WG}$',fontsize=tamletra, labelpad = 0)
plt.ylabel(r'$\omega/\omega_{WG}$',fontsize=tamletra, labelpad = 0)
plt.plot(list_y_air,list_omega_omegaWGA,'--',color = 'grey',ms = lw, label = 'air')
plt.plot(list_y_medium,list_omega_omegaWGA,'--',color = 'black',ms = lw, label = 'medium')
plt.xlim([-0.5,8.5])
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
plt.grid(1)
plt.tight_layout()
if save_plots == 1 :
    plt.savefig('disp_relation_Silicon_tot_d%inm.png' %(d*1e3))
## os.chdir(path_save)
#plt.savefig('disp_relation_Silica_Edu_solutions_d%inm.png' %(d*1e3))
#   
#list_mEv = np.array(list_omega*1e12*1e3*hb)
#%%

if save_txt == 1:
    
    header1 = 'omega/omega_{WG}     k_parallel [1/micrones]'  + ', '  + title + ', ' + name_this_py
    
    tabla_TM1A = np.array([list_omega_omegaWG_TM1A,zeros_eqTM1A])
    tabla_TM1A = np.transpose(tabla_TM1A)
    
    tabla_TM1B = np.array([list_omega_omegaWG_TM1B,zeros_eqTM1B])
    tabla_TM1B = np.transpose(tabla_TM1B)
    
    
    tabla_TM2A = np.array([list_omega_omegaWG_TM2A,zeros_eqTM2A])
    tabla_TM2A = np.transpose(tabla_TM2A)
    
    
    tabla_TM2B = np.array([list_omega_omegaWG_TM2B,zeros_eqTM2B])
    tabla_TM2B = np.transpose(tabla_TM2B)


    tabla_TE1A = np.array([list_omega_omegaWG_TE1A,zeros_eqTE1A])
    tabla_TE1A = np.transpose(tabla_TE1A)
    
    
    tabla_TE1B = np.array([list_omega_omegaWG_TE1B,zeros_eqTE1B])
    tabla_TE1B = np.transpose(tabla_TE1B)   


    tabla_TE2A = np.array([list_omega_omegaWG_TE2A,zeros_eqTE2A])
    tabla_TE2A = np.transpose(tabla_TE2A)
    
    
    tabla_TE2B = np.array([list_omega_omegaWG_TE2B,zeros_eqTE2B])
    tabla_TE2B = np.transpose(tabla_TE2B)

    
    np.savetxt('sol_TM1A_d%inm.txt' %(d*1e3), tabla_TM1A, fmt='%1.11e', delimiter='\t', header = header1)
    np.savetxt('sol_TM1B_d%inm.txt' %(d*1e3), tabla_TM1B, fmt='%1.11e', delimiter='\t', header = header1 )
    
    np.savetxt('sol_TM2A_d%inm.txt' %(d*1e3), tabla_TM2A, fmt='%1.11e', delimiter='\t', header = header1 )
    np.savetxt('sol_TM2B_d%inm.txt' %(d*1e3), tabla_TM2B, fmt='%1.11e', delimiter='\t', header = header1 )

    np.savetxt('sol_TE1A_d%inm.txt' %(d*1e3), tabla_TE1A, fmt='%1.11e', delimiter='\t', header = header1)
    np.savetxt('sol_TE1B_d%inm.txt' %(d*1e3), tabla_TE1B, fmt='%1.11e', delimiter='\t', header = header1 )

    np.savetxt('sol_TE2A_d%inm.txt' %(d*1e3), tabla_TE2A, fmt='%1.11e', delimiter='\t', header = header1 )
    np.savetxt('sol_TE2B_d%inm.txt' %(d*1e3), tabla_TE2B, fmt='%1.11e', delimiter='\t', header = header1 )

#%%

