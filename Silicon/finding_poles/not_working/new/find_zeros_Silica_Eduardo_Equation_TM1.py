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
from scipy.optimize import minimize,fsolve 
from mpmath import findroot

sns.set()

pedir_2_condiciones_rel_disp = 0

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


def k_parallel_WG_TM1(omega,epsilon2,d,k_parallel):
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

#%%
    
tamfig = (4.5,3.5)
tamlegend = 10
tamletra = 10
tamtitle = 10
tamnum = 9
loc2 = [0,1]
pad = -2
lw = 1.5
hp = 0.3
mk = 2
labelpady = 0
labelpadx = 0

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
N = 200
list_omega_omegaWGA = np.linspace(cota_inf_omega_omegaWG,cota_sup_omega_omegaWG,N)

#list_omega_omegaWG = np.linspace(4,0.005,N)
omegaWG = (np.pi*c/d)/(np.sqrt(epsilon2-1)) 
cte_new = c/omegaWG
        
list_k_parallel = np.linspace(0.001,10,N)

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
    
    print(k_air)
    
    list_y_air.append(k_air*cte_new)
    list_y_medium.append(k_medium*cte_new)

    init_condit = [k_air]
            
    def k_parallel_WG_TM1_1var(k_parallel):   
        return  np.abs(k_parallel_WG_TM1(omega,epsilon2,d,k_parallel))
    
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
        
#%%
        

#%%
"""
zeros_eqTM1B = []
init_condit_TM1B = []
list_omega_omegaWG_TM1B = []

k_TM1 = 0
for omega_omegaWG in list_omega_omegaWGB:
    omega = omega_omegaWG*omegaWG
    k_air = k_parallel_air(omega)
    k_medium = k_parallel_medium(omega,epsilon2)   
        
    init_condit = [k_air]
            
    def k_parallel_WG_TM1_1var(k_parallel):   
        return  np.abs(k_parallel_WG_TM1(omega,epsilon2,d,k_parallel))
    
    if k_TM1 == 0 :
        resTM1 = fsolve(k_parallel_WG_TM1_1var, [k_air])
    else:
        resTM1 = fsolve(k_parallel_WG_TM1_1var, [init_condit_TM1A[k_TM1-1]] )         
    resTM1_v = np.float(resTM1)
    if k_air < resTM1_v < k_medium:
        zeros_eqTM1B.append(resTM1_v)
        list_omega_omegaWG_TM1B.append(omega_omegaWG)
        k_TM1 = k_TM1 + 1
    #            if resTM1_v*cte_new < 7:
    #                init_condit_TM1A.append(resTM1_v)
    #            else:
    #                init_condit_TM1A.append(k_medium)
        
        if omega_omegaWG > 2.5:
            init_condit_TM1B.append(resTM1_v)
        else:
            init_condit_TM1B.append(k_medium)

#%%


def k_parallel_WG_TM1_2var(k_parallel_c_omegaWG,omega_omegaWG):   
    omega = omega_omegaWG*omegaWG
    k_parallel = k_parallel_c_omegaWG/cte_new
    return  np.abs(k_parallel_WG_TM1(omega,epsilon2,d,k_parallel))

#number_of_sol = 2*int(omega_omegaWG)
list_x = list_k_parallel
list_y = list_omega_omegaWGA
f = np.vectorize(k_parallel_WG_TM1_2var)
X, Y = np.meshgrid(list_x, list_y)
Z = f(X, Y)
maxi = np.max(Z)
Z = Z/maxi
limits = [min(list_x) , max(list_x), min(list_y) , max(list_y)]


labelx = r'$k_\parallel$c/$\omega_{WG}$'
labely = r'$\omega/\omega_{WG}$'
labelz = 'TM1'

graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
im = plt.imshow(Z, extent = limits, cmap=plt.cm.hot, interpolation='bilinear')
cbar = plt.colorbar(im)
cbar.ax.tick_params(labelsize = tamnum)
cbar.set_label(labelz,fontsize=tamlegend,labelpad = 1)
plt.plot(np.array(zeros_eqTM1A)*cte_new,list_omega_omegaWG_TM1A,'.',color = 'lightblue',ms = 4)
plt.plot(np.array(zeros_eqTM1B)*cte_new,list_omega_omegaWG_TM1B,'.',color = 'blue',ms = 4)
plt.tick_params(labelsize = tamnum,pad = pad)

"""

#%%
   
#
plt.figure(figsize=tamfig)    
plt.title(title,fontsize=tamtitle)
plt.plot(np.array(zeros_eqTM1A)*cte_new,list_omega_omegaWG_TM1A,'.',color = 'lightblue',ms = 4)
plt.xlabel(r'$k_\parallel$c/$\omega_{WG}$',fontsize=tamletra, labelpad = 0)
plt.ylabel(r'$\omega/\omega_{WG}$',fontsize=tamletra, labelpad = 0)
plt.plot(list_y_air,list_omega_omegaWGA,'--',color = 'grey',ms = lw, label = 'air')
plt.plot(list_y_medium,list_omega_omegaWGA,'--',color = 'black',ms = lw, label = 'medium')
plt.xlim([-0.5,8.5])
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
plt.grid(1)
plt.tight_layout()
plt.savefig('disp_relation_Silicon_TM1_d%inm.png' %(d*1e3))

## os.chdir(path_save)

#   
#list_mEv = np.array(list_omega*1e12*1e3*hb)
#%%
