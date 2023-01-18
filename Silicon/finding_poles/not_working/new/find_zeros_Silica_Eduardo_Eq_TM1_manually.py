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
    
    return np.tan(chi_prima/2) - epsilon2*chi/chi_prima

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

omegaWG = (np.pi*c/d)/(np.sqrt(epsilon2-1)) 
cte_new = c/omegaWG


list_y_air = []
list_y_medium = []


for omega_omegaWG in list_omega_omegaWGA:
    omega = omega_omegaWG*omegaWG
    k_air = k_parallel_air(omega)
    k_medium = k_parallel_medium(omega,epsilon2)   
    
    list_y_air.append(k_air*cte_new)
    list_y_medium.append(k_medium*cte_new)

#%%

list_k_parallel = np.linspace(0.0001,15,4000)

zeros_eqTM1A = []
init_condit_TM1A = []
list_omega_omegaWG_TM1A = []

            
for omega_omegaWG in list_omega_omegaWGA:
    omega = omega_omegaWG*omegaWG    
    ind = 0
    list_num_imag = []
    list_num_real = []
    for k_parallel in list_k_parallel:
    
        num = np.abs(k_parallel_WG_TM1(omega,epsilon2,d,k_parallel))
    
#    
#        if ind > 0:
#    if list_num_imag[ind]*list_num_imag[ind-1] < 0 :
    
        if num < 1 : 
            
            print('sol')
            
            zeros_eqTM1A.append(k_parallel)
            list_omega_omegaWG_TM1A.append(omega_omegaWG)
                
#        ind = ind + 1

#%%
   
#
plt.figure(figsize=tamfig)    
plt.title(title,fontsize=tamtitle)
plt.plot(np.array(zeros_eqTM1A)*cte_new,list_omega_omegaWG_TM1A,'.',color = 'lightblue',ms = 4)
plt.tick_params(labelsize = tamnum,pad = pad)
plt.xlabel(r'$k_\parallel$c/$\omega_{WG}$',fontsize=tamletra, labelpad = 0)
plt.ylabel(r'$\omega/\omega_{WG}$',fontsize=tamletra, labelpad = 0)
plt.plot(list_y_air,list_omega_omegaWGA,'--',color = 'grey',ms = lw, label = 'air')
plt.plot(list_y_medium,list_omega_omegaWGA,'--',color = 'black',ms = lw, label = 'medium')
plt.xlim([-0.5,8.5])
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
plt.grid(1)
plt.tight_layout()

# os.chdir(path_save)
#plt.savefig('disp_relation_Silica_Edu_solutions_d%inm.png' %(d*1e3))
#   

#%%
