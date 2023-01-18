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
from scipy.optimize import minimize   

sns.set()

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
#path_ctes =  path_basic.replace('/' + 'Silica','')
#print('Importar modulos necesarios para este codigo')
path_save = path_basic + '/' + 'disp_relation_Silica'
try:
    sys.path.insert(1, path_basic)
    from constants import constantes
except ModuleNotFoundError:
    print('constants.py no se encuentra en ' + path_basic)

pi,hb,c,alfac,mu1,mu2,mu3 = constantes()
aux = c*hb

#%%
    
def k_parallel_WG_TM1(omega,epsilon2,d,k_parallel):
    k0 = omega/c
    chi_prima = d*np.sqrt(epsilon2*k0**2 - k_parallel**2)
    chi = d*np.sqrt(k_parallel**2 - k0**2)   
    return np.abs(np.tan(chi_prima/2) - epsilon2*chi/chi_prima)


def k_parallel_WG_TM2(omega,epsilon2,d,k_parallel):
    k0 = omega/c
    chi_prima = d*np.sqrt(epsilon2*k0**2 - k_parallel**2)
    chi = d*np.sqrt(k_parallel**2 - k0**2)
    return np.abs(1/np.tan(chi_prima/2) + epsilon2*chi/chi_prima)

def k_parallel_WG_TE1(omega,d,k_parallel):
    k0 = omega/c
    chi_prima = d*np.sqrt(epsilon2*k0**2 - k_parallel**2)
    chi = d*np.sqrt(k_parallel**2 - k0**2)
    return np.abs(np.tan(chi_prima/2) - chi/chi_prima)


def k_parallel_WG_TE2(omega,d,k_parallel):
    k0 = omega/c
    chi_prima = d*np.sqrt(epsilon2*k0**2 - k_parallel**2)
    chi = d*np.sqrt(k_parallel**2 - k0**2)
    return np.abs(1/np.tan(chi_prima/2) + chi/chi_prima)


def k_parallel_air(omega):

    return omega/c

def k_parallel_medium(omega,epsilon2):

    return omega*np.sqrt(epsilon2)/c

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

#%%


 
epsilon2 = 12
d = 10*1e-3

hb = 6.582118989999999e-16   #electron volts/seg
alfac = 1/137
c = 299792458000000 #microns/seg

list_colours = ['black', 'gold', 'orange', 'darksalmon', 'brown', 'darkred']

title = '$\epsilon_2 = %.1f$, d = %inm' %(epsilon2,d*1e3)

list_omega_omegaWG = np.linspace(0.01,4,150)

list_y_WG_TM1A = []
list_y_WG_TE1A = []

list_y_WG_TM2A = []
list_y_WG_TE2A = []

list_y_WG_TM1B = []
list_y_WG_TE1B = []

list_y_WG_TM2B = []
list_y_WG_TE2B = []



list_y_air = []
list_y_medium = []
list_x = []

omegaWG = (np.pi*c/d)/np.sqrt(epsilon2-1) 
cte_new = c/omegaWG
    
for omega_omegaWG in list_omega_omegaWG:
    
    omega = omega_omegaWG*omegaWG
#    list_x.append(omega/omegaWG)


    def k_parallel_WG_TM_minimiza1(k_parallel):
        return k_parallel_WG_TM1(omega,epsilon2,d,k_parallel)

    def k_parallel_WG_TE_minimiza1(k_parallel):
        return k_parallel_WG_TE1(omega,d,k_parallel)


    def k_parallel_WG_TM_minimiza2(k_parallel):
        return k_parallel_WG_TM2(omega,epsilon2,d,k_parallel)

    def k_parallel_WG_TE_minimiza2(k_parallel):
        return k_parallel_WG_TE2(omega,d,k_parallel)
    
    init_cond = omega_omegaWG/cte_new

    k_air = k_parallel_air(omega)
    k_medium = k_parallel_medium(omega,epsilon2)
   
    res1A = minimize(k_parallel_WG_TM_minimiza1, k_air, method='Nelder-Mead', tol=1e-9, 
                   options={'maxiter':1500})
    
    res1B = minimize(k_parallel_WG_TM_minimiza1, k_medium, method='Nelder-Mead', tol=1e-9, 
                   options={'maxiter':1500})
        
    
    res2A = minimize(k_parallel_WG_TE_minimiza1, k_air, method='Nelder-Mead', tol=1e-9, 
                   options={'maxiter':1500})

    res2B = minimize(k_parallel_WG_TE_minimiza1, k_medium, method='Nelder-Mead', tol=1e-9, 
                   options={'maxiter':1500})    

    res3A = minimize(k_parallel_WG_TM_minimiza2, k_air, method='Nelder-Mead', tol=1e-9, 
                   options={'maxiter':1500})


    res3B = minimize(k_parallel_WG_TM_minimiza2, k_medium, method='Nelder-Mead', tol=1e-9, 
                   options={'maxiter':1500})
    
    res4A = minimize(k_parallel_WG_TE_minimiza2, k_air, method='Nelder-Mead', tol=1e-9, 
                   options={'maxiter':1500})
    
    res4B = minimize(k_parallel_WG_TE_minimiza2, k_medium, method='Nelder-Mead', tol=1e-9, 
                   options={'maxiter':1500})


    list_y_WG_TM1A.append(res1A.x[0]*cte_new)
    list_y_WG_TE1A.append(res2A.x[0]*cte_new)
    
    list_y_WG_TM2A.append(res3A.x[0]*cte_new)
    list_y_WG_TE2A.append(res4A.x[0]*cte_new)
    
    list_y_WG_TM1B.append(res1B.x[0]*cte_new)
    list_y_WG_TE1B.append(res2B.x[0]*cte_new)
    
    list_y_WG_TM2B.append(res3B.x[0]*cte_new)
    list_y_WG_TE2B.append(res4B.x[0]*cte_new)    

    
    list_y_air.append(k_air*cte_new)
    list_y_medium.append(k_medium*cte_new)
    
    
#list_mEv = np.array(list_omega*1e12*1e3*hb)
#%%
plt.figure(figsize=tamfig)    
plt.title(title,fontsize=tamtitle)
plt.plot(list_y_WG_TE1A,list_omega_omegaWG,'.-',color = 'red',ms = 3, label = 'WG TE')
plt.plot(list_y_WG_TE1B,list_omega_omegaWG,'.-',color = 'red',ms = 3)

plt.plot(list_y_WG_TM1A,list_omega_omegaWG,'.-',color = 'blue',ms = 3, label = 'WG TM')
plt.plot(list_y_WG_TM1B,list_omega_omegaWG,'.-',color = 'blue',ms = 3)

plt.plot(list_y_WG_TE2A,list_omega_omegaWG,'.-',color = 'red',ms = 3)
plt.plot(list_y_WG_TE2B,list_omega_omegaWG,'.-',color = 'red',ms = 3)

plt.plot(list_y_WG_TM2A,list_omega_omegaWG,'.-',color = 'blue',ms = 3)
plt.plot(list_y_WG_TM2B,list_omega_omegaWG,'.-',color = 'blue',ms = 3)

plt.plot(list_y_air,list_omega_omegaWG,'--',color = 'grey',ms = lw, label = 'air')
plt.plot(list_y_medium,list_omega_omegaWG,'--',color = 'black',ms = lw, label = 'medium')
plt.xlabel(r'$k_\parallel$c/$\omega_{WG}$',fontsize=tamletra, labelpad = 0)
plt.ylabel(r'$\omega/\omega_{WG}$',fontsize=tamletra, labelpad = 0)
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
plt.tick_params(labelsize = tamnum,pad = pad)
plt.grid(1)
plt.tight_layout()
# os.chdir(path_save)
plt.savefig('disp_relation_Silica_real_vs_mEv.png')

#%%

