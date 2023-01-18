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

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_graphene =  path_basic.replace('/' + 'plane','')
#print('Importar modulos necesarios para este codigo')

try:
    sys.path.insert(1, path_basic)
    from constants import constantes
except ModuleNotFoundError:
    print('constants_plane.py no se encuentra en ' + path_basic)

pi,hb,c,alfac,mu1,mu2 = constantes()
aux = c*hb

#%%
    
def omega_function(k_parallel,mu_c,gamma_in,epsilon1,epsilon2): 
    """
    Parameters
    ----------
    omega : frecuencia en Hz
    mu_c : chemical potential of graphene in eV
    gamma_in : frecuencia de colision en eV
    epsilon1 : permeabilidad electrica del medio de arriba
    epsilon2 : permeabilidad electrica del medio de abajo
    Returns
        relacion de dispersion para un plano de Ag
        (solucion analitica)
    -------
    """
    cte = (epsilon1 + epsilon2)/alfac
    
    f  = (hb*omega + 1j*gamma_in)/mu_c
    
    k0 = omega/c
    
    return cte*f*k0*0.25         

#%%
    
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.optimize import minimize   

sns.set()
path_save = path_basic + '/' + 'disp_relation_graphene'
colors = ['darkred','steelblue','coral','yellowgreen']

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

 
epsilon1,epsilon2 = 1,1

gamma_in = 0.0001      # collision frequency in eV
mu_c = 0.3
n  = 150

title1 = '$\epsilon_1 = %.1f$, $\epsilon_2 = %.1f$, $\gamma_{in}$ = %.2fmeV' %(epsilon1,epsilon2,gamma_in*1e3)
title2 = '$\mu_c$ = %.2feV' %(mu_c)

def k_parallel_ev(omega):
    return k_parallel(omega,mu_c,gamma_in,epsilon1,epsilon2)

list_Q = np.linspace(0.05,2.5,n)
list_y_re = []
list_y_im = []
for omega in list_omega:
    omegaTHz = omega*1e12
    valuey = k_parallel_ev(omegaTHz)
    list_y_re.append(valuey.real)
    list_y_im.append(valuey.imag)
    

plt.figure(figsize=tamfig)    
plt.title(title1 + ', ' + title2 + '\n',fontsize=tamtitle)
plt.plot(list_y_re,list_omega,'-',color = colors[0],ms = lw, label = 'SP')

for int_v in [1,0.1,0.2,0.3,0.4]:
    v = c*int_v
    
    def electron(v,omega):
        return omega/v
    
    list_y2_re = []
    list_y2_im = []
    for omega in list_omega:
        valuey_e = electron(v,omega)
        list_y2_re.append(valuey_e.real)
        list_y2_im.append(valuey_e.imag)  
    ##  

    plt.plot(list_y2_re,list_omega,'-',ms = lw, label = 'v= c/%i' %(int_v))
plt.xlabel('Re($k_\parallel$) [1/$\mu$m]',fontsize=tamletra, labelpad = 0)
plt.ylabel('$\omega$ [THz]',fontsize=tamletra, labelpad = 0)
plt.legend(loc='best',markerscale=2,fontsize=int(tamlegend*0.8))
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
plt.tick_params(labelsize = tamnum,pad = pad)
plt.grid(1)
plt.tight_layout()
os.chdir(path_save)
plt.savefig('disp_relation_graphene_real.png')

plt.figure(figsize=tamfig)      
plt.title(title1 + ', ' + title2,fontsize=tamtitle)
plt.plot(list_y_im,list_omega,'-',color = colors[0],ms = lw)
plt.xlabel('Im($k_\parallel$) [1/$\mu$m]',fontsize=tamletra, labelpad = 0)
plt.ylabel('$\omega$ [THz]',fontsize=tamletra, labelpad = 0)
plt.legend(loc='best',markerscale=2,fontsize=int(tamlegend*0.8))
plt.tick_params(labelsize = tamnum,pad = pad)
plt.grid(1)
plt.tight_layout()
os.chdir(path_save)
plt.savefig('disp_relation_graphene_imag.png')

#%%

