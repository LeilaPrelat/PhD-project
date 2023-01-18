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
    
def k_parallel(omega,mu_c,gamma_in,epsilon1,epsilon2): 
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
    
    f  = (hb*omega + 1j*gamma_in)/mu_c  #mu_c en eV
    
    k0 = omega/c
    
    return cte*f*k0*0.25         

def k_parallel_ev(omega):
    return k_parallel(omega,mu_c,gamma_in,epsilon1,epsilon2)

#%%
    
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.optimize import minimize   

sns.set()
path_save = path_basic + '/' + 'disp_relation_graphene'

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

hb = 6.582118989999999e-16   #electron volts/seg
alfac = 1/137
c = 299792458000000 #microns/seg

list_colours = ['black', 'gold', 'orange', 'darksalmon', 'brown', 'darkred']

title1 = '$\epsilon_1 = %.1f$, $\epsilon_2 = %.1f$, $\gamma_{in}$ = %.4feV' %(epsilon1,epsilon2,gamma_in)
title2 = '$\mu_c$ = %.2feV' %(mu_c)

list_omega = np.linspace(1,80,n)
list_y_re = []
list_y_im = []
for omega in list_omega:
    omegaTHz = omega*1e12

    cte = (epsilon1 + epsilon2)/alfac
    f  = (hb*omegaTHz + 1j*gamma_in)/mu_c  #mu_c en eV    
    k0 = omegaTHz/c        
    valuey = cte*f*k0*0.25 

    list_y_re.append(valuey.real)
    list_y_im.append(valuey.imag)
    

plt.figure(figsize=tamfig)    
plt.title(title1 + ', ' + title2 + '\n',fontsize=tamtitle)
plt.plot(list_y_re,list_omega,'-',color = 'darkgreen',ms = lw, label = 'DR')

j = 0
for int_v in [1,0.1,0.2,0.3,0.4,0.5]:
    v = c*int_v
    
    list_y2_re = []
    list_y2_im = []
    for omega in list_omega:
        omegaTHz = omega*1e12        
        valuey_e = omegaTHz/v
        list_y2_re.append(valuey_e.real)
        list_y2_im.append(valuey_e.imag)  
    ##  
    if int_v != 1:
        plt.plot(list_y2_re,list_omega,'-', color = list_colours[j],ms = lw, label = 'v= %.1fc' %(int_v))
    else:
        plt.plot(list_y2_re,list_omega,'-', color = list_colours[j], ms = lw, label = 'light')
        
    j = j + 1
        
plt.xlabel('$k_\parallel$ [1/$\mu$m]',fontsize=tamletra, labelpad = 0)
plt.ylabel('$\omega$ [THz]',fontsize=tamletra, labelpad = 0)
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
plt.tick_params(labelsize = tamnum,pad = pad)
plt.grid(1)
plt.tight_layout()
# os.chdir(path_save)
plt.savefig('disp_relation_graphene_real_vs_Omega.png')

plt.figure(figsize=tamfig)      
plt.title(title1 + ', ' + title2,fontsize=tamtitle)
plt.plot(list_y_im,list_omega,'-',ms = lw)
plt.xlabel('Im($k_\parallel$) [1/$\mu$m]',fontsize=tamletra, labelpad = 0)
plt.ylabel('$\omega$ [THz]',fontsize=tamletra, labelpad = 0)
plt.legend(loc='best',markerscale=2,fontsize=int(tamlegend*0.8))
plt.tick_params(labelsize = tamnum,pad = pad)
plt.grid(1)
plt.tight_layout()
# os.chdir(path_save)
plt.savefig('disp_relation_graphene_imag_vs_Omega.png')


list_vc = np.linspace(0.01,1,n)
list_omega = np.linspace(0.01,15,n)
list_y_re = []
list_y_im = []
for j in range(n):
    omega = list_omega[j]
    v = c*list_vc[j]
    omegaTHz = omega*1e12
    k0 = omegaTHz/v

    cte = (epsilon1 + epsilon2)/alfac
    f  = (hb*omegaTHz + 1j*gamma_in)/mu_c  #mu_c en eV          
    valuey = cte*f*k0*0.25 

    list_y_re.append(valuey.real)
    list_y_im.append(valuey.imag)
    
fig,ax = plt.subplots(figsize=tamfig)
ax.set_title(title1 + ', ' + title2 + '\n',fontsize=tamtitle)
ax.plot(list_vc,list_y_re,'-',color = 'darkgreen',ms = lw, label = 'DR')
ax.set_ylabel('Re($k_\parallel$) [1/$\mu$m]',fontsize=tamletra, labelpad = 0)
ax.set_xlabel("v/c",fontsize=tamletra, labelpad = 0)
ax.tick_params(labelsize = tamnum,pad = 0)

ax2 = ax.twinx()
ax2.tick_params(labelsize = tamnum,pad = pad)
ax2.set_ylabel('$\omega$ [THz]',color="blue",fontsize=14)
# plt.grid(1)
fig.tight_layout()
# os.chdir(path_save)
plt.savefig('disp_relation_graphene_real_vs_v.png')

#%%

