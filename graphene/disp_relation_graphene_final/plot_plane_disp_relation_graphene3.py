#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 19 20:56:58 2020

@author: leila

relacion de dispersion
solucion analitica
para un plano de Ag

"""
from scipy.optimize import fsolve 
import os 
import sys
import matplotlib.pyplot as plt
#import seaborn as sns
import numpy as np
#sns.set()

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_graphene =  path_basic.replace('/' + 'disp_relation_graphene_final','')
#print('Importar modulos necesarios para este codigo')


try:
    sys.path.insert(1, path_graphene)
    from graphene_sigma import sigma_DL
except ModuleNotFoundError:
    print('graphene_sigma.py no se encuentra en ' + path_graphene)


try:
    sys.path.insert(1, path_graphene)
    from constants import constantes
except ModuleNotFoundError:
    print('constants_plane.py no se encuentra en ' + path_graphene)

pi,hb,c,alfac,mu1,mu2 = constantes()
aux = c*hb

#%%
    
def k_parallel(omegaTHz,mu_c,gamma_in,epsilon1,epsilon2): 
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

#    
#    omegaTHz = omega*1e12   
    omega = omegaTHz*1e12
    
    cte = (epsilon1 + epsilon2)/alfac
    
    f  = (hb*omega + 1j*gamma_in)/mu_c  #mu_c en eV
    
    k0 = omega/c
    
    return cte*f*k0*0.25         



def k_parallel_num(hbw_mev,mu_c,gamma_in,epsi1,epsi2,k_parallel):
    omegac = hbw_mev*1e-3/(hb*c)
    u = k_parallel/omegac
    cond = 4*np.pi*alfac*sigma_DL(hbw_mev*1e-3,mu_c,gamma_in) #no es necesario dividir por c porque habria que multiplicar por c

    rp_den = epsi2*1j*u + epsi1*1j*u - cond*u**2
    
    return rp_den



#%%
n = 200
list_omega_THz = np.linspace(1,80,n)
list_mEv = np.array(list_omega_THz*1e12*1e3*hb)

epsi1,epsi2 = 1,1
gamma_in = 0.0001      # collision frequency in eV
mu_c = 0.3

def k_parallel_ev(omegaTHz):
    return k_parallel(omegaTHz,mu_c,gamma_in,epsi1,epsi2)


zeros_eqTM1B = []
list_energy_TM1B = []

k_TM1 = 0
for E in list_mEv:
    omega = E*1e-3/(hb)
    omegaTHz = omega*1e12
#    omegaTHz = omega*   
    
    
    if k_TM1 == 0:
        init_condit = [np.real(k_parallel(omegaTHz,mu_c,gamma_in,epsi1,epsi2))]
    else:
        init_condit = zeros_eqTM1B[k_TM1-1]
        
        
    def k_parallel_WG_TM1_1var(k_parallel):   
        return  np.abs(k_parallel_num(E,mu_c,gamma_in,epsi1,epsi2,k_parallel))
    

    resTM1 = fsolve(k_parallel_WG_TM1_1var, init_condit,maxfev = 1000 )         
    resTM1_v = np.float(resTM1)
    if resTM1_v < 0 :
        resTM1_v = -resTM1_v
    
    zeros_eqTM1B.append(resTM1_v)
    list_energy_TM1B.append(E)
    k_TM1 = k_TM1 + 1


#%%
    
path_save = path_basic + '/' + 'disp_relation_graphene'

tamfig = [2.5, 2]
tamletra = 7
tamtitle  = 8
tamnum = 6
tamlegend = 6
labelpady = 2
labelpadx = 3
pad = 2.5
mk = 1
ms = 1
hp = 0.5
length_marker = 1.5
dpi = 500

n  = 150

list_colours = ['gold', 'orange', 'darksalmon','black', 'brown', 'darkred']

title1 = '$\epsilon_1 = %.1f$, $\epsilon_2 = %.1f$, $\gamma_{in}$ = %.2fmeV' %(epsi1,epsi2,gamma_in*1e3)
title2 = '$\mu_c$ = %.2feV' %(mu_c)



list_y_re = []
list_y_im = []
for omegaTHz in list_omega_THz:
   
    valuey = k_parallel(omegaTHz,mu_c,gamma_in,epsi1,epsi2)

    list_y_re.append(valuey.real)
    list_y_im.append(valuey.imag)
    


plt.figure(figsize=tamfig)    
#plt.title(title1 + ', ' + title2 + '\n',fontsize=tamtitle)
plt.plot(list_y_re,list_mEv,'-',color = 'darkgreen',lw = ms, label = 'plasmon')



j = 0
for int_v in [0.1,0.15,0.3,1]:
    v = c*int_v
    
    list_y2_re = []
    list_y2_im = []
    for omegaTHz in list_omega_THz:
#        omegaTHz = omega*1e12        
        valuey_e = omegaTHz*1e12/v
        list_y2_re.append(valuey_e.real)
        list_y2_im.append(valuey_e.imag)  
    ##  
    if int_v != 1:
        plt.plot(list_y2_re,list_mEv,'-', color = list_colours[j],lw = ms, label = 'v= %.2fc' %(int_v))
    else:
        plt.plot(list_y2_re,list_mEv,'-', color = list_colours[j], lw = ms, label = 'light')
        
    j = j + 1
        
plt.xlabel('Parallel wave-vector $k_\parallel$ (1/$\mu$m)',fontsize=tamletra, labelpad = 2)
plt.ylabel('Plasmon energy $\hbar\omega$ (meV)',fontsize=tamletra, labelpad = 2)
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=False,handletextpad=hp, handlelength=1)
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in", pad = pad)
#    plt.title(title,fontsize=int(tamtitle*0.9))
#plt.grid(1)
plt.tight_layout()
# os.chdir(path_save)
plt.savefig('disp_relation_graphene_real_vs_mEv.png', dpi=dpi)

plt.figure(figsize=tamfig)      
plt.title(title1 + ', ' + title2,fontsize=tamtitle)
plt.plot(list_y_im,list_mEv,'-',ms = ms)
plt.xlabel('Im($k_\parallel$) [1/$\mu$m]',fontsize=tamletra, labelpad = 0)
plt.ylabel('$\hbar\omega$ [meV]',fontsize=tamletra, labelpad = 0)
plt.legend(loc='best',markerscale=2,fontsize=int(tamlegend*0.8))
plt.tick_params(labelsize = tamnum,pad = pad)
plt.grid(1)
plt.tight_layout()
# os.chdir(path_save)
plt.savefig('disp_relation_graphene_imag_vs_mEv.png')

#
#list_vc = np.linspace(0.01,1,n)
#list_omega = np.linspace(0.01,15,n)
#list_y_re = []
#list_y_im = []
#for j in range(n):
#    omega = list_omega[j]
#    v = c*list_vc[j]
#    omegaTHz = omega*1e12
#    k0 = omegaTHz/v
#
#    cte = (epsilon1 + epsilon2)/alfac
#    f  = (hb*omegaTHz + 1j*gamma_in)/mu_c  #mu_c en eV          
#    valuey = cte*f*k0*0.25 
#
#    list_y_re.append(valuey.real)
#    list_y_im.append(valuey.imag)
#    
#fig,ax = plt.subplots(figsize=tamfig)
#ax.set_title(title1 + ', ' + title2 + '\n',fontsize=tamtitle)
#ax.plot(list_vc,list_y_re,'-',color = 'darkgreen',ms = lw, label = 'DR')
#ax.set_ylabel('Re($k_\parallel$) [1/$\mu$m]',fontsize=tamletra, labelpad = 0)
#ax.set_xlabel("v/c",fontsize=tamletra, labelpad = 0)
#ax.tick_params(labelsize = tamnum,pad = 0)
#
#ax2 = ax.twinx()
#ax2.tick_params(labelsize = tamnum,pad = pad)
#ax2.set_ylabel('$\omega$ [THz]',color="blue",fontsize=14)
## plt.grid(1)
#fig.tight_layout()
## os.chdir(path_save)
#plt.savefig('disp_relation_graphene_real_vs_v.png')

#%%

