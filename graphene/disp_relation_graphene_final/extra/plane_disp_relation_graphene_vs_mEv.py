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

graficar = 1

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
#path_graphene =  path_basic.replace('/' + 'plane','')
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
    
    f  = (hb*omega + 1j*gamma_in)/mu_c
    
    k0 = omega/c
    
    return cte*f*k0         

#%%

if graficar == 1:
    
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
    
    title1 = '$\epsilon_1 = %.1f$, $\epsilon_2 = %.1f$, $\gamma_{in}$ = %.4feV' %(epsilon1,epsilon2,gamma_in)
    title2 = '$\mu_c$ = %.2feV' %(mu_c)

    def k_parallel_ev(meV):
        omega = meV*1e-3/hb 
        return k_parallel(omega,mu_c,gamma_in,epsilon1,epsilon2)
    
    list_omega = np.linspace(0.001,4,n)
    list_y_re = []
    list_y_im = []
    for omega in list_omega:
        omegaTHz = omega*1e12
        valuey = k_parallel_ev(omegaTHz)
        list_y_re.append(valuey.real)
        list_y_im.append(valuey.imag)
        

    int_v = 5
    v = c/int_v
    
    def electron(v,omega):
        return omega*1e12/v
    
    def resta_function(omega):
        elec = electron(v,omega)
        sp = k_parallel_ev(omega)
        return np.abs(elec - sp)

    res = minimize(resta_function, 0.5*1e13, method='Nelder-Mead', tol=1e-13, 
                   options={'maxiter':1150})
    
    omega_c = res.x[0]*1e-12
    print(res)

    list_meV = np.linspace(0.001,80,n)
    list_y2_re = []
    list_y2_im = []
    for energy in list_meV:
        
        
        valuey_e = electron(v,omega)
        list_y2_re.append(valuey_e.real)
        list_y2_im.append(valuey_e.imag)  
##  

    
    plt.figure(figsize=tamfig)    
    plt.title(title1 + ', ' + title2 + '\n' + '$\omega_c$ = 0.7 10$^{12}$Hz'%(omega_c),fontsize=tamtitle)
    plt.plot(list_y_re,list_omega,'-',color = colors[0],ms = lw, label = 'SP')
#    plt.plot(list_y2_re,list_omega,'-',color = colors[1],ms = lw, label = 'v= c/%i' %(int_v))
    plt.xlabel('Re($k_\parallel$) [1/$\mu$m]',fontsize=tamletra, labelpad = 0)
    plt.ylabel('E [meV]',fontsize=tamletra, labelpad = 0)
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
    plt.ylabel('E [meV]',fontsize=tamletra, labelpad = 0)
    plt.legend(loc='best',markerscale=2,fontsize=int(tamlegend*0.8))
    plt.tick_params(labelsize = tamnum,pad = pad)
    plt.grid(1)
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig('disp_relation_graphene_imag.png')

#%%

