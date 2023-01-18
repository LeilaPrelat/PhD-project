#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
conductividad de Ag
"""
import numpy as np
graficar = 0

#%%

def sigma_DL(omega,omega_bulk,gamma_in,d,epsilon_b): 
    """
    Parameters
    ----------
    omega : frecuencia en Hz
    omega_bulk : frecuencia de bulk Hz
    gamma_in : frecuencia de colision en eV
    d : espesor del plano Ag en micrometros
    epsilon_b : permeabilidad electrica del plano 
    Returns
        conductividad del Ag 
    -------

    """
    hb = 6.58211899*10**(-16)     ### Planck constant hbar in eV*s
    den = omega*(omega + 1j*gamma_in/hb)
    num = omega_bulk**2
    f = 1 - epsilon_b + num/den 

    return f*d*1j*omega/(4*np.pi)

#%%

if graficar == 1:

    import matplotlib.pyplot as plt
    import sys
    import seaborn as sns
    import os
    sns.set()
    name_this_py = os.path.basename(__file__)
    path = os.path.abspath(__file__) #path absoluto del .py actual
    path_basic = path.replace('/' + name_this_py,'')
    path_save = path_basic + '/' + 'sigma_Ag'
    
    try:
        sys.path.insert(1, path_basic)
        from constants import constantes
    except ModuleNotFoundError:
        print('constants.py no se encuentra en ' + path_basic)
    
    pi,hb,c,alfac,mu1,mu2 = constantes()
    
    aux = hb*c
    tamfig = (4.5,3.5)
    tamlegend = 10
    tamletra = 10
    tamtitle = 10
    tamnum = 9
    loc2 = [0,1]
    pad = -2
    lw = 1.5

    colors = ['darkred','steelblue','coral','yellowgreen']
     
    gamma_in = 0.0001      # collision frequency in eV
    Ebulk = 9.17
    omega_bulk = Ebulk/hb
    d = 0.0001 #micrones, 10 nm
    epsilon_b = 4
    hb = 6.58211899*10**(-16)     ### Planck constant hbar in eV*s
    
    n  = 150
    list_E = np.linspace(0.2,1,n)
    list_omega = []
    list_y_re = []
    list_y_im = []
    for E in list_E:
        omega = E/hb
        valuey = sigma_DL(omega,omega_bulk,gamma_in,d,epsilon_b)
        list_omega.append(omega)
        list_y_re.append(valuey.real)
        list_y_im.append(valuey.imag)
            
    
    title1 = '$\gamma_{in}$ = %.4feV, $\epsilon_b$ = %i' %(gamma_in,epsilon_b)
    title2 = '$E_{bulk}$ = %.2feV, d = %.4f$\mu$m' %(Ebulk,d)
    
    plt.figure(figsize=tamfig)    
    plt.title(title1 + '\n' + title2,fontsize=tamtitle)
    plt.plot(list_E,list_y_re,'-',color = colors[0],ms = lw)
    plt.ylabel('Re($\sigma$) Ag',fontsize=tamletra, labelpad = 0)
    plt.xlabel('energy [eV]',fontsize=tamletra, labelpad = 0)
    plt.tick_params(labelsize = tamnum,pad = pad)
    plt.grid(1)
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig('sigma_Ag_real.png')

    plt.figure(figsize=tamfig)      
    plt.title(title1 + '\n' + title2,fontsize=tamtitle)
    plt.plot(list_E,list_y_im,'-',color = colors[0],ms = lw)
    plt.ylabel('Im($\sigma$) Ag',fontsize=tamletra, labelpad = 0)
    plt.xlabel('energy [eV]',fontsize=tamletra, labelpad = 0)
    plt.tick_params(labelsize = tamnum,pad = pad)
    plt.grid(1)
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig('sigma_Ag_imag.png')

#%%
