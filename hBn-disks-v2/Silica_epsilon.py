#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
"""

import numpy as np
import sys
import os 
from scipy.interpolate import interp1d

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')

try:
    sys.path.insert(1, path_basic)
    from constants import constantes
except ModuleNotFoundError:
    print('constants.py no se encuentra en ' + path_basic)

pi,hb,c,alfac,mu1,mu2 = constantes()
aux = c*hb

graficar = 0


#%%

def epsilon_Silica(hbw): # min energy :  0.08684904004367106 max energy :  0.8064559363534132
    os.chdir(path_basic) 
    
    tabla_n = np.loadtxt('data_epsilon_Silica_n.txt',delimiter='\t', skiprows=1)
    tabla_n = np.transpose(tabla_n)
    [lambdda_micros, list_n] = tabla_n
    
    tabla_k = np.loadtxt('data_epsilon_Silica_k.txt',delimiter='\t', skiprows=1)
    tabla_k = np.transpose(tabla_k)
    [lambdda_micros, list_k] = tabla_k
    
    
    energy_eV = (hb*c)*2*np.pi/np.array(lambdda_micros)
    real_epsilon = np.array(list_n)**2 - np.array(list_k)**2
    imag_epsilon = 2*np.array(list_n)*np.array(list_k)
    
    try:
        f_real = interp1d(energy_eV, real_epsilon)
        f_imag = interp1d(energy_eV, imag_epsilon)
    
        return f_real(hbw) + 1j*f_imag(hbw)
    
    except ValueError: 
        print('El rango valido de energia es:', (np.min(energy_eV),np.max(energy_eV)),'eV',energy_eV )



def epsilon_Silica_extra(hbw): # min energy :  0.08684904004367106 max energy :  0.8064559363534132
    os.chdir(path_basic) 
    
    tabla_n = np.loadtxt('data_epsilon_Silica_n.txt',delimiter='\t', skiprows=1)
    tabla_n = np.transpose(tabla_n)
    [lambdda_micros, list_n] = tabla_n
    
    tabla_k = np.loadtxt('data_epsilon_Silica_k.txt',delimiter='\t', skiprows=1)
    tabla_k = np.transpose(tabla_k)
    [lambdda_micros, list_k] = tabla_k
    
    
    energy_eV = (hb*c)*2*np.pi/np.array(lambdda_micros)
    
    try:
        f_n = interp1d(energy_eV, list_n)
        f_k = interp1d(energy_eV, list_k)
    
        return f_n(hbw),f_k(hbw)
    
    except ValueError: 
        print('El rango valido de energia es:', (np.min(energy_eV),np.max(energy_eV)),'eV',energy_eV )
    
    

    

#%%    


if graficar ==1:
    
    path_save = path_basic + '/' + 'epsilon_silica'
    import matplotlib.pyplot as plt
    
    tamfig = [2.5, 2]
    tamletra = 7
    tamtitle  = 8
    tamnum = 6
    tamlegend = 6
    labelpady = 2
    labelpadx = 3
    pad = 2.5
    mk = 1
    ms = 1.5
    hp = 0.3
    length_marker = 1.5
    dpi = 500
    
    print('Graficar ε de Silica vs E')

    N = 1e3
    list_Energy_ev = np.linspace(0.087, 0.8, N)
    epsilon_real = np.real(epsilon_Silica(list_Energy_ev))
    epsilon_imag = np.imag(epsilon_Silica(list_Energy_ev))

    n_list,k_list = epsilon_Silica_extra(list_Energy_ev)


    plt.figure(figsize=tamfig)
    plt.plot(list_Energy_ev,epsilon_real,'.-',color = 'purple',ms=ms)
#    plt.title('Real part of ε from Silica',fontsize=tamtitle)
    plt.xlabel(r'$\hbar\omega$ (eV)',fontsize=tamletra,labelpad =labelpadx)
    plt.ylabel(r'Re(ε) of SiO2',fontsize=tamletra,labelpad =labelpady)
    plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in", pad = pad)
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 're_epsilon_silica' + '.png',bbox_inches='tight',pad_inches = 0.01, format='png', dpi=dpi)     
    
    plt.figure(figsize=tamfig)
    plt.plot(list_Energy_ev,epsilon_imag,'.-',color = 'purple',ms=ms)
#    plt.title('Imaginary part of ε from Silica',fontsize=tamtitle)
    plt.xlabel(r'$\hbar\omega$ (eV)',fontsize=tamletra,labelpad =labelpadx)
    plt.ylabel(r'Im(ε) of SiO2',fontsize=tamletra,labelpad =labelpady)
    plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in", pad = pad)
    plt.tight_layout()
    plt.savefig( 'im_epsilon_silica' + '.png',bbox_inches='tight',pad_inches = 0.01, format='png', dpi=dpi)     
    
    
    
    
    plt.figure(figsize=tamfig)
    plt.plot(list_Energy_ev,n_list,'.-',color = 'purple',ms=ms,label = 'n')
    plt.plot(list_Energy_ev,k_list,'.-',color = 'orange',ms=ms,label = 'k')
#    plt.title('Real part of ε from Silica',fontsize=tamtitle)
    plt.xlabel(r'$\hbar\omega$ (eV)',fontsize=tamletra,labelpad =labelpadx)
    plt.ylabel(r'$n$ and $k$ of SiO2',fontsize=tamletra,labelpad =labelpady)
    plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0.05,handletextpad=0.2, handlelength=length_marker)
    plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in", pad = pad)
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'n_k_silica' + '.png',bbox_inches='tight',pad_inches = 0.01, format='png', dpi=dpi)     
    
    
    

