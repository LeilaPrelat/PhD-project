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
#import seaborn as sns
import numpy as np


#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
#print('Importar modulos necesarios para este codigo')

try:
    sys.path.insert(1, path_basic)
    from constants import constantes
except ModuleNotFoundError:
    print('constants_plane.py no se encuentra en ' + path_basic)

pi,hb,c,alfac,mu1,mu2 = constantes()
aux = c*hb

#%%
    

def sigma_DL(hbw,epsi1,epsi3,d_nano):
    """
    Parameters
    ----------
    hbw : energia = hbar*omega en eV
    hbmu : potencial quimico del grafeno mu (en eV)
    hbgama : collision frequency in eV 

    Returns
    -------
    devuelve la conductividad de Silver
    del modelo drude lorentz (parte intra aproximada)
    en unidades de e^2/hbar = alfac*c (unidades de c en microns/seg)
    """

    d_micros = d_nano*1e-3
    
#    c = 3*10**(14)                ### light velocity in micron/seg
#    alfac = 1/137.0359            ### fine structure
    

    hb_omega_p = 9.17 #eV
    
    hb_gamma = 21*1e-3
    
    omega_p = hb_omega_p/hb
    
    hb_omega_D = hb_omega_p*omega_p*d_micros/(alfac*c)
    
    num1 = 1j*alfac*c*hb_omega_D


    num2 = hbw + 1j*hb_gamma
    
    return num1/num2

#%%


def k_parallel(hbw,epsi1,epsi3,d_nano): 
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

    omegac = hbw/(hb*c)
    num = (epsi1 + epsi3)*1j*omegac
    
    den  = 4*np.pi*sigma_DL(hbw,epsi1,epsi3,d_nano)/c
    

    return num/den


def k_parallel_v2(hbw,epsi1,epsi3,d_nano): 
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

    hb_omega_p = 9.17 #eV

    d_micros = d_nano*1e-3
    
    A = d_micros*(hb_omega_p**2)/(epsi1 + epsi3)
    
    return hbw**2/A

def epsilon_m(hbw):
    
    epsi_b = 4
    hbgamma = 21*1e-3
    
    
    num = (9.17)**2
    
    den = hbw*(hbw + 1j*hbgamma ) 
    
    return epsi_b - num/den


def k_parallel_Silver_guess(hbw,epsi1,epsi3,d_nano): 
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

    omegac = hbw/(hb*c)
    d_micros = d_nano*1e-3


    epsilon2 =  epsilon_m(hbw)  
    
    r12 = (epsilon2*np.sqrt(epsi1) - epsi1*np.sqrt(epsilon2))/(epsilon2*np.sqrt(epsi1) + epsi1*np.sqrt(epsilon2))
    
    r23 = (epsi3*np.sqrt(epsilon2) - epsilon2*np.sqrt(epsi3))/(epsi3*np.sqrt(epsilon2) + epsilon2*np.sqrt(epsi3))
    
    
    A = -np.real(np.log(r12*r23)*4/(1j*d_micros*np.sqrt(epsilon2)))
    
#    k_parallel = omegac/np.sqrt(A)
#    
    
    k_parallel = np.sqrt(4/((d_micros*epsilon2)**2) + omegac**2)
    
    return k_parallel


#%%
    
path_save = path_basic + '/' + 'disp_relation_Silver'

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

 
epsilon1, epsilon3 = 1,1
d_nano = 10 # thickness en nanometros
list_int_v1 = [1,5,10,15,20]
list_int_v2 = [1,25,50,75,100]
list_int_v3 = [0.1,0.2,0.3,1]
list_int_v4 = [0.07,0.1,0.2,1]

if d_nano == 100:
    list_int_v = list_int_v1
elif d_nano == 20:
    list_int_v = list_int_v2

elif d_nano == 10:
    list_int_v = list_int_v3

elif d_nano == 1:
    list_int_v = list_int_v4
    

list_colours = ['gold', 'orange', 'darksalmon','black', 'brown', 'darkred']

title1 = '$\epsilon_1 = %i$, $\epsilon_3 = %i$, d = %i nm, $\hbar\omega_{bulk}$ = 9.17 eV, $\hbar\gamma$ = 21 meV' %(epsilon1,epsilon3,d_nano)

list_E_ev = np.linspace(0.001,5.5,100)


list_y_re = []
list_y_im = []

for hb_omega in list_E_ev:
   
    
    valuey = k_parallel_v2(hb_omega,epsilon1,epsilon3,d_nano)

#    valuey_aprox = k_parallel_Silver_guess(hb_omega,epsi1,epsilon3,d_nano)
    
    list_y_re.append(valuey.real)
    list_y_im.append(valuey.imag)
    

plt.figure(figsize=tamfig)    
#plt.title(title1 ,fontsize=tamtitle)
plt.plot(list_y_re,list_E_ev,'-',color = 'darkgreen',lw = ms, label = 'plasmon')


j = 0
for int_v in list_int_v:
    if d_nano != 1:
        v = c*int_v
    else:
        v = c*int_v
    list_y2_re = []
    list_y2_im = []
    for hb_omega in list_E_ev:
        omega = hb_omega/hb        
        valuey_e = omega/v
        list_y2_re.append(valuey_e.real)
        list_y2_im.append(valuey_e.imag)  
    ##  
    if int_v != 1:
        plt.plot(list_y2_re,list_E_ev,'-', color = list_colours[j],lw = ms, label = 'v= %.2fc' %(int_v))
    else:
        plt.plot(list_y2_re,list_E_ev,'-', color = list_colours[j], lw = ms, label = 'light')
        
    j = j + 1
        
plt.xlabel('Parallel wave-vector $k_\parallel$ (1/$\mu$m)',fontsize=tamletra, labelpad = labelpadx)
plt.ylabel('Plasmon energy $\hbar\omega$ (eV)',fontsize=tamletra, labelpad = labelpady)
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=False,handletextpad=hp, handlelength=1)
plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in", pad = pad)
#plt.grid(1)
plt.tight_layout()
os.chdir(path_save)
plt.savefig('disp_relation_silver_real_vs_Ev_d%i.png' %(d_nano),dpi = dpi)

#%%

