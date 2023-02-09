#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: leila

graficar polarizabilidad
"""

import numpy as np
import sys
import os 
import matplotlib.pyplot as plt
import seaborn as sns

sns.set()

#%%

save_graphs = 1 #guardar los graficos 2D del campo
    
#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/polarizability' ,'')
path_save = path_basic + '/' + 'polarizability'

if save_graphs==1:    
    if not os.path.exists(path_save):
        print('Creating folder to save graphs')
        os.mkdir(path_save)

err = 'function_polarizability.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from function_polarizability import alpha_function
except ModuleNotFoundError:
    print(err)

try:
    sys.path.insert(1, path_basic)
    from function_polarizability import alpha_functionMIE
except ModuleNotFoundError:
    print(err)

try:
    sys.path.insert(1, path_constants)
    from constants import constantes
except ModuleNotFoundError:
    print('constants.py no se encuentra en ' + path_constants)

pi,hb,c,alfac,mu1,mu2 = constantes()
aux = c*hb

#%%

print('Definir parametros del problema')

#omega = 0.7*1e12

epsi1 = 1

omega0THz = 0.05
omega0 = omega0THz*1e12 
R = 10 # 10 nm en unidades de micrometros
kappa_factor_omega0 = 0.1
kappa_r_factor= 0.5

list_OmegaTHz = np.linspace(0.001,0.2,200) 

def f_alpha_function(omegac):
    return alpha_function(epsi1,omegac,omega0,kappa_factor_omega0,kappa_r_factor)

def f_alpha_functionMIE(omegac):
    return alpha_functionMIE(omegac,epsi1,R)

#%%

title1 = r'$\epsilon_1$ = %i, $\kappa$ = %.2f$\omega_0$, $\kappa_r$ = %.2f$\kappa$, $\omega_0$ =%.2fTHz' %(epsi1,kappa_factor_omega0,kappa_r_factor,omega0THz)   
title2 = r'$\epsilon_1$ = %i, R = %inm' %(epsi1,R*1e-3)   

labelx = r'$\omega$ [THz]'

tamfig = (4.5,3.5)
tamlegend = 10
tamletra = 11
tamtitle = 10
tamnum = 9
labelpady = 2
labelpadx = 2
pad = 0

def graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad):
    plt.figure(figsize=tamfig)
    plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
    plt.ylabel(labely,fontsize=tamletra,labelpad =labelpadx)
    plt.tick_params(labelsize = tamnum, pad = pad)
    plt.title(title,fontsize=int(tamtitle*0.9))

    return   

color_rgb1 = (1.0, 0.8, 0.6)
color_rgb2 = (0.99, 0.74, 0.71)
color_rgb3 = (0.529, 0.992, 0.945)
color_rgb4 = (0.41, 0.21, 0.61)

#%%

list_alpha_re = []
list_alphaMie_re = []
list_alpha_im = []
list_alphaMie_im = []

for OmegaTHz in list_OmegaTHz: 

    OmegaTHz = np.round(OmegaTHz,10)
    omegac = OmegaTHz*1e12/c

    rta1 = f_alpha_function(omegac)
    
    list_alpha_re.append(rta1.real)
    list_alpha_im.append(rta1.imag)

    rta2 = f_alpha_functionMIE(omegac)
    
    list_alphaMie_re.append(rta2.real)
    list_alphaMie_im.append(rta2.imag) 

#%%
print('Graficar el green tensor')
    
graph(title1,labelx,'Re(polarizability)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(list_OmegaTHz,list_alpha_re,'k-')
if save_graphs==1:
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'Re_polarizability' + '.png', format='png')

graph(title1,labelx,'Im(polarizability)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(list_OmegaTHz,list_alpha_im,'k-')
# plt.yscale('log')
if save_graphs==1:
    plt.tight_layout()
    plt.savefig( 'Im_polarizability' + '.png', format='png')

#%%

graph(title2,labelx,'Re(polarizability Mie)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(list_OmegaTHz,list_alphaMie_re,'k-')
if save_graphs==1:
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'Re_polarizabilityQE' + '.png', format='png')

graph(title2,labelx,'Im(polarizability Mie)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(list_OmegaTHz,list_alphaMie_im,'k-')
# plt.yscale('log')
if save_graphs==1:
    plt.tight_layout()
    plt.savefig( 'Im_polarizabilityQE' + '.png', format='png')

#%%
