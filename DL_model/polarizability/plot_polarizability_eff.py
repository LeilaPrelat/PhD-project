#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: leila

graficar polarizabilidad effectiva normalizada por alfa
vs omega
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
    from function_polarizability import alpha_function_eff
except ModuleNotFoundError:
    print(err)

try:
    sys.path.insert(1, path_basic)
    from function_polarizability import alpha_function_eff_QE
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
int_v = 400
v = c/int_v
b = -2 #electron position in z < 0 (arriba del plano)
xD = 2
epsi1 = 1
epsi2 = 1

hbmu,hbgama = 0.7,0.0001
zp = 2

omega0THz = 0.05
omega0 = omega0THz*1e12 
kappa_factor_omega0 = 0.1
kappa_r_factor= 0.5

list_OmegaTHz = np.linspace(0.001,0.2,200) 

def f_alpha_function(omegac):
    return alpha_function_eff(omegac,epsi1,epsi2,hbmu,hbgama,int_v,b,zp,xD,omega0,kappa_factor_omega0,kappa_r_factor)

def f_alpha_functionQE(omegac):
    return alpha_function_eff_QE(omegac,epsi1,epsi2,hbmu,hbgama,int_v,b,zp,xD,omega0,kappa_factor_omega0,kappa_r_factor)

#%%

title1 = r'$\epsilon_1$ = %i, v = c/%i $\mu$m/s, xD = %i$\mu$m, b = %.2f$\mu$m' %(epsi1,int_v,xD,b) 
labelx = r'$\omega$ [THz]'
labely = r'$\alpha_{eff}/\alpha$'
title2A = r'$\epsilon_2$ = %i, $\hbar\mu$ = %.2feV, $\hbar\gamma$ = %.4feV, $z_p=%.2f\mu$m' %(epsi2,hbmu,hbgama,zp) 
title3 = r'$\kappa$ = %.2f$\omega_0$, $\kappa_r$ = %.2f$\kappa$, $\omega_0$ =%.2fTHz' %(kappa_factor_omega0,kappa_r_factor,omega0THz)   
title = title1 + '\n'  + title2A  + '\n' + title3

tamfig = (4.5,3.5)
tamlegend = 10
tamletra = 11
tamtitle = 10
tamnum = 9
labelpady = 2
labelpadx = 2
pad = 0
ms = 4
mk = 2
hp = 0.3

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
list_alphaQE_re = []
list_alpha_im = []
list_alphaQE_im = []

j = 0
for OmegaTHz in list_OmegaTHz: 

    OmegaTHz = np.round(OmegaTHz,10)
    omegac = OmegaTHz*1e12/c

    rta1 = f_alpha_function(omegac)
    
    list_alpha_re.append(rta1.real)
    list_alpha_im.append(rta1.imag)

    rta2 = f_alpha_functionQE(omegac)
    
    list_alphaQE_re.append(rta2.real)
    list_alphaQE_im.append(rta2.imag) 

    j = j + 1
    print(j)

#%%
print('Graficar el green tensor')
    
graph(title,labelx,'Re ' + labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(list_OmegaTHz,list_alpha_re,'--',lw = 2,color = 'darkviolet',label = 'numerical')
plt.plot(list_OmegaTHz,list_alphaQE_re,'.',ms = ms,color = 'lightseagreen',label = 'numerical QE')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
if save_graphs==1:
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'Re_polarizability_eff' + '.png', format='png')

graph(title,labelx,'Im ' + labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(list_OmegaTHz,list_alpha_im,lw = 2,color = 'darkviolet',label = 'numerical')
plt.plot(list_OmegaTHz,list_alphaQE_im,'.',ms = ms,color = 'lightseagreen',label = 'numerical QE')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
if save_graphs==1:
    plt.tight_layout()
    plt.savefig( 'Im_polarizability_eff' + '.png', format='png')

#%%
