#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: leila

graficar green tensor Gnn con la misma convencion 
que el paper 370 (z hacia abajo, zplane>0, particula en r = 0, medio 1 arriba del plano, 
                  medio 2 abajo del plano)
"""

import numpy as np
import sys
import os 
import matplotlib.pyplot as plt

#%%

save_graphs = 1 #guardar los graficos 2D del campo

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_save = path_basic + '/' + 'green_tensor'

if save_graphs==1:    
    if not os.path.exists(path_save):
        print('Creating folder to save graphs')
        os.mkdir(path_save)

err = 'green_tensor_ICFO.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from green_tensor import green_tensor_function
except ModuleNotFoundError:
    print(err)

try:
    sys.path.insert(1, path_basic)
    from constants import constantes
except ModuleNotFoundError:
    print('constants.py no se encuentra en ' + path_basic)

pi,hb,c,alfac,mu1,mu2 = constantes()
aux = c*hb

#%%

print('Definir parametros del problema')

hbarmu = 0.3        # eV mu_c
# kx = 0.05
# energy = 0.0924 # eV
# energy = 0.003 # eV
# omegac = energy/aux # 1/micrones

# ref_ind = np.sqrt(epsi1*mu1)
# if omegac > kx : 
#     raise TypeError('omega outside light cone')

z_plane = 0.5
list_omegac = np.linspace(0.1,0.5,101)
px,py,pz = 1,1,1
epsi_h = 6.9
epsi1 = 1
epsi2 = 1
hbargama = 0.0001      # collision frequency in eV

#%%
    
labelx,labely1,labely2 = '$\omega/c$ [1/$\mu$m]', 'Re Green Tensor', 'Im Green Tensor'
labelpng = '_mu%.4f' %(hbarmu)

def graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad):
    plt.figure(figsize=tamfig)
    plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
    plt.ylabel(labely,fontsize=tamletra,labelpad =labelpadx)
    plt.tick_params(labelsize = tamnum, pad = pad)
    plt.title(title,fontsize=int(tamtitle*0.9))
    return   

tamfig = (4.5,3.5)
tamlegend = 10
tamletra = 11
tamtitle = 10
tamnum = 9
labelpady = -1.5
labelpadx = 0
pad = 0.5

#%%

title1 = r'$\mu_c$ = %.1feV, $\gamma_c$ = %.4f eV, $z_p$ = %.1f$\mu$m' %(hbarmu,hbargama,z_plane)
title2 = r'px = %i, py = %i, pz = %i' %(px,py,pz)
title3 = r'$\epsilon_1$ = %i, $\mu_1$ = %i, $\epsilon_2$ = %i, $\mu_2$ = %i, $\epsilon_h$ = %.1f' %(epsi1,mu1,epsi2,mu2,epsi_h)
title = title1 + '\n' + title2 + ', ' + title3

listy1 = []
listy2 = []
for omegac in list_omegac: 
    valuey = green_tensor_function(omegac,hbarmu,z_plane,epsi_h,epsi1,epsi2,hbargama,px,py,pz)
    listy1.append(valuey.real)
    listy2.append(valuey.imag)

print('Graficar el green tensor')
    
graph(title,labelx,labely1,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(list_omegac,listy1)

if save_graphs==1:
    plt.tight_layout(1)
    os.chdir(path_save)
    plt.savefig( 'green_tensor_Re' +  labelpng + '.png', format='png')

graph(title,labelx,labely2,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(list_omegac,listy2)

if save_graphs==1:
    plt.tight_layout(1)
    os.chdir(path_save)
    plt.savefig( 'green_tensor_Im' +labelpng + '.png', format='png')

#%%
