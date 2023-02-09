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
import seaborn as sns

sns.set()

#%%

save_graphs = 1 #guardar los graficos 2D del campo
    
#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('GreenTensor/' + 'direct','')
path_save = path_basic + '/' + 'green_tensor_ref_functions'

if save_graphs==1:    
    if not os.path.exists(path_save):
        print('Creating folder to save graphs')
        os.mkdir(path_save)

err = 'green_tensor_direct_numerical.py no se encuentra en ' + path_basic
err2 =  'green_tensor_direct_QE_numerical.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from green_tensor_ref_functions import functions_GreenTensor
except ModuleNotFoundError:
    print(err)

try:
    sys.path.insert(1, path_basic)
    from green_tensor_ref_QE_functions import functions_QE_GreenTensor
except ModuleNotFoundError:
    print(err2)

try:
    sys.path.insert(1, path_constants)
    from constants import constantes
except ModuleNotFoundError:
    print('constants.py no se encuentra en ' + path_constants)

pi,hb,c,alfac,mu1,mu2 = constantes()
aux = c*hb

#%%

print('Definir parametros del problema')

zD = 0
xD, yD, zD = 0,0,0
x,y,z = 1,1,1
#omega = 0.7*1e12
omegaTHz = 0.7
omegac = omegaTHz*1e12/c 
epsi1 = 1
epsi2 = 1
zp = 2 #micrones # posicion del plano (z apunta hacia abajo)
hbmu, hbgama = 0.4, 0.0001 #potencial quimico in ev, freq de collision in ev

list_alpha_parallel1 = np.linspace(0,0.95,99)
list_alpha_parallel2 = np.linspace(1.5,60,73)
list_alpha_parallel2 = np.linspace(1.5,1001.5,201)

list_alpha_parallel3 = np.linspace(1,1001,201) #for QE

def function(alpha_parallel):
    return functions_GreenTensor(omegac,epsi1,epsi2,hbmu,hbgama,alpha_parallel,x,y,z,xD,yD,zD,zp)

def functionQE(alpha_parallel):
    return functions_QE_GreenTensor(omegac,epsi1,epsi2,hbmu,hbgama,alpha_parallel,x,y,z,xD,yD,zD,zp)

#%%
    
labelx = r'$\alpha_\parallel$'

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

title1 = r'$\epsilon_1$ = %i, $\epsilon_2$ = %i, x = %i$\mu$m, y = %i$\mu$m, $x_D$ = %i$\mu$m, $y_D$ = %i$\mu$m' %(epsi1,epsi2,x,y,xD,yD)
title2 = r'$\omega$ = %.2f THz, z = %i$\mu$m, $z_D$ = %i$\mu$m' %(omegaTHz,z,zD)
title3 = r'$\hbar \mu$ = %.2feV, $\hbar \gamma$ = %.4feV, zp = %i$\mu$m' %(hbmu,hbgama,zp)
title = title1 + '\n' + title2 + '\n ' + title3

listItot_re = []
listItotQE_re = []
listItot_im = []
listItotQE_im = []
list_alpha_parallel_tot = []

for alpha_parallel in list_alpha_parallel1: 
    rta = function(alpha_parallel)
    
    listItot_re.append(rta.real)
    listItot_im.append(rta.imag)
    alpha_parallel = np.round(alpha_parallel,10)
    list_alpha_parallel_tot.append(alpha_parallel)

del rta, alpha_parallel  
for alpha_parallel in list_alpha_parallel2: 
    rta = function(alpha_parallel)
    
    listItot_re.append(rta.real)
    listItot_im.append(rta.imag) 
    alpha_parallel = np.round(alpha_parallel,10)
    list_alpha_parallel_tot.append(alpha_parallel)

for alpha_parallel in list_alpha_parallel3: 
    rta = functionQE(alpha_parallel)
    
    listItotQE_re.append(rta.real)
    listItotQE_im.append(rta.imag) 

#%%
print('Graficar el green tensor')
    
graph(title,labelx,'Re(Inside of Green Tensor)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(list_alpha_parallel_tot,listItot_re,'k-')
if save_graphs==1:
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'Re_green_tensor_ref' + '.png', format='png')

graph(title,labelx,'Im(Inside of Green Tensor)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(list_alpha_parallel_tot,listItot_im,'k-')
# plt.yscale('log')
if save_graphs==1:
    plt.tight_layout()
    plt.savefig( 'Im_green_tensor_ref' + '.png', format='png')

#%%

graph(title,labelx,'Re(Inside of Green Tensor QE)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(list_alpha_parallel3,listItotQE_re,'k-')
if save_graphs==1:
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'Re_green_tensorQE_ref' + '.png', format='png')

graph(title,labelx,'Im(Inside of Green Tensor QE)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(list_alpha_parallel3,listItotQE_im,'k-')
# plt.yscale('log')
if save_graphs==1:
    plt.tight_layout()
    plt.savefig( 'Im_green_tensorQE_ref' + '.png', format='png')

#%%
