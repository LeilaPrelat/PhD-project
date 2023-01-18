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

err = 'green_tensor_direct_numerical.py no se encuentra en ' + path_basic
err2 = 'green_tensor_direct_analytical.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from green_tensor_direct_numerical import green_tensor_functionDIR
except ModuleNotFoundError:
    print(err)

try:
    sys.path.insert(1, path_basic)
    from green_tensor_direct_analytical import green_tensor_functionDIR_2
except ModuleNotFoundError:
    print(err2)


try:
    sys.path.insert(1, path_basic)
    from constants import constantes
except ModuleNotFoundError:
    print('constants.py no se encuentra en ' + path_basic)

pi,hb,c,alfac,mu1,mu2 = constantes()
aux = c*hb

#%%

print('Definir parametros del problema')

zD = 0
xD, yD = 0,0
x,y = 1,1
#omega = 0.7*1e12
omegaTHz = 0.7
omegac = omegaTHz*1e12/c 
list_z = np.linspace(1,500,501)
epsi1 = 1

#%%
    
labelx = '$z$ [$\mu$m]'

tamfig = (4.5,3.5)
tamlegend = 10
tamletra = 11
tamtitle = 10
tamnum = 9
labelpady = -1.5
labelpadx = 0
pad = 0.5

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

title1 = r'$\epsilon_1$ = %i, x = %i$\mu$m, y = %i$\mu$m, $x_D$ = %i$\mu$m, $y_D$ = %i$\mu$m' %(epsi1,x,y,xD,yD)
title2 = r'$\omega$ = %.2f THz, $z_D$ = %i$\mu$m' %(omegaTHz,zD)
title = title1 + '\n' + title2 

listI05 = []
listI25 = []
listI06 = []
listI26 = []

listI05_v2 = []
listI25_v2 = []
listI06_v2 = []
listI26_v2 = []

for z in list_z: 
    I0_5, I2_5, I0_6, I2_6 = green_tensor_functionDIR(omegac,epsi1,x,y,z,xD,yD,zD)
    I0_5_v2, I2_5_v2, I0_6_v2, I2_6_v2 = green_tensor_functionDIR_2(omegac,epsi1,x,y,z,xD,yD,zD)
    
    listI05.append(I0_5)
    listI25.append(I2_5)
    listI06.append(I0_6)
    listI26.append(I2_6)

    listI05_v2.append(I0_5_v2)
    listI25_v2.append(I2_5_v2)
    listI06_v2.append(I0_6_v2)
    listI26_v2.append(I2_6_v2)

print('Graficar el green tensor')
    
graph(title,labelx,r'$I^{(0)}_5$',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(list_z,listI05,'.',color = color_rgb1,label = 'numerical')
plt.plot(list_z,listI05_v2,'-',color = 'black',label = 'analytical QE + asymptotic')
plt.legend(loc = 'best',markerscale=1,fontsize=tamlegend,frameon=True,handletextpad=0.1)
if save_graphs==1:
    plt.tight_layout(1)
    os.chdir(path_save)
    plt.savefig( 'green_tensor_I05' + '.png', format='png')

graph(title,labelx,r'$I^{(2)}_5$',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(list_z,listI25,'.',color = color_rgb2,label = 'numerical')
plt.plot(list_z,listI25_v2,'-',color = 'black',label = 'analytical QE')
plt.legend(loc = 'best',markerscale=1,fontsize=tamlegend,frameon=True,handletextpad=0.1)
if save_graphs==1:
    plt.tight_layout(1)
    plt.savefig( 'green_tensor_I25' + '.png', format='png')

graph(title,labelx,r'$I^{(0)}_6$',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(list_z,listI06,'.',color = color_rgb3,label = 'numerical')
plt.plot(list_z,listI06_v2,'-',color = 'black',label = 'analytical QE')
plt.legend(loc = 'best',markerscale=1,fontsize=tamlegend,frameon=True,handletextpad=0.1)
if save_graphs==1:
    plt.tight_layout(1)
    plt.savefig( 'green_tensor_I06' + '.png', format='png')

graph(title,labelx,r'$I^{(2)}_6$',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(list_z,listI26,'.',color = color_rgb4,label = 'numerical')
plt.plot(list_z,listI26_v2,'-',color = 'black',label = 'analytical QE')
plt.legend(loc = 'best',markerscale=1,fontsize=tamlegend,frameon=True,handletextpad=0.1)
if save_graphs==1:
    plt.tight_layout(1)
    plt.savefig( 'green_tensor_I26' + '.png', format='png')

#%%
