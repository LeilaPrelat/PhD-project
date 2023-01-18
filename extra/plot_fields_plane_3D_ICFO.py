#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: leila

graficar campos
seccion 4.4 del overleaf

"""

import numpy as np
import sys
import os 
import matplotlib.pyplot as plt

#%%

save_graphs = 1 #guardar los graficos 2D del campo

list_field = ['Etot', 'Htot'] 
type_field = list_field[1] 

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_save = path_basic + '/' + 'fields_plane_ICFO'

if save_graphs==1:    
    if not os.path.exists(path_save):
        print('Creating folder to save graphs')
        os.mkdir(path_save)

try:
    sys.path.insert(1, path_basic)
    from constants_plane import constantes
except ModuleNotFoundError:
    print('constants_plane.py no se encuentra en ' + path_basic)

pi,hb,c,alfac,hbargama,mu1,mu2,epsi1,epsi2 = constantes()
aux = c*hb

print('Definir parametros del problema')

hbarmu = 0.9        # eV mu_c
kx = 0.05
energy = 0.0924 # eV
energy = 0.003 # eV
omegac = energy/aux # 1/micrones

Ainc = 1
Binc = 1

ref_ind = np.sqrt(epsi1*mu1)
if omegac > kx : 
    raise TypeError('omega outside light cone')

z_plane = 0.5
cota = 2.5 # limite y, z en micrones

#%%

#print('Importar modulos necesarios para este codigo')
err = 'fields_plane_Fresnel_ICFO.py no se encuentra en ' + path_basic
if type_field == 'Etot':
    try:
        sys.path.insert(1, path_basic)
        from fields_plane_Fresnel_ICFO import fieldE_plane
    except ModuleNotFoundError:
        print(err)

    def fields(x,z):
        Ex,Ey,Ez = fieldE_plane(omegac,hbarmu,kx,z_plane,x,-z,Ainc,Binc)
        return np.abs(Ex)**2 + np.abs(Ey)**2 + np.abs(Ez)**2
    labelz = '|Etot|$^2$'
else: 
    try:
        sys.path.insert(1, path_basic)
        from fields_plane_Fresnel_ICFO import fieldH_plane
    except ModuleNotFoundError:
        print(err)

    def fields(x,z):
        Hx,Hy,Hz = fieldH_plane(omegac,hbarmu,kx,z_plane,x,-z,Ainc,Binc)
        return np.abs(Hx)**2 + np.abs(Hy)**2 + np.abs(Hz)**2
    labelz = '|Htot|$^2$'

#%%
    
n2 = 100
listx = np.linspace(-cota,cota,n2)
listy = np.linspace(-cota,cota,n2)
X, Y = np.meshgrid(listx, listy)
limits = [min(listx) , max(listx), min(listy) , max(listy)]

labelx,labely = 'x --> [$\mu$m]', '<--  z [$\mu$m]'
labelpng = '_kx%.4f_mu%.4f' %(kx,hbarmu)

def graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad):
    plt.figure(figsize=tamfig)
    plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
    plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
    plt.tick_params(labelsize = tamnum, pad = pad)
    plt.title(title,fontsize=int(tamtitle*0.9))
    plt.plot(listx,np.ones(n2)*z_plane, '-b',lw = 1)
    return   

tamfig = (4.5,3.5)
tamlegend = 10
tamletra = 11
tamtitle = 10
tamnum = 9
labelpady = -1.5
labelpadx = -0.5
pad = 0.5

#%%

title1 = r'kx = %.2f 1/$\mu$m, $\mu_c$ = %.1feV, E = %.2feV' %(kx,hbarmu,energy)
title2 = r'$z_p$ = %.1f$\mu$m' %(z_plane)
title3 = r'$\epsilon_1$ = %i, $\mu_1$ = %i, $\epsilon_2$ = %i, $\mu_2$ = %i' %(epsi1,mu1,epsi2,mu2)
title = title1 + '\n' + title2 + ', ' + title3

print('Graficar el campo ' + labelz + ' para el medio 1 y 2')
    
f1 = np.vectorize(fields)
Z1 = f1(X, Y)
    
graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
im = plt.imshow(Z1, extent = limits, cmap=plt.cm.hot, interpolation='bilinear')
cbar = plt.colorbar(im)
cbar.ax.tick_params(labelsize = tamnum)
cbar.set_label(labelz,fontsize=tamlegend,labelpad = 1)

if save_graphs==1:
    plt.tight_layout(1)
    os.chdir(path_save)
    plt.savefig(type_field + labelpng + '.png', format='png')

del Z1

#%%
