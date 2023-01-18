#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: leila

graficar campos

"""

import numpy as np
import sys
import os 
import matplotlib.pyplot as plt

#%%

save_graphs = 1 #guardar los graficos 2D del campo
modulo = 1 # if modulo == 1 : grafica |Ex| sino grafica Re(Ex)

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_save = path_basic + '/' + 'fields_Eind_ICFO'

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

hbarmu = 0.3        # eV mu_c
# kx = 0.05
# energy = 0.0924 # eV
# energy = 0.003 # eV
# omegac = energy/aux # 1/micrones

theta_inc = np.pi/4 # no cambiar

# ref_ind = np.sqrt(epsi1*mu1)
# if omegac > kx : 
#     raise TypeError('omega outside light cone')

z_electron = -0.5
z_plane = 0.5
cota = 2.5 # limite y, z en micrones
px,py,pz = 1,1,1
epsi_h = 6.9

x,y = 0,0 # no cambiar

#%%

#print('Importar modulos necesarios para este codigo')
err = 'fields_Eind_ICFO.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from fields_Eind_ICFO import fieldE_ind_dip_x
except ModuleNotFoundError:
    print(err)

def fields(omegac,z_electron):
    Ex = fieldE_ind_dip_x(omegac,hbarmu,theta_inc,z_plane,z_electron,x,y,epsi_h,px,py,pz)
    if modulo == 1:
        return np.abs(Ex)
    else:
        return Ex.real

#%%
    
n2 = 100
listx = np.linspace(-cota,cota,n2)
listy = np.linspace(-cota,cota,n2)
X, Y = np.meshgrid(listx, listy)
limits = [min(listx) , max(listx), min(listy) , max(listy)]

labelx,labely = '$\omega/c$ [1/$\mu$m]', '$z_{electron}$ [$\mu$m]'
labelpng = '_mu%.4f' %(hbarmu)
if modulo == 1:
    labelz = '|Ex|'
else:
    labelz = 'Re(Ex)'

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

title1 = r'$\theta_{inc}$ = $\pi/4$, $\mu_c$ = %.1feV, x = 0, z = 0' %(hbarmu)
title2 = r'$z_p$ = %.1f$\mu$m, px = %i, py = %i, pz = %i' %(z_plane,px,py,pz)
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
    plt.savefig(labelpng + '.png', format='png')

del Z1

#%%
