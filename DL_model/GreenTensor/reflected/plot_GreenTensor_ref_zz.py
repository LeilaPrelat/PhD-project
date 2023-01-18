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

#%%

list_type_plane = ['graphene','Ag'] ### revisar conductividad de Ag antes de hacer el green tensor reflected
type_plane = list_type_plane[0]

save_graphs = 1 #guardar los graficos 2D del campo


graphs_vs_z = 1 #graficar vs z
graphs_vs_omegaTHz = 0 #graficar vs omega en THz
sns.set()

if graphs_vs_z + graphs_vs_omegaTHz > 1:
    raise TypeError('Choose 1')

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('GreenTensor/' + 'reflected','')

if type_plane == 'graphene':
    err = 'green_tensor_ref_zz_graphene.py no se encuentra en ' + path_basic
    try:
        sys.path.insert(1, path_basic)
        from green_tensor_ref_zz_graphene import green_tensor_NUM_QE, green_tensor_NUM
    except ModuleNotFoundError:
        print(err)
    path_save = path_basic + '/' + 'green_tensor_ref_zz_graphene'  


else:
    err = 'green_tensor_ref_zz_Ag.py no se encuentra en ' + path_basic
    try:
        sys.path.insert(1, path_basic)
        from green_tensor_ref_zz_Ag import green_tensor_NUM_QE, green_tensor_NUM
    except ModuleNotFoundError:
        print(err)
    path_save = path_basic + '/' + 'green_tensor_ref_zz_Ag'

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
xD, yD = 0,0
x,y = 0,0
#omega = 0.7*1e12
epsi1 = 1
epsi2 = 1
zp = 2 #micrones # posicion del plano (z apunta hacia abajo)
b = -2
title1 = r'$\epsilon_1$ = %i, $\epsilon_2$ = %i, b = %i$\mu$m, x = %i$\mu$m, y = %i$\mu$m' %(epsi1,epsi2,b,x,y)

if type_plane == 'graphene':
    hbmu, hbgama = 0.4, 0.0001 #potencial quimico in ev, freq de collision in ev
    title3 = r'$\hbar \mu$ = %.2feV, $\hbar \gamma$ = %.4feV, zp = %i$\mu$m' %(hbmu,hbgama,zp)
    
    def functionQE(omegac,z):
        return green_tensor_NUM_QE(omegac,epsi1,epsi2,hbmu,hbgama,x,y,z,xD,yD,zD,b,zp)
    
    def function(omegac,z):    
        return green_tensor_NUM(omegac,epsi1,epsi2,hbmu,hbgama,x,y,z,xD,yD,zD,b,zp)

    if save_graphs==1:    
        path_save = path_save + '/' + 'mu_%.2f' %(hbmu) 
        if not os.path.exists(path_save):
            print('Creating folder to save graphs')
            os.mkdir(path_save)
        
else:
    energy_bulk = 9.17
    omega_bulk, hbar_gamma_in, d = energy_bulk/hb, 0.0001, 0.09 
    title3 = r'$E_{bulk}$ = %.2feV, $\hbar \gamma$ = %.4feV, d = %.2f$\mu$m, zp = %i$\mu$m' %(energy_bulk,hbar_gamma_in,d,zp)

    def functionQE(omegac,z):
        return green_tensor_NUM_QE(omegac,epsi1,epsi2,omega_bulk,hbar_gamma_in,d,x,y,z,xD,yD,zD,zp)
    
    def function(omegac,z):    
        return green_tensor_NUM(omegac,epsi1,epsi2,omega_bulk,hbar_gamma_in,d,x,y,z,xD,yD,zD,zp)    
    
########################
if graphs_vs_z == 1:
    omegaTHz = 0.7
    omegaTHz = 1.511
    
    omegaTHz = 0.01
    # omegaTHz = 1.663 - 0.152j
    
    omegac0 = omegaTHz*1e12/c 
    list_z = np.linspace(-5,5,161)

########################
if graphs_vs_omegaTHz == 1:
    z0 = 5
    list_OmegaTHz = np.linspace(0.01,2.01,251)

#%%
    
tamfig = (4.5,3.5)
tamlegend = 10
tamletra = 11
tamtitle = 10
tamnum = 9
labelpady = -1.5
labelpadx = 0
pad = 0.5
mk = 2
ms = 4
hp = 0.3
hspace = 0
wspace = 0.01
loc1 = [0.3,0.935] 

def graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad):
    plt.figure(figsize=tamfig)
    plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
    plt.ylabel(labely,fontsize=tamletra,labelpad =labelpadx)
    plt.tick_params(labelsize = tamnum, pad = pad)
    plt.title(title,fontsize=int(tamtitle*0.9))
    return   

color_rgb1 = (1.0, 0.01, 0.24)
color_rgb2 = (0.8, 0.58, 0.46)
color_rgb3 = (0.55, 0.71, 0)
color_rgb4 = (0.6, 0.73, 0.89)

#%%

listI06_re = []
listI06_im = []

listI06_v3_re = []
listI06_v3_im = []

if graphs_vs_z == 1:
    
    if omegaTHz.imag == 0:
        title2A = r'$\omega$ = %.2fTHz, $x_D$ = %i$\mu$m, $y_D$ = %i$\mu$m, $z_D$ = %i$\mu$m' %(omegaTHz.real,xD,yD,zD)
    elif omegaTHz.imag <0:
        title2A = r'$\omega$ = (%.2f %.2f)THz, $x_D$ = %i$\mu$m, $y_D$ = %i$\mu$m, $z_D$ = %i$\mu$m' %(omegaTHz.real,omegaTHz.imag,xD,yD,zD)
    else:
        title2A = r'$\omega$ = (%.2f + %.2f)THz, $x_D$ = %i$\mu$m, $y_D$ = %i$\mu$m, $z_D$ = %i$\mu$m' %(omegaTHz.real,omegaTHz.imag,xD,yD,zD)
 
    title = title1 + '\n' + title2A + '\n' + title3
    list_x = list_z
    labelx = '$z$ [$\mu$m]'
    label1 = '_vs_z'

    for z in list_z: 
        
        z = np.round(z,8)
        
        I0_6 = functionQE(omegac0,z)
        I0_6_v3 = function(omegac0,z)
        
        listI06_re.append(I0_6.real)
        listI06_im.append(I0_6.imag)
        
        listI06_v3_re.append(I0_6_v3.real)
        listI06_v3_im.append(I0_6_v3.imag)

elif graphs_vs_omegaTHz == 1 :

    title2B = r'z = %i$\mu$m, $x_D$ = %i$\mu$m, $y_D$ = %i$\mu$m, $z_D$ = %i$\mu$m' %(z0,xD,yD,zD)
    title = title1 + '\n' + title2B + '\n' + title3
    list_x = list_OmegaTHz
    labelx = '$\omega$ [THz]'
    label1 = '_vs_OmegaTHz'

    for OmegaTHz in list_OmegaTHz: 
        
        OmegaTHz = np.round(OmegaTHz,8)
        Omegac = OmegaTHz*1e12/c
        
        I0_6 = functionQE(Omegac,z0)
        I0_6_v3 = function(Omegac,z0)
        
        listI06_re.append(I0_6.real)
        listI06_im.append(I0_6.imag)
        
        listI06_v3_re.append(I0_6_v3.real)
        listI06_v3_im.append(I0_6_v3.imag)

#%%

print('Graficar el green tensor ' + label1)
    
#############################################################################################
# if graphs_vs_z == 1:
graph(title,labelx,r'Re(Reflected zz Green Tensor)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(list_x,listI06_re,'.',ms = ms,color = 'black',label = 'numerical QE')
plt.plot(list_x,listI06_v3_re,'-',color = 'green',label = 'numerical')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
if save_graphs==1:
    os.chdir(path_save)
    plt.tight_layout()
    plt.savefig( 'green_tensor_zz_Re'+ label1 + '.png', format='png')
 
graph(title,labelx,r'Im(Reflected zz Green Tensor)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(list_x,listI06_im,'.',ms = ms,color = 'black',label = 'numerical QE')
plt.plot(list_x,listI06_v3_im,'-',color = 'green',label = 'numerical')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
if save_graphs==1:
    plt.tight_layout()
    plt.savefig( 'green_tensor_zz_Im'+ label1 + '.png', format='png')
############################################################################################

#%%
