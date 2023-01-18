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

save_graphs = 1 #guardar los graficos 2D del campo
graphs_vs_z = 1 #graficar vs z
graphs_vs_omegaTHz = 0 #graficar vs omega en THz
sns.set()

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('GreenTensor/' + 'reflected','')
path_save = path_basic + '/' + 'green_tensor_reflected_Ag'

if save_graphs==1:    
    if not os.path.exists(path_save):
        print('Creating folder to save graphs')
        os.mkdir(path_save)

err = 'green_tensor_ref_Ag.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from green_tensor_ref_Ag import green_tensor_NUM_QE, green_tensor_NUM
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

zD = 0
xD, yD = 0,0
x,y = 0,0
#omega = 0.7*1e12
epsi1 = 1
epsi2 = 1
energy_bulk = 9.17
omega_bulk, hbar_gamma_in, d = energy_bulk/hb, 0.0001, 0.01 
#hbmu, hbgama = 0.3, 0.0001 #potencial quimico in ev, freq de collision in ev
zp = 2 #micrones # posicion del plano (z apunta hacia abajo)

########################
omegaTHz = 0.3
omegac = omegaTHz*1e12/c 
list_z = np.linspace(3,53,201)

########################
z0 = 20
list_OmegaTHz = np.linspace(0.5,1.5,201)

#%%
    
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

color_rgb1 = (1.0, 0.01, 0.24)
color_rgb2 = (0.8, 0.58, 0.46)
color_rgb3 = (0.55, 0.71, 0)
color_rgb4 = (0.6, 0.73, 0.89)

#%%

title1 = r'$\epsilon_1$ = %i, $\epsilon_2$ = %i, x = %i$\mu$m, y = %i$\mu$m' %(epsi1,epsi2,x,y)
title3 = r'$E_{bulk}$ = %.2feV, $\hbar \gamma$ = %.4feV, d = %.2f$\mu$m, zp = %i$\mu$m' %(energy_bulk, hbar_gamma_in,d,zp)
 
listI05_re = []
listI25_re = []
listI06_re = []
listI26_re = []

listI05_v3_re = []
listI25_v3_re = []
listI06_v3_re = []
listI26_v3_re = []

listI05_im = []
listI25_im = []
listI06_im = []
listI26_im = []

listI05_v3_im = []
listI25_v3_im = []
listI06_v3_im = []
listI26_v3_im = []

if graphs_vs_z == 1:
    
    title2A = r'$\omega$ = %.2fTHz, $x_D$ = %i$\mu$m, $y_D$ = %i$\mu$m, $z_D$ = %i$\mu$m' %(omegaTHz,xD,yD,zD)
    title = title1 + '\n' + title2A + '\n' + title3
    list_x = list_z
    labelx = '$z$ [$\mu$m]'
    label1 = '_vs_z'

    for z in list_z: 
        
        z = np.round(z,3)
        
        f1 = green_tensor_NUM_QE(omegac,epsi1,epsi2,omega_bulk,hbar_gamma_in,d,x,y,z,xD,yD,zD,zp)
        f2 = green_tensor_NUM(omegac,epsi1,epsi2,omega_bulk,hbar_gamma_in,d,x,y,z,xD,yD,zD,zp)
        
        I0_5, I2_5, I0_6, I2_6 = f1
        I0_5_v3, I2_5_v3, I0_6_v3, I2_6_v3 = f2
        
        listI05_re.append(I0_5.real)
        listI25_re.append(I2_5.real)
        listI06_re.append(I0_6.real)
        listI26_re.append(I2_6.real)
        
        listI05_im.append(I0_5.imag)
        listI25_im.append(I2_5.imag)
        listI06_im.append(I0_6.imag)
        listI26_im.append(I2_6.imag)
        
        listI05_v3_re.append(I0_5_v3.real)
        listI25_v3_re.append(I2_5_v3.real)
        listI06_v3_re.append(I0_6_v3.real)
        listI26_v3_re.append(I2_6_v3.real)
        
        listI05_v3_im.append(I0_5_v3.imag)
        listI25_v3_im.append(I2_5_v3.imag)
        listI06_v3_im.append(I0_6_v3.imag)
        listI26_v3_im.append(I2_6_v3.imag)

elif graphs_vs_omegaTHz == 1 :

    title2B = r'z = %i$\mu$m, $x_D$ = %i$\mu$m, $y_D$ = %i$\mu$m, $z_D$ = %i$\mu$m' %(z0,xD,yD,zD)
    title = title1 + '\n' + title2B + '\n' + title3
    list_x = list_OmegaTHz
    labelx = '$\omega$ [THz]'
    label1 = '_vs_OmegaTHz'

    for OmegaTHz in list_OmegaTHz: 
        
        OmegaTHz = np.round(OmegaTHz,3)
        Omegac = OmegaTHz*1e12/c
        
        f1 = green_tensor_NUM_QE(Omegac,epsi1,epsi2,omega_bulk,hbar_gamma_in,d,x,y,z0,xD,yD,zD,zp)
        f2 = green_tensor_NUM(Omegac,epsi1,epsi2,omega_bulk,hbar_gamma_in,d,x,y,z0,xD,yD,zD,zp)
        
        I0_5, I2_5, I0_6, I2_6 = f1
        I0_5_v3, I2_5_v3, I0_6_v3, I2_6_v3 = f2
        
        listI05_re.append(I0_5.real)
        listI25_re.append(I2_5.real)
        listI06_re.append(I0_6.real)
        listI26_re.append(I2_6.real)
        
        listI05_im.append(I0_5.imag)
        listI25_im.append(I2_5.imag)
        listI06_im.append(I0_6.imag)
        listI26_im.append(I2_6.imag)
        
        listI05_v3_re.append(I0_5_v3.real)
        listI25_v3_re.append(I2_5_v3.real)
        listI06_v3_re.append(I0_6_v3.real)
        listI26_v3_re.append(I2_6_v3.real)
        
        listI05_v3_im.append(I0_5_v3.imag)
        listI25_v3_im.append(I2_5_v3.imag)
        listI06_v3_im.append(I0_6_v3.imag)
        listI26_v3_im.append(I2_6_v3.imag)

#%%

print('Graficar el green tensor vs z')
    
graph(title,labelx,r'Re $I^{(0)}_7$',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(list_x,listI05_re,'.',ms = 7,color = color_rgb1,label = 'numerical QE')
# plt.yscale('log')
plt.plot(list_x,listI05_v3_re,'-',color = 'green',label = 'numerical')
plt.legend(loc = 'best',markerscale=1,fontsize=tamlegend,frameon=True,handletextpad=0.5)
if save_graphs==1:
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'green_tensor_Re_I07' + label1 + '.png', format='png')

graph(title,labelx,r'Im $I^{(0)}_7$',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(list_x,listI05_im,'.',ms = 7,color = color_rgb1,label = 'numerical QE')
# plt.yscale('log')
plt.plot(list_x,listI05_v3_im,'-',color = 'green',label = 'numerical')
plt.legend(loc = 'best',markerscale=1,fontsize=tamlegend,frameon=True,handletextpad=0.5)
if save_graphs==1:
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'green_tensor_Im_I07' + label1 + '.png', format='png')

#############################################################################################
graph(title,labelx,r'Re $I^{(2)}_7$',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(list_x,listI25_re,'.',ms = 7,color = color_rgb2,label = 'numerical QE')
# plt.yscale('log')
plt.plot(list_x,listI25_v3_re,'-',color = 'green',label = 'numerical')
plt.legend(loc = 'best',markerscale=1,fontsize=tamlegend,frameon=True,handletextpad=0.5)
if save_graphs==1:
    plt.tight_layout()
    plt.savefig( 'green_tensor_Re_I27' + label1 + '.png', format='png')

graph(title,labelx,r'Im $I^{(2)}_7$',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(list_x,listI25_im,'.',ms = 7,color = color_rgb2,label = 'numerical QE')
# plt.yscale('log')
plt.plot(list_x,listI25_v3_im,'-',color = 'green',label = 'numerical')
plt.legend(loc = 'best',markerscale=1,fontsize=tamlegend,frameon=True,handletextpad=0.5)
if save_graphs==1:
    plt.tight_layout()
    plt.savefig( 'green_tensor_Im_I27' + label1 + '.png', format='png')
#############################################################################################
graph(title,labelx,r'Re $I^{(0)}_8$',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(list_x,listI06_re,'.',ms = 7,color = color_rgb3,label = 'numerical QE')
plt.plot(list_x,listI06_v3_re,'-',color = 'green',label = 'numerical')
plt.legend(loc = 'best',markerscale=1,fontsize=tamlegend,frameon=True,handletextpad=0.5)
if save_graphs==1:
    plt.tight_layout()
    plt.savefig( 'green_tensor_Re_I08'+ label1 + '.png', format='png')

graph(title,labelx,r'Im $I^{(0)}_8$',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(list_x,listI06_im,'.',ms = 7,color = color_rgb3,label = 'numerical QE')
plt.plot(list_x,listI06_v3_im,'-',color = 'green',label = 'numerical')
plt.legend(loc = 'best',markerscale=1,fontsize=tamlegend,frameon=True,handletextpad=0.5)
if save_graphs==1:
    plt.tight_layout()
    plt.savefig( 'green_tensor_Im_I08'+ label1 + '.png', format='png')
############################################################################################
graph(title,labelx,r'Re $I^{(2)}_8$',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(list_x,listI26_re,'.',ms = 7,color = color_rgb4,label = 'numerical QE')
# plt.yscale('log')
plt.plot(list_x,listI26_v3_re,'-',color = 'green',label = 'numerical')
plt.legend(loc = 'best',markerscale=1,fontsize=tamlegend,frameon=True,handletextpad=0.5)
if save_graphs==1:
    plt.tight_layout()
    plt.savefig( 'green_tensor_Re_I28'+ label1 + '.png', format='png')

graph(title,labelx,r'Im $I^{(2)}_8$',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(list_x,listI26_im,'.',ms = 7,color = color_rgb4,label = 'numerical QE')
# plt.yscale('log')
plt.plot(list_x,listI26_v3_im,'-',color = 'green',label = 'numerical')
plt.legend(loc = 'best',markerscale=1,fontsize=tamlegend,frameon=True,handletextpad=0.5)
if save_graphs==1:
    plt.tight_layout()
    plt.savefig( 'green_tensor_Im_I28'+ label1 + '.png', format='png')

#%%
