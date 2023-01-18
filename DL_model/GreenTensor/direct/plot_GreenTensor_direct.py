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
plot_numerical_sinQE = 1
plot_vs_z = 0 #graficar vs z
plot_vs_omegaTHz = 0 #graficar vs omega en THz
plot_vs_xD = 1 # for the integral 

sns.set()    

if plot_vs_z + plot_vs_omegaTHz + plot_vs_xD > 1:
    raise TypeError('Choose one option')

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('GreenTensor/' + 'direct','')
path_save = path_basic + '/' + 'green_tensor_direct'

if save_graphs==1:    
    if not os.path.exists(path_save):
        print('Creating folder to save graphs')
        os.mkdir(path_save)

err = 'green_tensor_direct_numerical.py no se encuentra en ' + path_basic
err2 = 'green_tensor_direct_analytical.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from green_tensor_direct_numerical import green_tensor_NUM_QE, green_tensor_NUM
except ModuleNotFoundError:
    print(err)

try:
    sys.path.insert(1, path_basic)
    from green_tensor_direct_analytical import green_tensor_ANA_QE
except ModuleNotFoundError:
    print(err2)


try:
    sys.path.insert(1, path_constants)
    from constants import constantes
except ModuleNotFoundError:
    print('constants.py no se encuentra en ' + path_constants)

pi,hb,c,alfac,mu1,mu2 = constantes()
aux = c*hb
aux2 = 1e12/c

#%%

print('Definir parametros del problema')

zD = 0
yD = 0
x,y = 10,10
#omega = 0.7*1e12
epsi1 = 1

#%%

tamfig = (4.5,3.5)
tamlegend = 10
tamletra = 11
tamtitle = 10
tamnum = 9
labelpady = 2
labelpadx = 2
pad = 0
mk = 2
ms = 4
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

########################
if plot_vs_z == 1:
#    omegaTHz0 = 0.7
#    omegaTHz0 = 1.511
    omegaTHz0 = 0.70
    xD0 = 0
    # omegaTHz = 1.663 - 0.152j
    
    omegac0 = omegaTHz0*aux2
    list_z = np.linspace(1,100,100)
    list_x = list_z
    labelx = '$z$ [$\mu$m]'
    label1 = '_vs_z'
    title2 = r'$\omega$ = %.2f THz, $x_D$ = %i$\mu$m' %(omegaTHz0,xD0)
########################
if plot_vs_omegaTHz == 1:
    z0 = 5
    xD0 = 0
    
    list_OmegaTHz = np.linspace(0.01,2.01,251)
    list_x = list_OmegaTHz 
    labelx = '$\omega$ [THz]'
    label1 = '_vs_OmegaTHz'
    title2 = r'z = %.2f $\mu$m, $x_D$ = %i$\mu$m' %(z0,xD0)

if plot_vs_xD == 1:
    omegaTHz0  = 0.70
    omegac0 = omegaTHz0*aux2
    z0 = 10
    
    list_xD = np.linspace(-80,80,321)
    list_xD  = np.delete(list_xD , 160) #sacar el cero porque diverge
    list_x = list_xD
    labelx = r'$x_D$ [$\mu$m]'
    label1 = '_vs_xD'
    title2 = r'$\omega$ = %.2f THz, z = %.2f $\mu$m' %(omegaTHz0,z0)

#%%

title1 = r'$\epsilon_1$ = %i, x = %i$\mu$m, y = %i$\mu$m, $y_D$ = %i$\mu$m, $z_D$ = %i$\mu$m' %(epsi1,x,y,yD,zD)
title = title1 + '\n' + title2 

listI05 = []
listI25 = []
listI06 = []
listI26 = []
total1_re = []
total1_im = []

listI05_v2 = []
listI25_v2 = []
listI06_v2 = []
listI26_v2 = []
total2_re = []
total2_im = []

listI05_v3 = []
listI25_v3 = []
listI06_v3 = []
listI26_v3 = []
total3_re = []
total3_im = []

if plot_vs_z == 1:
    for z in list_z: 
        I0_5, I2_5, I0_6, I2_6 = green_tensor_NUM_QE(omegac0,epsi1,x,y,z,xD0,yD,zD)
        I0_5_v2, I2_5_v2, I0_6_v2, I2_6_v2 = green_tensor_ANA_QE(omegac0,epsi1,x,y,z,xD0,yD,zD)
        if plot_numerical_sinQE == 1:
            I0_5_v3, I2_5_v3, I0_6_v3, I2_6_v3 = green_tensor_NUM(omegac0,epsi1,x,y,z,xD0,yD,zD)
        
        listI05.append(I0_5)
        listI25.append(I2_5)
        listI06.append(I0_6)
        listI26.append(I2_6)
        tot1 = I0_5 + I2_5 + I0_6 + I2_6
        total1_re.append(tot1.real)
        total1_im.append(tot1.imag)
    
        listI05_v2.append(I0_5_v2)
        listI25_v2.append(I2_5_v2)
        listI06_v2.append(I0_6_v2)
        listI26_v2.append(I2_6_v2)
        tot2 = I0_5_v2 + I2_5_v2 + I0_6_v2 + I2_6_v2
        total2_re.append(tot2.real)
        total2_im.append(tot2.imag)
    
        if plot_numerical_sinQE == 1:
            listI05_v3.append(I0_5_v3)
            listI25_v3.append(I2_5_v3)
            listI06_v3.append(I0_6_v3)
            listI26_v3.append(I2_6_v3)
            tot3 = I0_5_v3 + I2_5_v3 + I0_6_v3 + I2_6_v3
            total3_re.append(tot3.real)
            total3_im.append(tot3.imag)
    del z
if plot_vs_omegaTHz == 1:
    for OmegaTHz in list_OmegaTHz:

        omegac = OmegaTHz*aux2   
        
        I0_5, I2_5, I0_6, I2_6 = green_tensor_NUM_QE(omegac,epsi1,x,y,z0,xD0,yD,zD)
        I0_5_v2, I2_5_v2, I0_6_v2, I2_6_v2 = green_tensor_ANA_QE(omegac,epsi1,x,y,z0,xD0,yD,zD)
        if plot_numerical_sinQE == 1:
            I0_5_v3, I2_5_v3, I0_6_v3, I2_6_v3 = green_tensor_NUM(omegac,epsi1,x,y,z0,xD0,yD,zD)
        
        listI05.append(I0_5)
        listI25.append(I2_5)
        listI06.append(I0_6)
        listI26.append(I2_6)
        tot1 = I0_5 + I2_5 + I0_6 + I2_6
        total1_re.append(tot1.real)
        total1_im.append(tot1.imag)
    
        listI05_v2.append(I0_5_v2)
        listI25_v2.append(I2_5_v2)
        listI06_v2.append(I0_6_v2)
        listI26_v2.append(I2_6_v2)
        tot2 = I0_5_v2 + I2_5_v2 + I0_6_v2 + I2_6_v2
        total2_re.append(tot2.real)
        total2_im.append(tot2.imag)
    
        if plot_numerical_sinQE == 1:
            listI05_v3.append(I0_5_v3)
            listI25_v3.append(I2_5_v3)
            listI06_v3.append(I0_6_v3)
            listI26_v3.append(I2_6_v3)
            tot3 = I0_5_v3 + I2_5_v3 + I0_6_v3 + I2_6_v3
            total3_re.append(tot3.real)
            total3_im.append(tot3.imag)

    del omegac
if plot_vs_xD == 1:
    
    for xD in list_xD:
        I0_5, I2_5, I0_6, I2_6 = green_tensor_NUM_QE(omegac0,epsi1,x,y,z0,xD,yD,zD)
        I0_5_v2, I2_5_v2, I0_6_v2, I2_6_v2 = green_tensor_ANA_QE(omegac0,epsi1,x,y,z0,xD,yD,zD)
        if plot_numerical_sinQE == 1:
            I0_5_v3, I2_5_v3, I0_6_v3, I2_6_v3 = green_tensor_NUM(omegac0,epsi1,x,y,z0,xD,yD,zD)
        
        listI05.append(I0_5)
        listI25.append(I2_5)
        listI06.append(I0_6)
        listI26.append(I2_6)
        tot1 = I0_5 + I2_5 + I0_6 + I2_6
        total1_re.append(tot1.real)
        total1_im.append(tot1.imag)
    
    
        listI05_v2.append(I0_5_v2)
        listI25_v2.append(I2_5_v2)
        listI06_v2.append(I0_6_v2)
        listI26_v2.append(I2_6_v2)
        tot2 = I0_5_v2 + I2_5_v2 + I0_6_v2 + I2_6_v2
        total2_re.append(tot2.real)
        total2_im.append(tot2.imag)
    
        if plot_numerical_sinQE == 1:
            listI05_v3.append(I0_5_v3)
            listI25_v3.append(I2_5_v3)
            listI06_v3.append(I0_6_v3)
            listI26_v3.append(I2_6_v3)
            tot3 = I0_5_v3 + I2_5_v3 + I0_6_v3 + I2_6_v3
            total3_re.append(tot3.real)
            total3_im.append(tot3.imag)

#%%
print('Graficar el green tensor')
    
graph(title,labelx,r'$I^{(0)}_5$',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(list_x,listI05,'.',ms = ms,color = color_rgb1,label = 'numerical QE')
plt.plot(list_x,listI05_v2,'-',color = 'black',label = 'analytical QE')
# plt.yscale('log')
if plot_numerical_sinQE == 1:
    plt.plot(list_x,listI05_v3,'--',color = 'green',label = 'numerical')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
if save_graphs==1:
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'green_tensor_I05' + label1 + '.png', format='png')

graph(title,labelx,r'$I^{(2)}_5$',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(list_x,listI25,'.',ms = ms,color = color_rgb2,label = 'numerical QE')
plt.plot(list_x,listI25_v2,'-',color = 'black',label = 'analytical QE')
# plt.yscale('log')
if plot_numerical_sinQE == 1:
    plt.plot(list_x,listI25_v3,'--',color = 'green',label = 'numerical')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
if save_graphs==1:
    plt.tight_layout()
    plt.savefig( 'green_tensor_I25' + label1 + '.png', format='png')

graph(title,labelx,r'$I^{(0)}_6$',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(list_x,listI06,'.',ms = ms,color = color_rgb3,label = 'numerical QE')
plt.plot(list_x,listI06_v2,'-',color = 'black',label = 'analytical QE')
if plot_numerical_sinQE == 1:
    plt.plot(list_x,listI06_v3,'--',color = 'green',label = 'numerical')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
if save_graphs==1:
    plt.tight_layout()
    plt.savefig( 'green_tensor_I06' + label1 + '.png', format='png')

graph(title,labelx,r'$I^{(2)}_6$',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(list_x,listI26,'.',ms = ms,color = color_rgb4,label = 'numerical QE')
plt.plot(list_x,listI26_v2,'-',color = 'black',label = 'analytical QE')
# plt.yscale('log')
if plot_numerical_sinQE == 1:
    plt.plot(list_x,listI26_v3,'--',color = 'green',label = 'numerical')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
if save_graphs==1:
    plt.tight_layout()
    plt.savefig( 'green_tensor_I26' + label1 + '.png', format='png')

########################################################################################

graph(title,labelx,r'Re(Direct Green Tensor)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(list_x,total1_re,'.',ms = ms,color = 'black',label = 'numerical QE')
plt.plot(list_x,total2_re,'-',color = 'black',label = 'analytical QE')
# plt.yscale('log')
if plot_numerical_sinQE == 1:
    plt.plot(list_x,total3_re,'--',color = 'green',label = 'numerical')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
if save_graphs==1:
    plt.tight_layout()
    plt.savefig( 'Re_green_tensor_tot' + label1 + '.png', format='png')

graph(title,labelx,r'Im(Direct Green Tensor)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(list_x,total1_im,'.',ms = ms,color = 'black',label = 'numerical QE')
plt.plot(list_x,total2_im,'-',color = 'black',label = 'analytical QE')
# plt.yscale('log')
if plot_numerical_sinQE == 1:
    plt.plot(list_x,total3_im,'--',color = 'green',label = 'numerical')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
if save_graphs==1:
    plt.tight_layout()
    plt.savefig( 'Im_green_tensor_tot' + label1 + '.png', format='png')

#%%
