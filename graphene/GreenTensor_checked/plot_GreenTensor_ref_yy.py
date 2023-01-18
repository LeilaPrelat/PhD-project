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

plot_vs_R = 0 #graficar vs z
plot_vs_E = 1 #graficar vs omega en THz

sns.set()

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('GreenTensor_checked/' + name_this_py,'')
path_save = path_basic + '/green_tensor_ref_yy'


err = 'green_tensor_ref_xx_graphene.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from green_tensor_ref_yy_graphene import green_tensor_ref_fresnel, green_tensor_ref_pole_aprox, green_tensor_PP2
except ModuleNotFoundError:
    print(err)
#    path_save = path_basic + '/' + 'green_tensor_ref_xx_graphene'  


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



#omega = 0.7*1e12
epsi1 = 1
epsi2 = 1

zp = 50*1e-3 #micrones # posicion del plano (z apunta hacia abajo)
phi = 0 
z = zp
ze = -10*1e-3

hbmu, hbgama = 0.3, 0.0001 #potencial quimico in ev, freq de collision in ev

#title1 = r'$\epsilon_1$ = %i, $\epsilon_2$ = %i, $\varphi$ = %i' %(epsi1,epsi2)
#title3 = r'$\hbar \mu$ = %.2feV, $\hbar \gamma$ = %.4feV' %(hbmu,hbgama)

title = r'$\varphi$ = %i, z = %i nm, $z_e$ = %i nm, $z_p$ = %i nm' %(phi,z*1e3,ze*1e3,zp*1e3)

def function_fresnel(energy0,R_nano):
    
    omegac0 = energy0*1e-3/aux 
    R = R_nano*1e-3
    return green_tensor_ref_fresnel(omegac0,epsi1,epsi2,hbmu,hbgama,R,phi,z,ze,zp)

def function_ana(energy0,R_nano):
    
    omegac0 = energy0*1e-3/aux 
    R = R_nano*1e-3
    return green_tensor_PP2(omegac0,epsi1,epsi2,hbmu,hbgama,R,phi,z,ze,zp)


def function_pole_aprox(energy0,R_nano):
    
    omegac0 = energy0*1e-3/aux 
    R = R_nano*1e-3
    return green_tensor_ref_pole_aprox(omegac0,epsi1,epsi2,hbmu,hbgama,R,phi,z,ze,zp)
        
#%%
N = 70
   
########################
if plot_vs_R == 1:
    E0 = 43
    listx = np.linspace(200,4000,N)

    labelx = 'R [nm]'
    label1 = '_vs_R'
    title = title +  r', $\hbar\omega$ = %i meV' %(E0)
    
########################
if plot_vs_E == 1:
#    omegaTHz0 = 0.7
#    omegaTHz0 = 1.511
    

    # omegaTHz = 1.663 - 0.152j
    
    R0 = 1000
    listx = np.linspace(15,65,N)

    labelx = '$\hbar\omega$ [meV]'
    label1 = '_vs_E'
    title = title +  r', R = %i nm' %(R0)
    
#%%

tamfig = (4.5,3.5)
tamlegend = 12
tamletra = 12
tamtitle = 11
tamnum = 10
labelpady = 2
labelpadx = 2
pad = 0
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

total_pole_aprox_re = []
total_pole_aprox_im = []
total_pole_aprox_abs = []


totalPP2_re = []
totalPP2_im = []
totalPP2_abs = []

total_num_re = []
total_num_im = []
total_num_abs = []

if plot_vs_R == 1:
    
    for value in listx: 
        
  
        num = function_fresnel(E0,value)
        ana = function_ana(E0,value)
        pole_aprox = function_pole_aprox(E0,value)
        
        total_num_re.append(np.real(num))
        total_num_im.append(np.imag(num))
        total_num_abs.append(np.abs(num))
        
        
        total_pole_aprox_re.append(np.real(pole_aprox))
        total_pole_aprox_im.append(np.imag(pole_aprox))
        total_pole_aprox_abs.append(np.abs(pole_aprox))


        totalPP2_re.append(np.real(ana))
        totalPP2_im.append(np.imag(ana))
        totalPP2_abs.append(np.abs(ana))


elif plot_vs_E == 1 :


    for value in listx: 
        
  
        num = function_fresnel(value,R0)
        ana = function_ana(value,R0)
        pole_aprox = function_pole_aprox(value,R0)
        
        total_num_re.append(np.real(num))
        total_num_im.append(np.imag(num))
        total_num_abs.append(np.abs(num))
        
        
        total_pole_aprox_re.append(np.real(pole_aprox))
        total_pole_aprox_im.append(np.imag(pole_aprox))
        total_pole_aprox_abs.append(np.abs(pole_aprox))


        totalPP2_re.append(np.real(ana))
        totalPP2_im.append(np.imag(ana))
        totalPP2_abs.append(np.abs(ana))

#%%

print('Graficar el green tensor ' + label1)

############################################################################################

#if plot_vs_E == 1:
#    from scipy import signal
#    sos = signal.butter(4, 50, 'hp', fs=1000, output='sos')
#    
#    filtered_re = signal.sosfilt(sos, total_num_re)
#    filtered_im = signal.sosfilt(sos, total_num_im)

graph(title,labelx,r'Re($G_{ref}$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(listx,total_num_re,'.-',ms = ms,color = 'lightseagreen',label = 'full numerical')
plt.plot(listx,total_pole_aprox_re,'.-',ms = ms,color = 'darkred',label = 'numerical PP')
plt.plot(listx,totalPP2_re,'.-',ms = ms,color = 'purple',label =  'analytical PP')
# plt.yscale('log')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
if save_graphs==1:
    os.chdir(path_save)
    plt.tight_layout()
    plt.savefig( 'green_tensor_ref_Re' + label1 + '.png', format='png')



graph(title,labelx,r'Im($G_{ref}$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(listx,total_num_im,'.-',ms = ms,color = 'lightseagreen',label = 'full numerical')
plt.plot(listx,total_pole_aprox_im,'.-',ms = ms,color = 'darkred',label =  'numerical PP')
plt.plot(listx,totalPP2_im,'.-',ms = ms,color = 'purple',label =  'analytical PP')
# plt.yscale('log')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
if save_graphs==1:
    plt.tight_layout()
    plt.savefig( 'green_tensor_ref_Im' + label1 + '.png', format='png')
    
    

graph(title,labelx,r'|$G_{ref}$|',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(listx,total_num_abs,'.-',ms = ms,color = 'lightseagreen',label = 'full numerical')
plt.plot(listx,total_pole_aprox_abs,'.-',ms = ms,color = 'darkred',label =  'numerical PP')
plt.plot(listx,totalPP2_abs,'.-',ms = ms,color = 'purple',label =  'analytical PP')
# plt.yscale('log')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
if save_graphs==1:
    plt.tight_layout()
    plt.savefig( 'green_tensor_ref_abs' + label1 + '.png', format='png')
    

#%%
