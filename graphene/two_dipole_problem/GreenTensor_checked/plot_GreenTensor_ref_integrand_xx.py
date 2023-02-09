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


sns.set()

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('GreenTensor_checked/' + name_this_py,'')
path_save = path_basic + '/green_tensor_ref_integrand_xx'


err = 'green_tensor_ref_xx_graphene.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from integrand_green_tensor_ref_xx_graphene import green_tensor_NUM_p_fresnel_QE_integrand
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

ze = -10*1e-3
xe, ye = 0,0
x,y,z = 5000*1e-3,5000*1e-3,5000*1e-3
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
E0 = 30
R0 = 50

title1 = r'$\varphi$ = %i, z = %i nm, $z_e$ = %i nm, $z_p$ = %i nm' %(phi,z*1e3,ze*1e3,zp*1e3)
title2 = r'$\hbar\omega$ = %i meV, R = %i nm' %(E0,R0)
title = title1 + '\n' + title2

def function_numQE(alpha_parallel):
    omegac = E0/aux
    return green_tensor_NUM_p_fresnel_QE_integrand(omegac,epsi1,epsi2,hbmu,hbgama,R0,phi,z,ze,zp,alpha_parallel)    
    
########################



list_alpha_parallel = np.linspace(0.01,0.5,251)
list_x = list_alpha_parallel
labelx = r'$\alpha_\parallel = k_\parallel/k$'
label1 = '_vs_alpha_parallel'

    
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

total_num_re = []
total_num_im = []

for alpha_parallel in list_alpha_parallel: 
        
    num_z = function_numQE(alpha_parallel)

    total_num_re.append(np.real(num_z))
    total_num_im.append(np.imag(num_z))

#%%

print('Graficar el green tensor ' + label1)

############################################################################################


graph(title,labelx,r'Re(F($\alpha_\parallel$,k)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(list_x,total_num_re,'.-',ms = ms,color = 'black',label = 'numerical QE')
# plt.yscale('log')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
if save_graphs==1:
    os.chdir(path_save)
    plt.tight_layout()
    plt.savefig( 'green_tensor_ref_Re' + label1 + '.png', format='png')


graph(title,labelx,r'Im(F($\alpha_\parallel$,k)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(list_x,total_num_im,'.-',ms = ms,color = 'black',label = 'numerical QE')
# plt.yscale('log')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
if save_graphs==1:
    plt.tight_layout()
    plt.savefig( 'green_tensor_ref_Im' + label1 + '.png', format='png')

#%%
