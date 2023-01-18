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
#plot_numerical_sinQE = 1
plot_vs_R = 0 #graficar vs z
plot_vs_E = 1 #graficar vs omega en THz

sns.set()    

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/GreenTensor_checked','')
path_save = path_basic + '/' + 'green_tensor_direct_yy'

if save_graphs==1:    
    if not os.path.exists(path_save):
        print('Creating folder to save graphs')
        os.mkdir(path_save)

err = 'green_tensor_direct.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from green_tensor_direct_yy import green_tensor_NUM_QE, green_tensor_ANA_QE
except ModuleNotFoundError:
    print(err)
try:
    sys.path.insert(1, path_constants)
    from constants import constantes
except ModuleNotFoundError:
    print('constants.py no se encuentra en ' + path_constants)

pi,hb,c,alfac,mu1,mu2 = constantes()
aux = c*hb
aux2 = 1e12/c

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


N = 70

epsi1 = 1

cte_phi = 4
phi = np.pi/cte_phi
ze = -0.01
z = 0.05

title = r'$\varphi$ = $\pi$/%i, z = %i nm, $z_e$ = %i nm' %(cte_phi,z*1e3,ze*1e3)

########################
if plot_vs_E == 1:
#    omegaTHz0 = 0.7
#    omegaTHz0 = 1.511
    

    # omegaTHz = 1.663 - 0.152j
    
    R0 = 100
    listx = np.linspace(15,65,N)

    labelx = '$\hbar\omega$ [meV]'
    label1 = '_vs_E'
    title = title +  r', R = %i nm' %(R0)


elif plot_vs_R == 1:
#    omegaTHz0 = 0.7
#    omegaTHz0 = 1.511
    

    # omegaTHz = 1.663 - 0.152j
    
    E0 = 20
    listx = np.linspace(50,200,N)

    labelx = 'R [nm]'
    label1 = '_vs_R'
    title = title +  r', $\hbar\omega$ = %i meV' %(E0)
    
def green_tensor_ANA_QE_function(energy0,R_nano):
    
    omegac0 = energy0*1e-3/aux 
    R = R_nano*1e-3
    
    return green_tensor_ANA_QE(omegac0,epsi1,R,phi,z,ze)

def green_tensor_NUM_QE_function(energy0,R_nano):
    
    omegac0 = energy0*1e-3/aux 
    R = R_nano*1e-3
    
    return green_tensor_NUM_QE(omegac0,epsi1,R,phi,z,ze)

#%%


total1_re = []
total1_im = []


total2_re = []
total2_im = []

if plot_vs_R == 1:
    for value in listx: 
        
        tot2 = green_tensor_ANA_QE_function(E0,value)
        total2_re.append(tot2.real)
        total2_im.append(tot2.imag)
    

            
        tot1 = green_tensor_NUM_QE_function(E0,value)
   
        total1_re.append(tot1.real)
        total1_im.append(tot1.imag)
    

elif plot_vs_E == 1:
    for value in listx:
        
        tot2 = green_tensor_ANA_QE_function(value,R0)
        total2_re.append(tot2.real)
        total2_im.append(tot2.imag)
    

            
        tot1 = green_tensor_NUM_QE_function(value,R0)
   
        total1_re.append(tot1.real)
        total1_im.append(tot1.imag)


#%%
print('Graficar el green tensor')


########################################################################################
graph(title,labelx,r'Re($G_{dir}$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(listx,total2_re,'-',color = 'darkred',label = 'analytical')
# plt.yscale('log')
plt.plot(listx,total1_re,'.',ms = ms,color = 'green',label = 'numerical')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
if save_graphs==1:
    os.chdir(path_save)
    plt.tight_layout()
    plt.savefig( 'Re_green_tensor_tot' + label1 + '.png', format='png')

graph(title,labelx,r'Im($G_{dir}$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(listx,total2_im,'-',color = 'darkred',label = 'analytical')
# plt.yscale('log')
plt.plot(listx,total1_im,'.',ms = ms,color = 'green',label = 'numerical')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
if save_graphs==1:
    plt.tight_layout()
    plt.savefig( 'Im_green_tensor_tot' + label1 + '.png', format='png')

#%%
