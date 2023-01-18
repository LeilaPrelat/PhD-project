
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: leila

graficar el campo externo directo con la convencion de z hacia abajo
en z = 0
graficar mapa de color x,y
"""

import numpy as np
import sys
import os 
import matplotlib.pyplot as plt

from scipy.signal import find_peaks

paper = 1

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_save = path_basic + '/' + 'rp_coefficient'
if not os.path.exists(path_save):
    print('Creating folder to save graphs')
    os.mkdir(path_save)

err = 'rp_coefficient.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from rp_coefficient import rp_fresnel_num, rp_pole_aprox
except ModuleNotFoundError:
    print(err)

try:
    sys.path.insert(1, path_basic)
    from constants import constantes
except ModuleNotFoundError:
    print('constants.py no se encuentra en ' + path_basic)

pi,hb,c,alfac,mu1,mu2 = constantes()
aux = c*hb
aux2 = 1e12/c

#%%

print('Definir parametros del problema')

epsi1,epsi3 = 1,1



d_nano = 10

 
title = r'$\epsilon_1$ = %i, $\epsilon_3$ = %i, d = %i nm' %(epsi1,epsi3,d_nano)


N = 150

    # z0 = 0.06*1e3
labelx = r'$k_\parallel$ [1/nm]'
labely =  r'$\hbar\omega$ [eV]'
label1 = 'rp_coeff_3D'

x1 = 0.09260651629072682 
x2 = 0.10112781954887218
x3 = 0.17030075187969923
x4 = 0.19937343358395992
listx = np.linspace(0.001,1,N)

rama = 1

if rama == 1:
    listy = np.linspace(0.09,0.1025,N) ## primera rama
elif rama == 2: 
    listy = np.linspace(0.169,0.201,N) ## segunda rama 

X, Y = np.meshgrid(listx, listy)


def function_num(u,energy0):
    omegac0 = energy0/aux 

    rp = rp_fresnel_num(omegac0,epsi1,epsi3,d_nano,u)
    # print(rta)
   

    return np.imag(rp)

def function_num_find_peaks(omegac0):

    list_y = []
    for u in listx:
        rp = rp_fresnel_num(omegac0,epsi1,epsi3,d_nano,u)
    
        list_y.append(np.imag(rp))
    
    
    
    peaks, _ = find_peaks(list_y, height=0)
    
    
    
    return listx[peaks]


def function_pole_approx(u,energy0):
    omegac0 = energy0/aux 

    rp  = rp_pole_aprox(omegac0,epsi1,epsi3,d_nano,u)
    
    return np.imag(rp)


#%%
    
tamfig = (4.5,3.5)
tamlegend = 12
tamletra = 12
tamtitle = 11
tamnum = 10
labelpady = -1.5
labelpadx = 2
pad = 0
mk = 2
ms = 4
hp = 0.3
length_marker = 1

def graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad):
    plt.figure(figsize=tamfig)
    plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
    plt.ylabel(labely,fontsize=tamletra,labelpad =labelpadx)
    plt.tick_params(labelsize = tamnum, pad = pad)
    plt.title(title,fontsize=int(tamtitle*0.9))

    return   

#%%

f_num = np.vectorize(function_num)
Z_num = f_num(X, Y)


f_pole_approx = np.vectorize(function_pole_approx)
Z_pole_approx = f_pole_approx(X, Y)


#%%

E0 = 0.18
omegac0 = E0/aux 

maxis = function_num_find_peaks(omegac0)

#%%

ejey_aux1 = np.linspace(np.min([np.min(Z_num),np.min(Z_num)]), np.max([np.max(Z_num),np.max(Z_num)]) , 10)
limits = [np.min(listx) , np.max(listx), np.min(listy) , np.max(listy)]
graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#for x in [x1,x2,x3,x4]:
#    plt.plot(ejey_aux1,x*np.ones(10),'--',color = 'grey' )
im = plt.imshow(Z_num, extent = limits, cmap=plt.cm.hot, aspect='auto', interpolation = 'none',origin = 'lower') 
#plt.plot(maxis,0.18*np.ones(len(maxis)),'o',color = 'green' )
cbar = plt.colorbar(im, extend='both', fraction=0.046, pad=0.04)
cbar.ax.tick_params(labelsize = tamnum)
if paper == 0:
    cbar.set_label('Im{$r_p$}',fontsize=tamlegend,labelpad = 1)
plt.tight_layout()
os.chdir(path_save)
plt.savefig( 'Fresnel_' +  label1 + '_rama%i.png' %(rama), format='png')   

#%%
ejey_aux1 = np.linspace(np.min([np.min(Z_pole_approx),np.min(Z_pole_approx)]), np.max([np.max(Z_pole_approx),np.max(Z_pole_approx)]) , 10)
limits = [np.min(listx) , np.max(listx), np.min(listy) , np.max(listy)]
graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#for x in [x1,x2,x3,x4]:
#    plt.plot(ejey_aux1,x*np.ones(10),'--',color = 'grey' )
im = plt.imshow(Z_pole_approx, extent = limits, cmap=plt.cm.hot, aspect='auto', interpolation = 'none',origin = 'lower')  
cbar = plt.colorbar(im, extend='both', fraction=0.046, pad=0.04)
cbar.ax.tick_params(labelsize = tamnum)
if paper == 0:
    cbar.set_label('Im{$r_p$}',fontsize=tamlegend,labelpad = 1)
plt.tight_layout()
os.chdir(path_save)
plt.savefig( 'PP_' + label1 + '_rama%i.png'%(rama), format='png')   


#%%
