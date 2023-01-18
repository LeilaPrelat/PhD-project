
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
import seaborn as sns
sns.set()

#%%
paper = 1

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/potential_field/potential_and_electric_field/potential_and_electric_field_with_dipole_moment_formula','')
path_save = path_basic + '/' + 'optimum_zp'
if not os.path.exists(path_save):
    print('Creating folder to save graphs')
    os.mkdir(path_save)

err = 'decay_rate_film3.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from decay_rate_film3_resonance import EELS_film_ana_f
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

print('Definir parametros del problema')

#v = c/int_v
#omega = 0.7*1e12

#v = c/int_v
#omega = 0.7*1e12

#v = c/int_v
#omega = 0.7*1e12

epsi1,epsi2 = 1,1
hbmu,hbgama = 0.3,0.0001

b = -0.01


 
#title1 = r'$\kappa$ = %.2f$\omega_0$, $\kappa_r$ = %.2f$\kappa$, $E_0$=%i meV' %(kappa_factor_omega0, kappa_r_factor, energy0_pol)     
#title2 = r'$\hbar\mu$ = %.2feV, $\hbar\gamma$ = %.4feV' %(hbmu,hbgama) 
#title3 = r'$z_p$=%inm, px=%i, py=%i, pz=%i' %(zp*1e3,px,py,pz)
title = r'b = %i nm, $\hbar\mu_c$ = %.2f eV' %(b*1e3,hbmu)
labelp = r'_b%inm' %(b*1e3)

N = 50


def function_imag_ana(energy0,v_sobre_c): ## devuelve el zp optimo en nanometros
    omegac0 = energy0*1e-3/aux 
    int_v = 1/v_sobre_c
    listx = np.linspace(300,4000,500)
    listy = []
    for zp_nano in listx:
        zp = zp_nano*1e-3
        rta = EELS_film_ana_f(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,b,zp)
        listy.append(rta)
#    print(energy0,v_sobre_c)
    
    peaks, _ = find_peaks(listy, height=0)
    maxi = listx[peaks]
    
    return maxi

labelx = r'$\hbar\omega$ [meV]'
labely = r'v/c'
labelz1 = r'$z_p$ [nm]'

listx = np.linspace(36,55,N)
listy = np.linspace(0.01,0.4,N)

#%%

X, Y = np.meshgrid(listx, listy, sparse=True)
f_opt = np.vectorize(function_imag_ana)
Z_opt = f_opt(X, Y)

#%%
    
tamfig = (4.5,3.5)
tamlegend = 10
tamletra = 11
tamtitle = 10
tamnum = 9
labelpady = 0
labelpadx = 0
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

limits = [np.min(listx) , np.max(listx), np.min(listy) , np.max(listy)]
graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)

im = plt.imshow(Z_opt, extent = limits, cmap=plt.cm.hot, aspect=50) 
cbar = plt.colorbar(im, extend='both', fraction=0.046, pad=0.04)
cbar.ax.tick_params(labelsize = tamnum)
if paper == 0:
    cbar.set_label(labelz1,fontsize=tamlegend,labelpad = 1)
plt.tight_layout()
os.chdir(path_save)
plt.savefig( 'zp_optimum_for_decay_rate_resonance3D' + labelp + '.png', format='png')   
 


#%%




