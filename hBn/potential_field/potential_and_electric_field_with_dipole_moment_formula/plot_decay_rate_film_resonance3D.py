
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

#%%
paper = 1

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/potential_field/potential_and_electric_field_with_dipole_moment_formula','')
path_save = path_basic + '/' + 'decay_rate_film'
if not os.path.exists(path_save):
    print('Creating folder to save graphs')
    os.mkdir(path_save)

err = 'decay_rate_film3.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from decay_rate_film_resonance import EELS_film_ana_f , EELS_dir_ana_f
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

epsi1,epsi3 = 1,1
zp = 0.05
b = -0.01

d_nano = 10

int_v = 10

x1 = 0.09260651629072682 
x2 = 0.10112781954887218
x3 = 0.17030075187969923
x4 = 0.19937343358395992

#title1 = r'$\kappa$ = %.2f$\omega_0$, $\kappa_r$ = %.2f$\kappa$, $E_0$=%i meV' %(kappa_factor_omega0, kappa_r_factor, energy0_pol)     
#title2 = r'$\hbar\mu$ = %.2feV, $\hbar\gamma$ = %.4feV' %(hbmu,hbgama) 
#title3 = r'$z_p$=%inm, px=%i, py=%i, pz=%i' %(zp*1e3,px,py,pz)
title4 = r'b = %i nm, v = c/%i' %(b*1e3,int_v)
labelp = r'_res_v%i_b%inm_d%inm' %(int_v,b*1e3,d_nano)

N = 300

labelx = r'Plasmon energy, $\hbar\omega$ (eV)'
labely = r'Surface-dipole distance, $z_{\rm p}$ ($\mu$m)'


def function_imag_ana(energy0,zp_nano):
    omegac0 = energy0/aux 
    zp = zp_nano

    rta1 = EELS_film_ana_f(omegac0,epsi1,epsi3,d_nano,int_v,b,zp)
    rta2 = EELS_dir_ana_f(omegac0,epsi1,epsi3,d_nano,int_v,b,zp)
    
    return rta1/(rta1 + rta2)


title =  title4 


#%%
    
tamfig = [2.5, 2.25]
tamletra = 7
tamtitle  = 8
tamnum = 7
tamlegend = 6
labelpady = 2.5
labelpadx = 3
pad = 2
mk = 2
ms = 4
hp = 0.3
length_marker = 1
dpi = 500


def graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad):
    plt.figure(figsize=tamfig)
    plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
    plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
    plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in", pad = pad)
 #   plt.title(title,fontsize=int(tamtitle*0.9))

    return   
    
#%%

listx = np.linspace(x1*1.01,x4*0.99,N)
listy = np.linspace(50*1e-3,0.8,N)

X, Y = np.meshgrid(listx, listy, sparse=True)
f_opt = np.vectorize(function_imag_ana)
Z_opt = f_opt(X, Y)


#%%
ticks_z = [0,0.2,0.4,0.6,0.8,1]
ticks_z = [0.1,0.3,0.5,0.7,0.9]

limits = [np.min(listx) , np.max(listx), np.min(listy) , np.max(listy)]
graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)

im = plt.imshow(Z_opt, extent = limits, cmap='RdBu',  aspect='auto',origin = 'lower') 
cbar = plt.colorbar(im, fraction=0.046, pad=0.04, orientation = 'vertical')
plt.clim(np.min(ticks_z),np.max(ticks_z))
cbar.ax.tick_params(labelsize = tamnum, width=0.5, length = 2,pad = 1.5)
cbar.ax.set_title(r'$\Gamma_{\rm SP}/\Gamma$',fontsize=tamletra)
plt.tight_layout()
cbar.set_ticks(ticks_z)
cbar.set_ticklabels(ticks_z)
os.chdir(path_save)
plt.savefig( 'decay_rate_3D' + labelp + '.png', format='png', dpi=dpi)   
 

#%%
ticks_z = [0,0.2,0.4,0.6,0.8,1]
import matplotlib.colors as mcolors
limits = [np.min(listx) , np.max(listx), np.min(listy) , np.max(listy)]
graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)


norm = mcolors.DivergingNorm(vmin=np.min(ticks_z), vmax = np.max(ticks_z),vcenter = 0.5)
im = plt.imshow(Z_opt, extent = limits, cmap='RdBu',  aspect='auto',origin = 'lower',norm=norm) 
cbar = plt.colorbar(im, fraction=0.046, pad=0.04, orientation = 'vertical')
cbar.ax.tick_params(labelsize = tamnum, width=0.5, length = 2,pad = 1.5)
cbar.ax.set_title(r'$\Gamma_{\rm SP}/\Gamma$',fontsize=tamletra)
plt.tight_layout()
cbar.set_ticks(ticks_z)
cbar.set_ticklabels(ticks_z)
os.chdir(path_save)
plt.savefig( 'decay_rate_3D_v2' + labelp + '.png', format='png', dpi=dpi)   

#import matplotlib as mpl
#
#fig, ax = plt.subplots(figsize=(tamfig))
#
#cmap = mpl.cm.RdBu
#bounds = ticks_z
#norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
#plt.imshow(Z_opt, extent = limits, cmap='RdBu',  aspect='auto',origin = 'lower') 
#fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
#             cax=ax, orientation='vertical',
#             label=r'$\Gamma_{\rm SP}/\Gamma$')



#%%
#import matplotlib.colors as colors
#
#
#limits = [np.min(listx) , np.max(listx), np.min(listy) , np.max(listy)]
#graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#
#im = plt.imshow(Z_opt, extent = limits, cmap='RdBu',  aspect='auto',origin = 'lower') 
#
#cmap = plt.cm.RdBu
#vmin, vmax = np.min(ticks_z), np.max(ticks_z)
#pcm = plt.pcolormesh(X, Y, Z_opt,
#                     norm = colors.BoundaryNorm(bounds, cmap.N ) )


#pcm.ax.tick_params(labelsize = tamnum, width=0.5, length = 2,pad = 1.5)
#pcm.set_label(r'$\Gamma_{\rm SP}/\Gamma$',fontsize=tamletra,labelpad = 1.5)
#plt.tight_layout()
#pcm.set_ticks(ticks_z)
#pcm.set_ticklabels(ticks_z)































