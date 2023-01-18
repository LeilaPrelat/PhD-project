
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
    from decay_rate_film_resonance import EELS_film_ana_f , EELS_dir_ana_f,EELS_film_ana_f_sin_integrar
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
b = -0.01

d_nano = 20

int_v = 5



#title1 = r'$\kappa$ = %.2f$\omega_0$, $\kappa_r$ = %.2f$\kappa$, $E_0$=%i meV' %(kappa_factor_omega0, kappa_r_factor, energy0_pol)     
#title2 = r'$\hbar\mu$ = %.2feV, $\hbar\gamma$ = %.4feV' %(hbmu,hbgama) 
#title3 = r'$z_p$=%inm, px=%i, py=%i, pz=%i' %(zp*1e3,px,py,pz)
title4 = r'b = %i nm, v = c/%i' %(b*1e3,int_v)
labelp = r'_res_v%i_b%inm_d%inm' %(int_v,b*1e3,d_nano)



labelx = r'Plasmon energy, $\hbar\omega$ (eV)'
labely = r'Surface-dipole distance, $z_{\rm p}$ ($\mu$m)'


def function_imag_ana(energy0,zp_micro):
    omegac0 = energy0/aux 
    zp = zp_micro

    rta1 = EELS_film_ana_f_sin_integrar(omegac0,epsi1,epsi3,d_nano,int_v,b,zp)
    rta2 = EELS_dir_ana_f(omegac0,epsi1,epsi3,d_nano,int_v,b,zp)
    print( rta2 + rta1 , rta1)
    
    return rta1/(rta1+rta2)


title =  title4 


#%%
    
tamfig = [2.5, 2.25]
tamletra = 7
tamtitle  = 8
tamnum = 6
tamlegend = 6
labelpady = 2
labelpadx = 3
pad = 2.5
mk = 1
ms = 1
hp = 0.5
length_marker = 1.5
dpi = 500

def graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad):
    plt.figure(figsize=tamfig)
    plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
    plt.ylabel(labely,fontsize=tamletra,labelpad =labelpadx)
    plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in", pad = pad)
 #   plt.title(title,fontsize=int(tamtitle*0.9))

    return   
    
#%%
N = 50
listx = np.linspace(0.25,1.5,N)
listy = np.linspace(100,200,N)
listy = np.linspace(10,120,N)
if d_nano == 10:
    listy = np.linspace(80,200,N)
elif d_nano == 100:
    listy = np.linspace(50,800,N)
listy = np.array(np.linspace(0.001,100,N))
listy = np.array(np.linspace(0.00005,0.05,N))
X, Y = np.meshgrid(listx, listy, sparse=True)
f_opt = np.vectorize(function_imag_ana)
Z_opt = f_opt(X, Y)
#Z_opt = Z_opt/np.max(Z_opt)

#%%
ticks_z = [0,0.2,0.4,0.6,0.8,1]
import matplotlib.colors as mcolors
norm = mcolors.DivergingNorm(vmin=np.min(ticks_z), vmax = np.max(ticks_z),vcenter = 0.5)
limits = [np.min(listx) , np.max(listx), np.min(listy) , np.max(listy)]
graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)

im = plt.imshow(Z_opt, extent = limits, cmap='RdBu',  aspect='auto',origin = 'lower',norm=norm) 
cbar = plt.colorbar(im, fraction=0.046, pad=0.04, orientation = 'vertical')
cbar.ax.tick_params(labelsize = tamnum, width=0.5, length = 2,pad = 1.5)
cbar.ax.set_title(r'$\Gamma_{\rm SP}/\Gamma$',fontsize=tamlegend)
plt.tight_layout()
cbar.set_ticks(ticks_z)
cbar.set_ticklabels(ticks_z)
os.chdir(path_save)
plt.savefig( 'decay_rate_3D' + labelp + '.png', format='png',bbox_inches='tight',pad_inches = 0,dpi = dpi)   

#%%