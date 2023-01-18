
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
import seaborn as sns
sns.set()

#%%
paper = 1

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/potential_field/potential_and_electric_field_with_dipole_moment_formula','')
path_save = path_basic + '/' + 'optimum_zp'
if not os.path.exists(path_save):
    print('Creating folder to save graphs')
    os.mkdir(path_save)


    
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

epsi1,epsi3 = 1,1
d_nano = 1

b = -0.01

d_nano = 1
d_micro = d_nano*1e-3

energy0_pol = 0.5 ## eV
omega0 = energy0_pol/hb 
#R = 10 # 10 nm en unidades de micrometros
kappa_factor_omega0 = 0.1
kappa_r_factor= 0.5

int_v = 3

dir_dip_moment1 = 'z'
dir_dip_moment2 = 'z'

#title1 = r'$\kappa$ = %.2f$\omega_0$, $\kappa_r$ = %.2f$\kappa$, $E_0$=%i meV' %(kappa_factor_omega0, kappa_r_factor, energy0_pol)     
#title2 = r'$\hbar\mu$ = %.2feV, $\hbar\gamma$ = %.4feV' %(hbmu,hbgama) 
#title3 = r'$z_p$=%inm, px=%i, py=%i, pz=%i' %(zp*1e3,px,py,pz)
title1 = r'$\kappa$ = %.2f$\omega_0$, $\kappa_r$ = %.2f$\kappa$, $\hbar\omega_o$=%.2f eV' %(kappa_factor_omega0, kappa_r_factor, energy0_pol) 
title4 = title1 + '\n' + r'b = %i nm, d = %i nm, hBN disks' %(b*1e3,d_nano)

labelp1 = r'_d%inm_p%s' %(d_nano,dir_dip_moment1)
labelp2 = r'_d%inm_p%s' %(d_nano,dir_dip_moment2)


labelx = r'$\hbar\omega$ [eV]'
labely = r'optimal $z_p$ [nm]'

N = 60

#%%

tamfig = (4.5,3.5)
tamlegend = 11
tamletra = 11
tamtitle = 10
tamnum = 9
labelpady = 0.5
labelpadx = 0.5
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

os.chdir(path_save)
tabla1 = np.loadtxt('zp_optimum_for_decay_rate' + labelp1 + '.txt', delimiter='\t', skiprows=1)
tabla1 = np.transpose(tabla1)
[listx1,listy1] = tabla1


tabla2 = np.loadtxt('zp_optimum_for_decay_rate' + labelp2 + '.txt', delimiter='\t', skiprows=1)
tabla2 = np.transpose(tabla2)
[listx2,listy2] = tabla2


#%%

graph(title4,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(listx1,listy1,'.-',ms = ms,color = 'purple',label = r'$\mathbf{p}$ $\parallel$ $\hat{z}$')
plt.plot(listx2,listy2,'.-',ms = ms,color = 'lightseagreen',label = r'$\mathbf{p}$ $\parallel$ $\hat{x}$')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=0.1,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
os.chdir(path_save)
#plt.yscale('log')
plt.savefig( 'zp_optimum_for_decay_rate_dif_p' + '.png', format='png')  


#%%

