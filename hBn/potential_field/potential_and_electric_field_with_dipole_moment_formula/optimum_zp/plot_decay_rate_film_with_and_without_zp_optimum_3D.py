
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
path_basic2 = path_basic.replace('/optimum_zp' ,'')
path_constants =  path_basic.replace('/potential_field/potential_and_electric_field/potential_and_electric_field_with_dipole_moment_formula/optimum_zp','')
path_save = path_basic + '/' + 'decay_rate_with_optimum_zp'
if not os.path.exists(path_save):
    print('Creating folder to save graphs')
    os.mkdir(path_save)

err = 'decay_rate_film3.py no se encuentra en ' + path_basic2
try:
    sys.path.insert(1, path_basic2)
    from decay_rate_film3 import EELS_film_ana_f
except ModuleNotFoundError:
    print(err)


try:
    sys.path.insert(1, path_basic2)
    from decay_rate_film3_resonance import EELS_film_ana_f as EELS_film_ana_f_res
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
energy0_pol = 0.18
omega0 = energy0_pol/hb 
#R = 10 # 10 nm en unidades de micrometros
kappa_factor_omega0 = 0.1
kappa_r_factor= 0.5

int_v = 10
 
title1 = r'$\kappa$ = %.2f$\omega_0$, $\kappa_r$ = %.2f$\kappa$, $\hbar\omega_o$=%.2f eV' %(kappa_factor_omega0, kappa_r_factor, energy0_pol)     
#title2 = r'$\hbar\mu$ = %.2feV, $\hbar\gamma$ = %.4feV' %(hbmu,hbgama) 
#title3 = r'$z_p$=%inm, px=%i, py=%i, pz=%i' %(zp*1e3,px,py,pz)
title4 = r'b = %i nm, v = c/%i, hBn' %(b*1e3,int_v)
labelp = r'_E0%imeV' %(energy0_pol*1e3)
title = title1  + '\n'  + title4 

def function_ana(energy0,zp_nano): ## devuelve el zp optimo en nanometros
    omegac0 = energy0/aux
    zp = zp_nano*1e-3
    rta = EELS_film_ana_f(omegac0,epsi1,epsi3,d_nano,int_v,b,zp,omega0,kappa_factor_omega0,kappa_r_factor)    
    return rta

#%%
    
os.chdir(path_basic)

d_nano = 1
tabla1 = np.loadtxt('zp_optimum_for_decay_rate_hBn_d%inm.txt'%(d_nano), delimiter='\t', skiprows=1)
tabla1 = np.transpose(tabla1)
[listx1,listy1,listz1] = tabla1

d_nano = 10
tabla2 = np.loadtxt('zp_optimum_for_decay_rate_hBn_d%inm.txt'%(d_nano), delimiter='\t', skiprows=1)
tabla2 = np.transpose(tabla2)
[listx2,listy2,listz2] = tabla2


list_zp0_1 = np.linspace(np.min(listx1), np.max(listx1),np.len(listx1))
list_zp0_1 = np.linspace(np.min(listx1), np.max(listx1),np.len(listx2))

labelx = r'$\hbar\omega$ [eV]'
labely = r'$\Gamma(z^{opt}_p)/\Gamma(%i nm)$'%(zp0)


list_decay1 = []
for j in range(len(listx1)):
    energy = listx1[j]
    zp_opt  = listy1[j]  
    sol1 = function_ana(energy,zp_opt)
    sol2 = function_ana(energy,zp0)
    
    list_decay1.append(sol1/sol2)
    
    
list_decay2 = []    
for j in range(len(listx2)):
    energy = listx2[j]
    zp_opt  = listy2[j]  
    sol1 = function_ana(energy,zp_opt)
    sol2 = function_ana(energy,zp0)
    
    list_decay2.append(sol1/sol2)

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

graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(listx1,list_decay1,'.-',ms = ms,color = 'purple',label = 'd = %i nm' %(1))
plt.plot(listx2,list_decay2,'.-',ms = ms,color = 'lightseagreen',label = 'd = %i nm'%(10))
plt.legend(loc = 'best',markerscale=1.5,fontsize=tamlegend,frameon=0.1,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
plt.yscale('log')
os.chdir(path_save)
plt.savefig( 'decay_rate_with_and_without_zp_optimum' + labelp + '.png', format='png')  


#%%
###