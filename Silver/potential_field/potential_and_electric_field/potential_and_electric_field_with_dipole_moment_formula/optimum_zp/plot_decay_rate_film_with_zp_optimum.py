
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

#v = c/int_v
#omega = 0.7*1e12

#v = c/int_v
#omega = 0.7*1e12

#v = c/int_v
#omega = 0.7*1e12

epsi1,epsi2 = 1,1
hbmu,hbgama = 0.3,0.0001

b = -0.01

energy0_pol = 43
omega0 = energy0_pol*1e-3/hb 
#print(omega0*1e-12)
#R = 10 # 10 nm en unidades de micrometros
kappa_factor_omega0 = 0.1
kappa_r_factor= 0.5

int_v = 10
 
title1 = r'$\kappa$ = %.2f$\omega_0$, $\kappa_r$ = %.2f$\kappa$, $\hbar\omega_o$=%i meV' %(kappa_factor_omega0, kappa_r_factor, energy0_pol)     
#title2 = r'$\hbar\mu$ = %.2feV, $\hbar\gamma$ = %.4feV' %(hbmu,hbgama) 
#title3 = r'$z_p$=%inm, px=%i, py=%i, pz=%i' %(zp*1e3,px,py,pz)
title2 = r'b = %i nm, v = c/%i, $\hbar\mu_c$ = %.2f eV' %(b*1e3,int_v,hbmu)
title = title1 + '\n' + title2
labelp = r'_b%inm_E0%imeV' %(b*1e3,energy0_pol)
labelp_res = r'_b%inm' %(b*1e3)

def function_ana(energy0_meV,zp_nano): ## devuelve el zp optimo en nanometros
    omegac0 = energy0_meV*1e-3/aux
    zp = zp_nano*1e-3
    rta = EELS_film_ana_f(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,b,zp,omega0,kappa_factor_omega0,kappa_r_factor)    
    return rta


def function_ana_res(energy0_meV,zp_nano): ## devuelve el zp optimo en nanometros
    omegac0 = energy0_meV*1e-3/aux
    zp = zp_nano*1e-3
    rta = EELS_film_ana_f_res(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,b,zp)    
    return rta

#%%
    
os.chdir(path_basic)
tabla_res = np.loadtxt('zp_optimum_for_decay_rate_graphene_resonance_b%inm.txt'%(b*1e3), delimiter='\t', skiprows=1)
tabla_res = np.transpose(tabla_res)
[listx_res,listy_res] = tabla_res


tabla = np.loadtxt('zp_optimum_for_decay_rate_graphene_b%inm_E0%imeV.txt'%(b*1e3,energy0_pol), delimiter='\t', skiprows=1)
tabla = np.transpose(tabla)
[listx,listy] = tabla

labelx = r'$\hbar\omega$ [meV]'
labely = r'$\Gamma_{decay}$'

list_decay = []
list_decay_res = []
list_decay_res_v2 = []
for j in range(len(listx)):
    energy = listx[j]
    zp_opt  = listy[j]    
    list_decay.append(function_ana(energy,zp_opt))

for j in range(len(listx_res)):
    energy = listx_res[j]
    zp_opt  = listy_res[j]    
    list_decay_res_v2.append(function_ana(energy,zp_opt))


for j in range(len(listx_res)):
    energy = listx_res[j]
    zp_opt  = listy_res[j]  
    list_decay_res.append(function_ana_res(energy,zp_opt))

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

graph(title2,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(listx_res,list_decay_res,'.',ms = ms,color = 'lightseagreen',label = r'$\Gamma_{res} (z^{res}_p)$')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=0.1,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
os.chdir(path_save)
plt.savefig( 'decay_rate_res_with_zp_optimum' + labelp_res + '.png', format='png')  


graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(listx,list_decay,'.',ms = ms,color = 'purple',label = r'$\Gamma (z_p)$')
plt.plot(listx_res,list_decay_res_v2,'.',ms = ms,color = 'lightseagreen',label = r'$\Gamma (z^{res}_p)$')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=0.1,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
plt.yscale('log')
os.chdir(path_save)
plt.savefig( 'decay_rate_graphene_with_zp_optimum' + labelp + '.png', format='png')  


#%%
###
