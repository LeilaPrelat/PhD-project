
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
    from decay_rate_film_resonance import EELS_film_ana_f
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

epsi1,epsi3 = 1,1
d_nano = 1

b = -0.01


int_v = 10
 
#title2 = r'$\hbar\mu$ = %.2feV, $\hbar\gamma$ = %.4feV' %(hbmu,hbgama) 
#title3 = r'$z_p$=%inm, px=%i, py=%i, pz=%i' %(zp*1e3,px,py,pz)
title2 = r'b = %i nm, v = c/%i, d = %i nm' %(b*1e3,int_v,d_nano)
title = title2
labelp = r'_b%inm' %(b*1e3)
labelp_res = r'_b%inm' %(b*1e3)

#def function_ana(energy0_meV,zp_nano): ## devuelve el zp optimo en nanometros
#    omegac0 = energy0_meV*1e-3/aux
#    zp = zp_nano*1e-3
#    rta = EELS_film_ana_f(omegac0,epsi1,epsi3,d_nano,int_v,b,zp)    
#    return rta


def function_ana_res(energy0_meV,zp_nano): ## devuelve el zp optimo en nanometros
    omegac0 = energy0_meV*1e-3/aux
    zp = zp_nano*1e-3
    rta = EELS_film_ana_f(omegac0,epsi1,epsi3,d_nano,int_v,b,zp)    
    return rta

#%%
    
os.chdir(path_basic)
tabla_res = np.loadtxt('zp_optimum_for_decay_rate_resonance_d%inm_v%i.txt'%(d_nano,int_v), delimiter='\t', skiprows=1)
tabla_res = np.transpose(tabla_res)
[listx_res,listy_res,listz] = tabla_res



labelx = r'$\hbar\omega$ [meV]'
labely = r'$\Gamma(z^{opt}_p)/\Gamma(0 nm)$'

list_decay = []

list_rate_decay_rate = []
for j in range(len(listx_res)):
    energy = listx_res[j]
    zp_opt  = listy_res[j]  
    sol1 = function_ana_res(energy,zp_opt)
    sol2 = function_ana_res(energy,0)
    
    list_decay.append(sol1/sol2)
    
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
plt.plot(listx_res,list_decay,'.-',ms = ms,color = 'purple')
plt.tight_layout()
plt.yscale('log')
os.chdir(path_save)
plt.savefig( 'decay_rate_graphene_with_and_without_zp_optimum' + labelp + '.png', format='png')  


#%%
###
