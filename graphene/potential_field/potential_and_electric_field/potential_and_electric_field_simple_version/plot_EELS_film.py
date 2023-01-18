
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

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/potential_field/potential_and_electric_field/potential_and_electric_field_simple_version','')
path_save = path_basic + '/' + 'EELS_film3'
if not os.path.exists(path_save):
    print('Creating folder to save graphs')
    os.mkdir(path_save)

err = 'EELS_film.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from EELS_film3 import EELS_film_ana_f
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
cota = 75 #nanometros
epsi1,epsi2 = 1,1
hbmu,hbgama = 0.3,0.0001
zp = 0.05
b = -0.01

px = 1
py = 0
pz = 0

title2 = r'$\hbar\mu$ = %.2feV, $\hbar\gamma$ = %.4feV' %(hbmu,hbgama) 
title4 = r'$z_p$=%inm, px = %i, py = %i, pz = %i' %(zp*1e3,px,py,pz)
labelp = r'_px%i_py%i_pz%i' %(px,py,pz)

labelp = r'_px%i_py%i_pz%i' %(px,py,pz)

N = 400


    
labelx = r'$\hbar\omega$ [meV]'   
label1 = 'vs_E' + labelp
#    listx = np.linspace(0.0001,2,N)
listx = np.linspace(15,65,N)


title = title2 + '\n' + title4

def function_imag_ana(energy0):
    omegac0 = energy0*1e-3/aux 

    rta = EELS_film_ana_f(omegac0,epsi1,epsi2,hbmu,hbgama,px,py,pz,b,zp)
    
    return rta

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

listy_im_ana = []

for value in listx: 

    y_im_ana = function_imag_ana(value)        
    listy_im_ana.append(y_im_ana)

#%%



graph(title,labelx,'$\Gamma_{film}$',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(listx,listy_im_ana,'.-',ms = ms,color = 'purple')
plt.legend(loc = 'best',markerscale=1.5,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
if plot_vs_c == 1:
    plt.yscale('log')
os.chdir(path_save)
plt.savefig( 'EELS' + label1 + '.png', format='png')   


    