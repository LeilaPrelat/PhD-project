
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


energy0 = 2 # eV
d_nano = 1
hbar_omega_bulk = 9.17
hbar_gamma_in = 21*1e-3
 
title = r'$\epsilon_1$ = %i, $\epsilon_3$ = %i, d = %i nm, $\hbar\omega_{bulk}$ = %.2f eV, $\hbar\gamma_{in}$ = %i meV' %(epsi1,epsi3,d_nano,hbar_omega_bulk,hbar_gamma_in*1e3)
labelp = r'_E0%ieV_d%inm' %(energy0,d_nano)

N = 100

    # z0 = 0.06*1e3
labelx = r'$k_\parallel/\omega/c$'   
label1 = 'vs_alfa_parallel' + labelp
listx = np.linspace(0.01,200,N)


def function_num(energy0,u):
    omegac0 = energy0/aux 

    rp = rp_fresnel_num(omegac0,epsi1,epsi3,d_nano,u)
    # print(rta)
   

    return rp


def function_pole_approx(energy0,u):
    omegac0 = energy0/aux 

    rp  = rp_pole_aprox(omegac0,epsi1,epsi3,d_nano,u)
    
    return rp


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



listy_re_num = []
listy_re_pp = []

listy_im_num = []
listy_im_pp = []

for value in listx: 

    
    y_re_num = function_num(energy0,value)
    y_re_pole = function_pole_approx(energy0,value)


    listy_im_num.append(np.imag(y_re_num))
    listy_im_pp.append(np.imag(y_re_pole))


    listy_re_num.append(np.real(y_re_num))
    listy_re_pp.append(np.real(y_re_pole))



#%%
graph(title,labelx,r'Im{$r_p$}',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#    plt.plot(listx,listy_re_ana,'.',ms = ms,color = 'purple',label = 'analytical')
plt.plot(listx,listy_im_pp,'.-',ms = ms,color = 'purple',label = 'pole approx')
plt.plot(listx,listy_im_num,'.-',ms = ms,color = 'lightseagreen',label = 'Fresnel')
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0.05,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
#plt.yscale('log')
os.chdir(path_save)
plt.savefig( 'Im_rp_' + label1 + '.png', format='png')   

#%%

graph(title,labelx,r'Re{$r_p$}',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#    plt.plot(listx,listy_re_ana,'.',ms = ms,color = 'purple',label = 'analytical')
plt.plot(listx,listy_re_pp,'.-',ms = ms,color = 'purple',label = 'pole approx')
plt.plot(listx,listy_re_num,'.-',ms = ms,color = 'lightseagreen',label = 'Fresnel')
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0.05,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
#plt.yscale('log')
os.chdir(path_save)
plt.savefig( 'Re_rp_' + label1 + '.png', format='png')   