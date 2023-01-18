
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
    from rp_coefficient import rp_fresnel_num, rp_pole_aprox, rp_pole_aprox_v2
except ModuleNotFoundError:
    print(err)
    
    

try:
    sys.path.insert(1, path_basic)
    from Silica_epsilon import epsilon_Silica
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




energy0 = 0.2 # eV  np.linspace(0.087, 0.8, N) ### for higher freq de 
#energy0 = 0.2

d_nano = 0.4
zp_nano = 10
 
title = r'$\hbar\omega$ = %.3f eV, $d$ = %.2f nm, $z_0$ = %i nm' %(energy0,d_nano,zp_nano)
labelp = r'd%.1fnm_E0%imeV' %(d_nano,energy0*1e3)

N = 100

    # z0 = 0.06*1e3
labelx = r'$k_\parallel$ [nm]'   
label1 = 'vs_alfa_parallel_' + labelp
listx = np.linspace(0.01,3000,N)


def function_num(energy0,u):
    omegac0 = energy0/aux 

    rp = rp_fresnel_num(omegac0,epsilon_Silica,d_nano,zp_nano,u)
    # print(rta)
   

    return rp


def function_pole_approx(energy0,u):
    omegac0 = energy0/aux 

    rp  = rp_pole_aprox(omegac0,epsilon_Silica,d_nano,zp_nano,u)
    
    return rp



def function_pole_approx2(energy0,u):
    omegac0 = energy0/aux 

    rp  = rp_pole_aprox_v2(omegac0,epsilon_Silica,d_nano,zp_nano,u)
    
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
listy_re_pp2 = []

listy_im_num = []
listy_im_pp = []
listy_im_pp2 = []

for value in listx: 
    
    print(value)

    
    y_re_num = function_num(energy0,value)
    y_re_pole = function_pole_approx(energy0,value)
    y_re_pole2 = function_pole_approx2(energy0,value)

    listy_im_num.append(np.imag(y_re_num))
    listy_im_pp.append(np.imag(y_re_pole))
    listy_im_pp2.append(np.imag(y_re_pole2))

    listy_re_num.append(np.real(y_re_num))
    listy_re_pp.append(np.real(y_re_pole))
    listy_re_pp2.append(np.real(y_re_pole2))


#%%
graph(title,labelx,r'Im{$r_p$}',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(listx,listy_im_pp2,'.-',ms = ms,color = 'purple',label =  r'$r_{\rm p} = k_{\rm p}/(k_\parallel - k_{\rm p})$')
plt.plot(listx,listy_im_pp,'.',ms = ms+1,color = 'blue',label = r'$r_{\rm p} = k_\parallel/(k_\parallel - k_{\rm p})$')
plt.plot(listx,listy_im_num,'.-',ms = ms,color = 'lightseagreen',label = 'Fresnel')
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0.05,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
#plt.yscale('log')
os.chdir(path_save)
plt.savefig( 'Im_rp_' + label1 + '.png', format='png')   


graph(title,labelx,r'Re{$r_p$}',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(listx,listy_re_pp2,'.-',ms = ms,color = 'purple',label =  r'$r_{\rm p} = k_{\rm p}/(k_\parallel - k_{\rm p})$')
plt.plot(listx,listy_re_pp,'.',ms = ms+1,color = 'blue',label =  r'$r_{\rm p} = k_\parallel/(k_\parallel - k_{\rm p})$')
plt.plot(listx,listy_re_num,'.-',ms = ms,color = 'lightseagreen',label = 'Fresnel')
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0.05,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
#plt.yscale('log')
os.chdir(path_save)
plt.savefig( 'Re_rp_' + label1 + '.png', format='png')   
