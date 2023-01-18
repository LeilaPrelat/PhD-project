
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
    from decay_rate_film3 import EELS_film_ana_f
except ModuleNotFoundError:
    print(err)
    
try:
    sys.path.insert(1, path_constants)
    from graphene_sigma import sigma_DL
except ModuleNotFoundError:
    print('graphene_sigma.py no se encuentra en ' + path_basic)

    
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

energy0_pol = 20
omega0 = energy0_pol*1e-3/hb 
print(omega0*1e-12)
#R = 10 # 10 nm en unidades de micrometros
kappa_factor_omega0 = 0.1
kappa_r_factor= 0.5

int_v = 10
 
title1 = r'$\kappa$ = %.2f$\omega_0$, $\kappa_r$ = %.2f$\kappa$, $\hbar\omega_o$=%i meV' %(kappa_factor_omega0, kappa_r_factor, energy0_pol)     
#title2 = r'$\hbar\mu$ = %.2feV, $\hbar\gamma$ = %.4feV' %(hbmu,hbgama) 
#title3 = r'$z_p$=%inm, px=%i, py=%i, pz=%i' %(zp*1e3,px,py,pz)
title = title1 + '\n' + r'b = %i nm, v = c/%i, $\hbar\mu_c$ = %.2f eV' %(b*1e3,int_v,hbmu)
labelp = r'_b%inm_E0%imeV' %(b*1e3,energy0_pol)

N = 100


def function_imag_ana(energy0): ## devuelve el zp optimo en nanometros
    omegac0 = energy0*1e-3/aux
    listx = np.linspace(50,900,300)
    listy = []
    for zp_nano in listx:
        zp = zp_nano*1e-3
        rta = EELS_film_ana_f(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,b,zp,omega0,kappa_factor_omega0,kappa_r_factor)
        listy.append(rta)
#    print(energy0,v_sobre_c)
    
    peaks, _ = find_peaks(listy, height=0)
    maxi = listx[peaks]
    print(energy0, maxi)
    if len(maxi ) > 1 :
        listy = np.array(listy)
        list_maxis_y = listy[peaks]
        
        maxi_ind = np.argmax(list_maxis_y)
        maxi =  listx[peaks[maxi_ind]]
        print(maxi)
    
    return float(maxi)


def lambda_p(energy0):
    
    E = energy0*1e-3
    
    omegac = E/aux
    
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    alfa_p = 1j*(epsi1 + epsi2)/(cond)
    kp = alfa_p*omegac

    return 2*np.pi/kp

labelx = r'$\hbar\omega$ [meV]'
labely = r'optimal $z_p$ [nm]'

listx = np.linspace(20,65,N)
listy = []
list_lambda_p = []
for x in listx:
    listy.append(function_imag_ana(x))
    list_lambda_p.append(np.real(lambda_p(x))*1e3)
    
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
plt.plot(listx,listy,'.',ms = ms,color = 'purple',label = r'opt $z_p$')
plt.plot(listx,list_lambda_p,'.',ms = ms,color = 'lightseagreen',label = r'$\lambda_p$')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=0.1,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
os.chdir(path_save)
plt.savefig( 'zp_optimum_for_decay_rate' + labelp + '.png', format='png')  

graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(listx,listy,'.-',ms = ms,color = 'purple',label = r'opt $z_p$')
#plt.plot(listx,list_lambda_p,'.',ms = ms,color = 'lightseagreen',label = r'$\lambda_p$')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=0.1,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
os.chdir(path_save)
plt.savefig( 'zp_optimum_for_decay_rate_v2' + labelp + '.png', format='png')  


#%%
### change some numbers 

tabla = np.array([listx,listy])
tabla = np.transpose(tabla)
info_title = title1 + ', ' + r'b = %i nm, v = c/%i, $\hbar\mu_c$ = %.2f eV' %(b*1e3,int_v,hbmu)
info = 'zp_optimum_for_decay_rate'  
header1 = 'E [meV]     zp [nm]' + info + ', '  + info_title + ', ' + name_this_py
np.savetxt('zp_optimum_for_decay_rate' + labelp + '.txt', tabla, fmt='%1.11e', delimiter='\t', header = header1)

#%%

tabla = np.loadtxt('zp_optimum_for_decay_rate.txt', delimiter='\t', skiprows=1)
tabla = np.transpose(tabla)
[listx,listy] = tabla

#%%

graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(listx,listy,'.-',ms = ms,color = 'purple',label = r'opt $z_p$')
plt.plot(listx,list_lambda_p,'.-',ms = ms,color = 'lightseagreen',label = r'$\lambda_p$')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=0.1,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
os.chdir(path_save)
plt.savefig( 'zp_optimum_for_decay_rate' + labelp + '.png', format='png')  


