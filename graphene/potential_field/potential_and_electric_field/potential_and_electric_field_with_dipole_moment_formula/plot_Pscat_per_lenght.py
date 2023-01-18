
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

plot_vs_E = 1
plot_vs_c = 0

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/potential_field/potential_and_electric_field/potential_and_electric_field_with_dipole_moment_formula','')
path_save = path_basic + '/' + 'Pscat_per_L'
if not os.path.exists(path_save):
    print('Creating folder to save graphs')
    os.mkdir(path_save)

err = 'EELS_film.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from Pscat_per_lenght import Pscat_inf_dipoles_per_lenght
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
a = 0.05

z = zp


energy0_pol = 43
omega0 = energy0_pol*1e-3/hb 
#R = 10 # 10 nm en unidades de micrometros
kappa_factor_omega0 = 0.1
kappa_r_factor= 0.5
 
title1 = r'$\kappa$ = %.2f$\omega_0$, $\kappa_r$ = %.2f$\kappa$, $E_0$=%i meV' %(kappa_factor_omega0, kappa_r_factor, energy0_pol)     
#title2 = r'$\hbar\mu$ = %.2feV, $\hbar\gamma$ = %.4feV' %(hbmu,hbgama) 
#title3 = r'$z_p$=%inm, px=%i, py=%i, pz=%i' %(zp*1e3,px,py,pz)
title4 = r'b = %i nm, a = %i nm, $z_p$ = %i nm' %(b*1e3,a*1e3,zp*1e3)
labelp = r'_E0%i' %(energy0_pol)

N = 400

if plot_vs_c == 1:
    E0 = 44 # meV
    # z0 = 0.06*1e3
    labelx = r'v/c'   
    title4 = title4 + ', ' + '$\hbar\omega$ = %i meV' %(E0)
    label1 = 'vs_v' + labelp
#    listx = np.linspace(0.0001,2,N)
    listx = np.linspace(0.01,0.6,N)

if plot_vs_E ==1:
    # z0 = 0.06*1e3
    int_v0 = 10
    
    labelx = r'$\hbar\omega$ [meV]'   
    title4 = title4 + ', ' + 'v = c/%i' %(int_v0)
    label1 = 'vs_E' + labelp
#    listx = np.linspace(0.0001,2,N)
    listx = np.linspace(15,65,N)


title = title1  + '\n'  + title4 

def function_imag_ana(energy0,int_v):
    omegac0 = energy0*1e-3/aux 

    rta = Pscat_inf_dipoles_per_lenght(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,a,b,zp,omega0,kappa_factor_omega0,kappa_r_factor)
    
    return np.real(rta)

#%%
    
tamfig = (4.5,3.5)
tamlegend = 10
tamletra = 11
tamtitle = 10
tamnum = 9
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


    

if plot_vs_E == 1: 

    listy_im_ana = []
    
    for value in listx: 

        y_im_ana = function_imag_ana(value,int_v0)        
        listy_im_ana.append(y_im_ana)


elif plot_vs_c == 1:       

    
    listy_im_ana = []
    
    for value in listx: 
        
        value2 = 1/value
#        print(value2)

        y_im_ana = function_imag_ana(E0,value2)        
        listy_im_ana.append(y_im_ana)
    
#%%



graph(title,labelx,'$P_{scat}$/L',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(listx,listy_im_ana,'.-',ms = ms,color = 'purple')
plt.legend(loc = 'best',markerscale=1.5,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
#if plot_vs_c == 1:
#    plt.yscale('log')
os.chdir(path_save)
plt.savefig( 'Pscat' + label1 + '.png', format='png')   


    