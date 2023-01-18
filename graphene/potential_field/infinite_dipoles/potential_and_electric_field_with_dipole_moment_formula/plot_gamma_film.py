
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

plot_vs_v = 0
plot_vs_a = 1
plot_vs_zp = 0

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/potential_field/infinite_dipoles/potential_and_electric_field_with_dipole_moment_formula','')
path_save = path_basic + '/' + 'Gamma_film'
if not os.path.exists(path_save):
    print('Creating folder to save graphs')
    os.mkdir(path_save)

try:
    sys.path.insert(1, path_constants)
    from gamma_film import Gamma_film
except ModuleNotFoundError:
    print('dipole_moment.py no se encuentra en ' + path_constants)


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
cota = 75 #nanometros
epsi1,epsi2 = 1,1
hbmu,hbgama = 0.3,0.0001



energy0_pol = 43
omega0 = energy0_pol*1e-3/hb 
#R = 10 # 10 nm en unidades de micrometros
kappa_factor_omega0 = 0.1
kappa_r_factor= 0.5

#Ndipoles = 50
#list_xD = np.linspace(-0.01,0.01,Ndipoles)
b = -0.01
x = y = z = 0

title1 = r'$\kappa$ = %.2f$\omega_0$, $\kappa_r$ = %.2f$\kappa$, $E_0$=%i meV' %(kappa_factor_omega0, kappa_r_factor, energy0_pol)    
title2 = r'b = %i nm' %(b*1e3) 



N = 150
if plot_vs_v == 1: 
    
    a0 = 0.05*1e3
    zp0 =  0.05*1e3
    
    labelx = 'v/c'  
    label1 = 'vs_v' 
    title = title1 + '\n' + title2 + r', a = %i nm, $z_p$ = % i nm' %(a0,zp0)

    
    listx = np.linspace(0.001,0.2,N) ##v/c


elif plot_vs_a == 1:
    
    int_v0 = 10
    zp0 =  0.05*1e3    
    
    labelx = 'a [nm]' 
    label1 = 'vs_a' 
    title = title1 + '\n' + title2 + r', v = c/%i, $z_p$ = % i nm' %(int_v0,zp0)


    listx = np.linspace(1, 1000,N)

elif plot_vs_zp == 1:

    int_v0 = 10
    a0 =  0.05*1e3 
    
    labelx = '$z_p$ [nm]'   
    label1 = 'vs_zp' 
    title = title1 + '\n' + title2 + ', v = c/%i, a = % i nm' %(int_v0,a0)

    list_n = [1,2,3]
    listx = np.linspace(4, 16,N)

#%%
    
    
def Gamma_film_f(int_v,zp_nm,a_nm,n):

    
    def omega_n(int_v,a_nm,n):
    
        a = a_nm*1e-3
        
        aux = alfac*int_v*hbmu/(hb)
        rta = aux + 0.5*np.sqrt((2*aux)**2 + 16*alfac*c*hbmu*n*np.pi/(a*hb))
        
        return rta


    omegac = omega_n(int_v,a_nm,n)/c
    
    zp = zp_nm*1e-3
    a = a_nm*1e-3
    
    gamma = Gamma_film(omegac,epsi1,epsi2,hbmu,hbgama,int_v,zp,a,b,x,y,z,omega0,kappa_factor_omega0,kappa_r_factor,n)
    
    return np.imag(gamma)



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
    
labely = r'$\Gamma_{film}$ [eV$^{-1}$]'

listy_im_ana_tot = []


if plot_vs_v == 1: 
    for n in list_n:
        listy_im_ana = []       
        for value in listx:    
            y_im_ana = Gamma_film_f(1/value,zp0,a0,n)        
            listy_im_ana.append(y_im_ana)           
        listy_im_ana_tot.append(listy_im_ana)
        

elif plot_vs_a == 1: 
    for n in list_n:
        listy_im_ana = []       
        for value in listx:    
            y_im_ana = Gamma_film_f(int_v0,zp0,value,n)      
            listy_im_ana.append(y_im_ana)           
        listy_im_ana_tot.append(listy_im_ana)
        
elif plot_vs_zp == 1: 
    for n in list_n:
        listy_im_ana = []       
        for value in listx:    
            y_im_ana = Gamma_film_f(int_v0,value,a0,n)     
            listy_im_ana.append(y_im_ana)           
        listy_im_ana_tot.append(listy_im_ana)
        
    
#%%       
colors = ['darkred','steelblue','yellowgreen','coral']  
graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
k = 0
for modo in list_n:
    plt.plot(listx,np.array(listy_im_ana_tot[k])/hb,'.-',color = colors[k],ms = ms, label = 'n = %i' %(modo))
    k = k + 1

plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
if plot_vs_v == 1:
    plt.yscale('log')
os.chdir(path_save)
plt.savefig( 'EELS_film' + label1 + '.png', format='png')   
   
#%%
