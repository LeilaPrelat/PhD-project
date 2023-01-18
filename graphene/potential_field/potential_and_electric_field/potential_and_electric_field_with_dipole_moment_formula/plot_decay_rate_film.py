
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
from scipy import special
import seaborn as sns
sns.set()

#%%
plot_vs_E = 0
plot_vs_c = 0
plot_vs_zp = 1

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/potential_field/potential_and_electric_field/potential_and_electric_field_with_dipole_moment_formula','')
path_save = path_basic + '/' + 'decay_rate_film'
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
cota = 75 #nanometros
epsi1,epsi2 = 1,1
hbmu,hbgama = 0.3,0.0001
zp = 0.05
b = -0.01


energy0_pol = 20
omega0 = energy0_pol*1e-3/hb 
#R = 10 # 10 nm en unidades de micrometros
kappa_factor_omega0 = 0.1
kappa_r_factor= 0.5
 
title1 = r'$\kappa$ = %.2f$\omega_0$, $\kappa_r$ = %.2f$\kappa$, $\hbar\omega_o$=%i meV' %(kappa_factor_omega0, kappa_r_factor, energy0_pol)     
#title2 = r'$\hbar\mu$ = %.2feV, $\hbar\gamma$ = %.4feV' %(hbmu,hbgama) 
#title3 = r'$z_p$=%inm, px=%i, py=%i, pz=%i' %(zp*1e3,px,py,pz)
title4 = r'b = %i nm' %(b*1e3)
labelp = r'_E0%i' %(energy0_pol)

N = 200





def function_imag_ana(energy0,int_v,zp_nano):
    omegac0 = energy0*1e-3/aux 
    zp = zp_nano*1e-3

    rta = EELS_film_ana_f(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,b,zp,omega0,kappa_factor_omega0,kappa_r_factor)
    
    return rta


def minimum_function(energy0,int_v):
    
    omegac = energy0*1e-3/aux 
    arg = np.abs(b)*omegac*int_v
    K1 = special.kn(1,arg)
    K0 = special.kn(0,arg)

    E = omegac*aux
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    Rp = 2*epsi1/(epsi1 + epsi2)
    alfa_p = 1j*(epsi1 + epsi2)/(cond)
    kp = alfa_p*omegac


    omegav = omegac*int_v
    omegav_2 = omegav**2
    
    kp_2 = kp**2
    Rp_2 = Rp**2
    
    expo_2_menos = np.exp(-2*kp*np.abs(b))
    expo_2_mas = np.exp(2*kp*np.abs(b))
    
    num1 = -4*kp*Rp*omegav*K1 

    num_aux = -9*expo_2_menos*kp_2*(np.abs(K0)**2 + 6*np.abs(K1)**2)/(kp_2 - omegav_2) + 16*np.abs(K1)**2                 
    num2 = np.sqrt(kp_2*Rp_2*omegav_2*num_aux)

    kp_4 = kp**4

    tot = -(num1 + num2)*expo_2_mas*(kp - omegav)*(kp + omegav)/(18*kp_4*np.pi*Rp_2)


    return -np.log(tot)/(2*kp)


if plot_vs_c == 1 :
    E0 = 44 # meV
    # z0 = 0.06*1e3
    zp0 = 0.05*1e3 
    
    labelx = r'v/c'   
    title4 = title4 + ', ' + r'$\hbar\omega$ = %i meV, $z_p$ = %i nm' %(E0,zp0)
    label1 = 'vs_v' + labelp
#    listx = np.linspace(0.0001,2,N)
    listx = np.linspace(0.005,0.12,N)

if plot_vs_E ==1 :
    # z0 = 0.06*1e3
    int_v0 = 10
    zp0 = 0.05*1e3 
    
    labelx = r'$\hbar\omega$ [meV]'   
    title4 = title4 + ', ' + r'v = c/%i, $z_p$ = %i nm' %(int_v0,zp0)
    label1 = 'vs_E' + labelp
#    listx = np.linspace(0.0001,2,N)
    listx = np.linspace(15,65,N)

if plot_vs_zp == 1 : 
    E0 = 40 # meV
    int_v0 = 10

    labelx = r'$z_p$ [nm]'   
    title4 = title4 + ', ' + r'v = c/%i, $\hbar\omega$ = %i meV' %(int_v0,E0)
    label1 = 'vs_zp' + labelp
#    listx = np.linspace(0.0001,2,N)
    listx = np.linspace(10,6000,N)
    
    print(minimum_function(E0,int_v0)*1e3)
    print(np.abs(minimum_function(E0,int_v0))*2*1e3)



title = title1  + '\n'  + title4 


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

        y_im_ana = function_imag_ana(value,int_v0,zp0)        
        listy_im_ana.append(y_im_ana)


elif plot_vs_c == 1:       

    
    listy_im_ana = []
    
    for value in listx: 
        
        value2 = 1/value
#        print(value2)

        y_im_ana = function_imag_ana(E0,value2,zp0)        
        listy_im_ana.append(y_im_ana)
        
elif plot_vs_zp == 1:
    
    listy_im_ana = []
    
    for value in listx: 

        y_im_ana = function_imag_ana(E0,int_v0,value)        
        listy_im_ana.append(y_im_ana)
    
#%%

graph(title,labelx,'$\Gamma_{film}$ [seg/nm$^2$]',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(listx,listy_im_ana,'.-',ms = ms,color = 'purple')
plt.legend(loc = 'best',markerscale=1.5,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
#if plot_vs_c == 1:
#    plt.yscale('log')
os.chdir(path_save)
plt.savefig( 'EELS' + label1 + '.png', format='png')   

