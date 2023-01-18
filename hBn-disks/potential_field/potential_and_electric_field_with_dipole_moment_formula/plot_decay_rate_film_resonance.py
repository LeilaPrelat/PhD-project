
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

plot_vs_E = 0
plot_vs_c = 0
plot_vs_zp = 1

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/potential_field/potential_and_electric_field_with_dipole_moment_formula','')
path_save = path_basic + '/' + 'decay_rate_film'
if not os.path.exists(path_save):
    print('Creating folder to save graphs')
    os.mkdir(path_save)

err = 'decay_rate_film3.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from decay_rate_film_resonance import EELS_film_ana_f, EELS_film_num_f, EELS_film_pole_aprox_f
except ModuleNotFoundError:
    print(err)


try:
    sys.path.insert(1, path_constants)
    from hBn_PP import hBn_lambda_p
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


b = -0.01

d_nano = 1
d_micro = d_nano*1e-3


def lambda_p(energy0):
    
    E = energy0
    

    
#    d_micros = d_nano*1e-3
#    Rp = Silver_Rp(E,epsi1,epsi3)
    lambda_p_v = hBn_lambda_p(E,epsi1,epsi3)*d_nano


    return lambda_p_v


dir_dip_moment = 'all'
#energy0_pol = 43
#omega0 = energy0_pol*1e-3/hb 
##R = 10 # 10 nm en unidades de micrometros
#kappa_factor_omega0 = 0.1
#kappa_r_factor= 0.5
 

title4 =   r'b = %i nm, d = %i nm' %(b*1e3,d_nano) 

N = 30


#%%

def function_imag_ana(energy0,int_v,zp_nano):
    omegac0 = energy0/aux 
    zp = zp_nano*1e-3

    rta = EELS_film_ana_f(omegac0,epsi1,epsi3,d_nano,int_v,b,zp,dir_dip_moment)
    
    return rta


def function_imag_num(energy0,int_v,zp_nano):
    omegac0 = energy0/aux 
    zp = zp_nano*1e-3

    rta = EELS_film_num_f(omegac0,epsi1,epsi3,d_nano,int_v,b,zp,dir_dip_moment)
    
    return rta


def function_imag_pole_aprox(energy0,int_v,zp_nano):
    omegac0 = energy0/aux 
    zp = zp_nano*1e-3

    rta = EELS_film_pole_aprox_f(omegac0,epsi1,epsi3,d_nano,int_v,b,zp,dir_dip_moment)
    
    return rta



#%%
    

if plot_vs_c == 1 :
    E0 = 44 # meV
    # z0 = 0.06*1e3
    zp0 = 0.05*1e3 
    
    labelx = r'v/c'   
    title4 = title4 + ', ' + r'$\hbar\omega$ = %.2f eV, $z_p$ = %i nm' %(E0,zp0)
    label1 = 'vs_v' 
#    listx = np.linspace(0.0001,2,N)
    listx = np.linspace(0.005,0.12,N)

if plot_vs_E ==1 :
    # z0 = 0.06*1e3
    int_v0 = 10
    zp0 = 0.05*1e3 
    
    labelx = r'$\hbar\omega$ [eV]'   
    title4 = title4 + ', ' + r'v = c/%i, $z_p$ = %i nm' %(int_v0,zp0)
    label1 = 'vs_E' 
#    listx = np.linspace(0.0001,2,N)
    listx = np.linspace(15,65,N)


x1 = 0.09260651629072682 
x2 = 0.10112781954887218
x3 = 0.17030075187969923
x4 = 0.19937343358395992




if plot_vs_zp == 1 : 
    E0 = 0.171 # eV
#    E0 = 0.171
#    E0 = 0.1875
    
    int_v0 = 10
    lambbda_p = np.real(lambda_p(E0))

    labelx = r'$z_p$ [nm]'   
    title4 = title4 + ', ' + r'v = c/%i, $\hbar\omega$ = %.3f eV, $\lambda_p$ = %.1f nm' %(int_v0,E0,lambbda_p)
    label1 = 'vs_zp_E%imeV' %(E0*1e3) 
#    listx = np.linspace(0.0001,2,N)
    if d_nano == 1:
        if E0 <= 0.187:
            listx = np.linspace(250,750,N)
        else:
            listx = np.linspace(1,100,N)

    else:
        listx = np.linspace(250,800,N)
        
    listx = np.linspace(1,600,N)
    
#    listx = np.linspace(250,750,N)
#    
    
    
title =  title4 

#%%
    
tamfig = (4.5,3.5)
tamlegend = 12
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


listy_im_ana = []
listy_im_num = []
listy_im_pole_aprox = []

if plot_vs_E == 1: 

    for value in listx: 

        y_im_ana = function_imag_ana(value,int_v0,zp0)        
        listy_im_ana.append(y_im_ana)
        
        y_im_num = function_imag_num(value,int_v0,zp0)        
        listy_im_num.append(y_im_num)
        
        y_im_pole_aprox = function_imag_pole_aprox(value,int_v0,zp0)        
        listy_im_pole_aprox.append(y_im_pole_aprox)
        
               
elif plot_vs_c == 1:       

    
    for value in listx: 
        
        value2 = 1/value
#        print(value2)

        y_im_ana = function_imag_ana(E0,value2,zp0)        
        listy_im_ana.append(y_im_ana)

        y_im_num = function_imag_num(E0,value2,zp0)        
        listy_im_num.append(y_im_num)
        
        y_im_pole_aprox = function_imag_pole_aprox(E0,value2,zp0)        
        listy_im_pole_aprox.append(y_im_pole_aprox)

        
elif plot_vs_zp == 1:

    
    for value in listx: 

        y_im_ana = function_imag_ana(E0,int_v0,value)        
        listy_im_ana.append(y_im_ana)
        
        y_im_num = function_imag_num(E0,int_v0,value)        
        listy_im_num.append(y_im_num)
        
        y_im_pole_aprox = function_imag_pole_aprox(E0,int_v0,value)        
        listy_im_pole_aprox.append(y_im_pole_aprox)


#%%
        
peaks, _ = find_peaks(listy_im_ana)
maxi = listx[peaks]
listy_aux  = np.linspace(np.min(listy_im_ana), np.max(listy_im_ana), 10)
print(maxi)
if len(maxi ) > 1 :
    listy_im_ana = np.array(listy_im_ana)
    list_maxis_y = listy_im_ana[peaks]
    
    maxi_ind = np.argmax(list_maxis_y)
    maxi =  listx[peaks[maxi_ind]]
    print(E0,maxi)

graph(title,labelx,'$\Gamma_{film}$ [$\mu$s]',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(listx,listy_im_ana,'.-',ms = ms,color = 'purple', label = 'PP analytical')
plt.plot(listx,listy_im_num,'.-',ms = ms,color = 'lightseagreen',label = 'full numerical')
plt.plot(listx,listy_im_pole_aprox,'.-',ms = ms,color = 'darkred',label = 'PP numerical')
#plt.plot(np.ones(10)*maxi, listy_aux,'-k')
plt.legend(loc = 'best',markerscale=1.5,fontsize=tamlegend,frameon=0.1,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
#if plot_vs_c == 1:
#    plt.yscale('log')
os.chdir(path_save)
plt.savefig( 'EELS' + label1 + '.png', format='png')   



#%%


