
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
load_data = 0
save_data = 1


primer_intervalo = 0

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/potential_field/potential_and_electric_field_with_dipole_moment_formula','')
path_save = path_basic + '/' + 'optimum_zp'
if not os.path.exists(path_save):
    print('Creating folder to save graphs')
    os.mkdir(path_save)

err = 'decay_rate_film3.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from decay_rate_film import EELS_film_ana_f,EELS_film_num_f
except ModuleNotFoundError:
    print(err)
    
try:
    sys.path.insert(1, path_constants)
    from hBn_PP import hBn_lambda_p
except ModuleNotFoundError:
    print('hBn_PP.py no se encuentra en ' + path_basic)

    
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

d_nano = 10
d_micro = d_nano*1e-3

energy0_pol = 0.2 ## eV
omega0 = energy0_pol/hb 
#R = 10 # 10 nm en unidades de micrometros
kappa_factor_omega0 = 0.1
kappa_r_factor= 0.5

int_v = 10

dir_dip_moment = 'all'


#title1 = r'$\kappa$ = %.2f$\omega_0$, $\kappa_r$ = %.2f$\kappa$, $E_0$=%i meV' %(kappa_factor_omega0, kappa_r_factor, energy0_pol)     
#title2 = r'$\hbar\mu$ = %.2feV, $\hbar\gamma$ = %.4feV' %(hbmu,hbgama) 
#title3 = r'$z_p$=%inm, px=%i, py=%i, pz=%i' %(zp*1e3,px,py,pz)
title1 = r'$\kappa$ = %.2f$\omega_0$, $\kappa_r$ = %.2f$\kappa$, $\hbar\omega_o$=%.2f eV' %(kappa_factor_omega0, kappa_r_factor, energy0_pol) 
title4 = title1 + '\n' + r'b = %i nm, d = %i nm, v = c/%i, hBN disks, D = 120 nm' %(b*1e3,d_nano,int_v)
labelp = r'_d%inm_p%s' %(d_nano,dir_dip_moment)

x1 = 0.09260651629072682 
x2 = 0.10112781954887218

x3 = 0.17030075187969923
x4 = 0.19937343358395992

N = 25
listx = np.linspace(0.09,0.195,N)  
if primer_intervalo == 1:
    listx = np.linspace(x1*1.01,x2*0.99,N)
else:
    listx = np.linspace(0.171,0.195,N)

#%%
    
def function_imag_ana(energy0): ## devuelve el zp optimo en nanometros
    
    omegac0 = energy0/aux
    listx = np.linspace(0.01,80,800)
    
    if d_nano == 10:
        listx = np.linspace(1,40,200)
        
    listy = []
    for zp_nano in listx:
        zp = zp_nano*1e-3
        rta = EELS_film_ana_f(omegac0,epsi1,epsi3,d_nano,int_v,b,zp,omega0,kappa_factor_omega0,kappa_r_factor,dir_dip_moment)
        listy.append(rta)
#    print(energy0,v_sobre_c)
    
    peaks, _ = find_peaks(listy)
    maxi = listx[peaks]
    print(energy0, maxi)
    if len(maxi ) > 1 :
        listy = np.array(listy)
        list_maxis_y = listy[peaks]
        
        maxi_ind = np.argmax(list_maxis_y)
        maxi =  listx[peaks[maxi_ind]]
        print(maxi)
#        if listx[]
        
        
    
    return float(maxi)


def lambda_p(energy0):
    
    E = energy0
    

    
#    d_micros = d_nano*1e-3
#    Rp = Silver_Rp(E,epsi1,epsi3)
    lambda_p_v = hBn_lambda_p(E,epsi1,epsi3)*d_nano


    return lambda_p_v

labelx = r'$\hbar\omega$ [eV]'
labely = r'optimal $z_p$ [nm]'
    
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

if save_data == 1:

    listy = []
    list_lambda_p = []
    for x in listx:
        listy.append(function_imag_ana(x))
        list_lambda_p.append(np.real(lambda_p(x)))
    
    graph(title4,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    plt.plot(listx,listy,'.-',ms = ms,color = 'purple',label = r'opt $z_p$')
    plt.plot(listx,list_lambda_p,'.-',ms = ms,color = 'lightseagreen',label = r'$\lambda_p$')
    plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=0.1,handletextpad=0.2, handlelength=length_marker)
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'zp_optimum_for_decay_rate' + labelp + '.png', format='png')  


    graph(title4,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    plt.plot(listx,listy,'.-',ms = ms,color = 'purple',label = r'opt $z_p$')
#    plt.plot(listx,list_lambda_p,'.-',ms = ms,color = 'lightseagreen',label = r'$\lambda_p$')
    plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=0.1,handletextpad=0.2, handlelength=length_marker)
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'zp_optimum_for_decay_rate_v2' + labelp + '.png', format='png')    
    
    ### change some numbers 
    
    info = title1 + ', ' + title4
    
    tabla = np.array([listx,listy])
    tabla = np.transpose(tabla)
    info = 'zp_optimum_for_decay_rate_resonance' + labelp
    header1 = 'E [eV]     zp [nm]' + ', ' + info + ', ' + name_this_py
    np.savetxt('zp_optimum_for_decay_rate' + labelp + '.txt', tabla, fmt='%1.11e', delimiter='\t', header = header1)


#%%

if load_data == 1:

    tabla = np.loadtxt('zp_optimum_for_decay_rate_hBn_d%inm' %(d_nano) + '.txt', delimiter='\t', skiprows=1)
    tabla = np.transpose(tabla)
    [listx,listy] = tabla
    
    
    graph(title4,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    plt.plot(listx,listy,'.-',ms = ms,color = 'purple',label = r'opt $z_p$')
    plt.plot(listx,list_lambda_p,'.-',ms = ms,color = 'lightseagreen',label = r'$\lambda_p$')
    plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=0.1,handletextpad=0.2, handlelength=length_marker)
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'zp_optimum_for_decay_rate' + labelp + '.png', format='png')  
    
    
    
    
    


