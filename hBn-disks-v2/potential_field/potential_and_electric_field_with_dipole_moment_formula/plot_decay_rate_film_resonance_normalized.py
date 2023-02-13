
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
from scipy.signal import find_peaks
#import seaborn as sns
#sns.set()

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
    from decay_rate_film_resonance import EELS_film_ana_f_div_gamma0,EELS_film_ana_f_div_gamma0_v3
except ModuleNotFoundError:
    print(err)
try:
    sys.path.insert(1, path_constants)
    from hBn_PP import hBn_lambda_p,epsilon_x
except ModuleNotFoundError:
    print(err)

try:
    sys.path.insert(1, path_basic)
    from Silica_epsilon import epsilon_Silica
except ModuleNotFoundError:
    print('Silica_epsilon.py no se encuentra en ' + path_basic)

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

epsi1 = 1
zp = 0.05
b = -0.01


d_nano_film = 1

D_disk_nano = 100
d_thickness_disk_nano = 1
#title1 = r'$\kappa$ = %.2f$\omega_0$, $\kappa_r$ = %.2f$\kappa$, $E_0$=%i meV' %(kappa_factor_omega0, kappa_r_factor, energy0_pol)     
#title2 = r'$\hbar\mu$ = %.2feV, $\hbar\gamma$ = %.4feV' %(hbmu,hbgama) 
#title3 = r'$z_p$=%inm, px=%i, py=%i, pz=%i' %(zp*1e3,px,py,pz)
title4 = r'b = %i nm' %(b*1e3)
labelp = r'_res_dfilm%.2fnm_ddisk%.2fnm_D%inm' %(d_nano_film,d_thickness_disk_nano,D_disk_nano) 

N = 350



def function_imag_ana(energy0,int_v,zp_nano):
    omegac0 = energy0/aux 
    zp = zp_nano*1e-3

    rta1 = EELS_film_ana_f_div_gamma0_v3(omegac0,epsilon_Silica, d_nano_film, d_thickness_disk_nano, D_disk_nano,int_v,b,zp)
#    rta2 = EELS_dir_ana_f(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,b,zp)
    
#    print(rta1)
    return rta1
#
#def function_imag_ana(energy0,int_v,zp_nano):
#    omegac0 = energy0*1e-3/aux 
#    zp = zp_nano*1e-3
#
#    rta = EELS_film_ana_f(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,b,zp)
#    
#    return rta


#def function_parallel_ana(energy0,int_v,zp_nano):
#    omegac0 = energy0*1e-3/aux 
#    zp = zp_nano*1e-3 
#    rta = EELS_parallel_ana_f(omegac0,epsi1,epsi2,hbmu,hbgama,L_nano,int_v,b,zp)
#    
#    return rta    
    

def lambda_p(energy0):
    
#    d_micros = d_nano*1e-3
    lambda_p_v = hBn_lambda_p(energy0,epsilon_Silica(energy0),epsilon_Silica(energy0))*d_nano_film

    return lambda_p_v ## en nano


def lambda_p_v2(energy0): ## creo que esta no tiene mucho que ver 
    
    epsi_x = epsilon_x(energy0)
    epsi_HBN_par = epsi_x
    epsi_silica = epsilon_Silica(energy0)        
    
    omegac = energy0/aux
#    d_micros = d_nano*1e-3
    d_micro = d_nano_film*1e-3
    alfa_p = epsi_silica*2/(omegac*d_micro*(epsi_HBN_par-1))
    kp = alfa_p*omegac
      

    return (2*np.pi/kp)*1e3 ## en micro

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
    
    labelx = r'$\hbar\omega$ [eV]'   
    title4 = title4 + ', ' + r'v = c/%i, $z_p$ = %i nm' %(int_v0,zp0)
    label1 = 'vs_E' + labelp
#    listx = np.linspace(0.0001,2,N)
    listx = np.linspace(15,65,N)

if plot_vs_zp == 1 : 
    int_v0 = 10 ## deberia ser 150 (disp relation) pero funciona con 10 <--- problema con la relacion de dispersion
    E0 = 0.173 # eV
 #   E0 = 0.095

    E0 = 0.173 
#    E0 = 0.175
#    E0 = 0.195
#    
    labelx = r'Surface-dipole distance, $z_{\rm 0}$/$\lambda_{\rm p}$'   
    title4 = title4 + ', ' + r'v = c/%i, $\hbar\omega$ = %i eV' %(int_v0,E0)
    label1 = 'vs_zp' + labelp + '_E%imeV' %(E0*1e3)
#    listx = np.linspace(0.0001,2,N)
    if d_nano_film == 1:
        if E0 <= 0.187:
            listx = np.linspace(10,650,N)
            listx = np.linspace(10,950,N)
        else:
            listx = np.linspace(25,400,N)

    else:
        listx = np.linspace(0.0005,3,N)
    
    if D_disk_nano == 50:
        
        listx = np.linspace(0.001,3.5,N)
    else:
        listx = np.linspace(1,10,N)#    print(minimum_function(E0,int_v0)*1e3)
#    print(np.abs(minimum_function(E0,int_v0))*2*1e3)
    if E0 == 0.171:
        listx = np.linspace(3,14,N)
        
    elif E0 == 0.175:
        listx = np.linspace(7.3,9.3,N)


    elif E0 == 0.195:
        listx = np.linspace(1,6,N)
    
    else:
        listx = np.linspace(2,9,N) ## mas general 
    
    lambda_p_value = lambda_p(E0)
    


title =  title4 


#%%
    
tamfig = [2.5, 2]
tamletra = 9
tamtitle  = 8
tamnum = 7
tamlegend = 8
labelpady = 3
labelpadx = 3
pad = 3
mk = 1
ms = 2
hp = 0.3
length_marker = 0
dpi = 500

def graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad):
    plt.figure(figsize=tamfig)
    plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
    plt.ylabel(labely,fontsize=tamletra,labelpad =labelpadx)
    plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in", pad = pad)
#    plt.title(title,fontsize=int(tamtitle*0.9))

    return   


#%%


listy_im_ana = []
list_ana_parallel = []


if plot_vs_E == 1: 

    for value in listx: 

        y_im_ana = function_imag_ana(value,int_v0,zp0)        
        listy_im_ana.append(y_im_ana)
        


elif plot_vs_c == 1:       

    for value in listx: 
        
        value2 = 1/value
#        print(value2)

        y_im_ana = function_imag_ana(E0,value2,zp0)        
        listy_im_ana.append(y_im_ana)



elif plot_vs_zp == 1:
    
    for value in listx: 

        y_im_ana = function_imag_ana(E0,int_v0,value)        
        listy_im_ana.append(y_im_ana)
        
     
    
#%%
        
peaks, _ = find_peaks(listy_im_ana, height=0)
maxi = listx[peaks]
listy_aux  = np.linspace(np.min(listy_im_ana), np.max(listy_im_ana), 10)
listx_2 = np.array(listx)/lambda_p_value
maxi2 = maxi/lambda_p_value


omegaD_silver = 170.1*1e-3
omega_omega_D = E0/(omegaD_silver)

graph(title,labelx,r'$\Gamma_{\rm SP}/\Gamma_{\rm EELS}$',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(listx_2,np.array(listy_im_ana),'-',ms = ms,color = 'purple')
#plt.plot(listx,list_ana_parallel,'.-',ms = ms,color = 'darkred',label = r'$\Gamma_{\parallel}$')
#plt.plot(np.ones(10)*maxi2, np.array(listy_aux)*1e-12,'-k',label = r'$z^{\rm opt}_{\rm o}$/$\lambda_{\rm p}$')
plt.plot([],[],'-w',label = r'$\omega/\omega_{\parallel}$=%.2f'%(omega_omega_D))
#plt.text(0.01,0.01,r'$\omega/\omega_{\rm D}$=%.2f'%(omega_omega_D),fontsize=tamlegend,transform=plt.transAxes)
plt.legend(loc = 'best',markerscale=0,fontsize=tamlegend,frameon=False,handletextpad=0, handlelength=0)
plt.tight_layout()
#if plot_vs_c == 1:
#    plt.yscale('log')
os.chdir(path_save)
plt.savefig( 'EELS_film_' + label1 + '.png', format='png',bbox_inches='tight',pad_inches = 0.01, dpi=dpi)   


