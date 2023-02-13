
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
path_constants =  path_basic.replace('/potential_field/potential_and_electric_field/potential_and_electric_field_with_dipole_moment_formula','')
path_save = path_basic + '/' + 'decay_rate_film'
if not os.path.exists(path_save):
    print('Creating folder to save graphs')
    os.mkdir(path_save)

err = 'decay_rate_film3.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from decay_rate_film3_resonance import EELS_film_ana_f_div_gamma0,EELS_film_ana_f_div_gamma0_v3
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

epsi1,epsi2 = 1,1
hbmu,hbgama = 0.3,0.0001
zp = 0.05
b = -0.01


 
#title1 = r'$\kappa$ = %.2f$\omega_0$, $\kappa_r$ = %.2f$\kappa$, $E_0$=%i meV' %(kappa_factor_omega0, kappa_r_factor, energy0_pol)     
#title2 = r'$\hbar\mu$ = %.2feV, $\hbar\gamma$ = %.4feV' %(hbmu,hbgama) 
#title3 = r'$z_p$=%inm, px=%i, py=%i, pz=%i' %(zp*1e3,px,py,pz)
title4 = r'b = %i nm' %(b*1e3)
labelp = r'_res' 

N = 200



def function_imag_ana(energy0,int_v,zp_nano):
    omegac0 = energy0*1e-3/aux 
    zp = zp_nano*1e-3

    rta1 = EELS_film_ana_f_div_gamma0_v3(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,b,zp)
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


def lambda_p(energy0):
    
    E = energy0*1e-3
    
    omegac = E/aux
    
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    alfa_p = 1j*(epsi1 + epsi2)/(cond)
    kp = alfa_p*omegac

    return (2*np.pi/kp)*1e3 ## en nano


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
    E0 = 45 # meV
    int_v0 = 10

    labelx = r'Surface-dipole distance, $z_{\rm 0}$/$\lambda_{\rm p}$'   
    title4 = title4 + ', ' + r'v = c/%i, $\hbar\omega$ = %i meV' %(int_v0,E0)
    label1 = 'vs_zp' + labelp + '_E%imeV' %(E0)
#    listx = np.linspace(0.0001,2,N)
    listx = np.linspace(10,4000,N)
    
    print(minimum_function(E0,int_v0)*1e3)
    print(np.abs(minimum_function(E0,int_v0))*2*1e3)
    lambda_p_value = lambda_p(E0)


title =  title4 


#%%
    
tamfig = [2.5, 2]
tamletra = 9
tamtitle  = 8
tamnum = 7
tamlegend = 8
labelpady = 2
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

omega_omega_D = E0*1e-3/hbmu

#if E0 == 45 and int_v0 == 10:
#    zp_crit_lambda_p_value =  0.7361932072984527  ## for     E0 = 45 # meV  int_v0 = 10
#    zp_max_lambda_p_value = 0.63902942
#elif E0 == 43 and int_v0 == 10: 
#    zp_crit_lambda_p_value = 0.7268042897014039
#elif E0 == 35 and int_v0 == 10:     
#    zp_crit_lambda_p_value = 0.6759424385432607

#title1 = r'$z^{\rm crit}_{\rm 0}$/$\lambda_{\rm p}$ = %.2f'%(zp_crit_lambda_p_value)  ### from plot_decay_rate ... simpler 
#title2 = r'$z^{\rm opt}_{\rm 0}$/$\lambda_{\rm p}$ = %.2f' %(maxi2[0])
#
#title = title1 + ', ' + title2

graph(title,labelx,r'$\Gamma_{\rm SP}/\Gamma_{\rm EELS}$',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#plt.title(title,fontsize=int(tamtitle*0.9))
plt.plot(listx_2,np.array(listy_im_ana),'-',ms = ms,color = 'purple')
#plt.plot(listx,list_ana_parallel,'.-',ms = ms,color = 'darkred',label = r'$\Gamma_{\parallel}$')
#plt.plot(np.ones(10)*zp_crit_lambda_p_value, np.array(listy_aux)*1e-4,'--k')
#plt.plot([],[],'-w',label = r'$\omega/\omega_{\rm D}$=%.2f'%(omega_omega_D))
plt.text(0.8,0.041,r'$\omega/\omega_{\rm D}$ = %.2f'%(omega_omega_D),fontsize=tamlegend)
#plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=False,handletextpad=hp, handlelength=0)
plt.tight_layout()
#if plot_vs_c == 1:
#    plt.yscale('log')
os.chdir(path_save)
plt.savefig( 'EELS_film_' + label1 + 'omega_omega_D%.2f'%(omega_omega_D) + '.png', format='png',bbox_inches='tight',pad_inches = 0.01, dpi=dpi)   


tabla = np.array([np.real(np.array(listx_2)),np.real(np.array(listy_im_ana))])
tabla = np.transpose(tabla)
info_title = title4
info = 'decay rate normalized grafeno' 
header1 = 'z0/lambda_p     Gamma_SP/Gamma_EELS' + info_title + ', ' + ', '  +  info + ', ' + name_this_py
np.savetxt('decay_rate_normalized_grafeno_' + label1 + '.txt', tabla, fmt='%1.11e', delimiter='\t', header = header1)

