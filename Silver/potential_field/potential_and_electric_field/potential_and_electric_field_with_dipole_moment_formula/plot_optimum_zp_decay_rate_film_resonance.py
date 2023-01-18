
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
#import seaborn as sns
#sns.set()

#%%
paper = 1
create_data = 0
load_data = 1

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/potential_field/potential_and_electric_field/potential_and_electric_field_with_dipole_moment_formula','')
path_save = path_basic + '/' + 'optimum_zp'
if not os.path.exists(path_save):
    print('Creating folder to save graphs')
    os.mkdir(path_save)

err = 'decay_rate_film.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from decay_rate_film_resonance import  EELS_film_ana_f, EELS_film_num_f, EELS_film_pole_aprox_f
except ModuleNotFoundError:
    print(err)
    
try:
    sys.path.insert(1, path_constants)
    from Silver_PP import Silver_lambda_p
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

hbgamma_DL = 21*1e-3
hbomega_DL = 9.17 
epsilon_b = 4

epsi1,epsi3 = 1,1

b = -0.01

d_nano = 1
int_v = 10
 
#title1 = r'$\kappa$ = %.2f$\omega_0$, $\kappa_r$ = %.2f$\kappa$, $E_0$=%i meV' %(kappa_factor_omega0, kappa_r_factor, energy0_pol)     
#title2 = r'$\hbar\mu$ = %.2feV, $\hbar\gamma$ = %.4feV' %(hbmu,hbgama) 
#title3 = r'$z_p$=%inm, px=%i, py=%i, pz=%i' %(zp*1e3,px,py,pz)
title1 = r'b = %i nm, v = c/%i, d = %i nm, $\hbar\omega$ = %.2f eV' %(b*1e3,int_v, d_nano,hbomega_DL)
title2 = r'$\hbar\gamma_{in}$ = %i meV, $\epsilon_b$ = %i' %(hbgamma_DL*1e3,epsilon_b)
labelp = r'_d%inm_v%i' %(d_nano,int_v)
title = title1 + '\n' + title2


N = 100
list_freq = np.linspace(1,3.5,100)
list_zp = np.linspace(0.01,250,500)

#%%


def function_imag_ana(energy0,list_zp_nano): ## devuelve el zp optimo en nanometros
    omegac0 = energy0/aux
    
#    if energy0 < 0.4:
#        list_zp = np.linspace(0.1,14,200)
#    else:
#        list_zp = np.linspace(0.1,8,200)
        
    
    listy = []
    for zp_nano in list_zp_nano:
        zp = zp_nano*1e-3
        rta = EELS_film_ana_f(omegac0,epsi1,epsi3,d_nano,int_v,b,zp)
        listy.append(rta)
#    print(energy0,v_sobre_c)
    
    peaks, _ = find_peaks(listy, height=0)
    maxi = list_zp_nano[peaks]
    print(energy0, maxi)
    if len(maxi ) > 1 :
        listy = np.array(listy)
        list_maxis_y = listy[peaks]
        
        maxi_ind = np.argmax(list_maxis_y)
        maxi =  list_zp_nano[peaks[maxi_ind]]
        print(maxi)
    print('')
    
        
        
    return float(maxi)

    
#%%    
    
def lambda_p(energy0):
    
    E = energy0

#    d_micros = d_nano*1e-3
#    Rp = Silver_Rp(E,epsi1,epsi3)
    lambda_p_v = Silver_lambda_p(E,epsi1,epsi3)*d_nano


    return lambda_p_v

    
#%%
    

labelx = r'Plasmon energy $\hbar\omega$ (eV)'
labely = r'optimal $z_p$ [nm]'
    
tamfig = [2.5, 2]
tamletra = 7
tamtitle  = 8
tamnum = 6
tamlegend = 6
labelpady = 2
labelpadx = 3
pad = 2.5
mk = 1
ms = 2
hp = 0.3
length_marker = 1.5
dpi = 500

def graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad):
    plt.figure(figsize=tamfig)
    plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
    plt.ylabel(labely,fontsize=tamletra,labelpad =labelpadx)
    plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in", pad = pad)
#    plt.title(title,fontsize=int(tamtitle*0.9))

    return   

#%%

if create_data == 1:

    listy = []
    for x in list_freq:
        listy.append(function_imag_ana(x,list_zp))


    list_lambda_p = []
    for x in list_freq:
        list_lambda_p.append(np.real(lambda_p(x)))
    
    graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    plt.plot(list_freq,listy,'.-',ms = ms,color = 'purple',label = r'opt $z_p$')
    plt.plot(list_freq,list_lambda_p,'.-',ms = ms,color = 'lightseagreen',label = r'$\lambda_p$')
    plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=0.1,handletextpad=0.2, handlelength=length_marker)
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'zp_optimum_for_decay_rate_Silver_resonance' + labelp + '.png', format='png')  



    if d_nano!= 1:

        hspace = 0.12
        wspace = 0.1
        loc2 = [0.3,0.88]    # R grande modos 1 
        
        fig,axs = plt.subplots(2,1, sharex=True, facecolor='w', figsize = (5,3.5))
        plt.subplots_adjust(hspace =hspace,wspace = wspace)

        axs[0].plot(list_freq,list_lambda_p,'.-',ms = ms,color = 'lightseagreen',label =  r'$\lambda_p$')
        axs[1].plot(list_freq,listy,'.-',ms = ms,color = 'purple',label = r'opt $z_p$') 
        
        for i in [0,1]:
            axs[i].minorticks_on()
            axs[i].tick_params(labelsize = tamnum,pad = pad)
            axs[i].set_ylabel(labely,fontsize = tamletra,labelpad=labelpady) 
        
    
        axs[1].set_xlabel(labelx,fontsize = int(tamletra*1.1),labelpad=labelpadx)
        
        fig.legend(loc = loc2, ncol = 4,markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.5)
#        fig.tight_layout()
        plt.savefig( 'zp_optimum_for_decay_rate_resonance_sep' + labelp + '.png', format='png')  

    tabla = np.array([list_freq,listy,list_lambda_p])
    tabla = np.transpose(tabla)
    header1 = 'E [meV]     zp [nm]     Re(lambda_p) [nm]' + ', ' + title + ', ' + name_this_py
    np.savetxt('zp_optimum_for_decay_rate_resonance' + labelp + '.txt', tabla, fmt='%1.11e', delimiter='\t', header = header1)

#%%

if load_data == 1:
    from scipy.signal import savgol_filter
    
    os.chdir(path_save)
    tabla = np.loadtxt('zp_optimum_for_decay_rate_resonance' + labelp + '.txt' , delimiter='\t', skiprows=1)
    tabla = np.transpose(tabla)
    [listx,listy,list_lambda_p] = tabla

    labely = 'Surface-dipole distance \n for power scattered \n improvement ($\mu$m)'
    labely = 'Optimal dipole-surface \n separation ($\mu$m)'


    list_zp_div_lambda_p = np.array(listy)/np.array(list_lambda_p)
    
    labely = 'Surface-dipole' + '\n' +  r'distance  $z_{\rm 0}$/$\lambda_{\rm p}$'
    labelx = r'Frequency $\omega$/$\omega_{\rm D}$'
    
    omega_D = 9.17  ## eV sin el hbar porque se cancela con el listx 
    
    listx_2 = np.array(listx)/omega_D ## sin el hbar porque se cancela con el listx 
    
    list_zp_div_lambda_p = savgol_filter(list_zp_div_lambda_p, 71, 3)

    graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    plt.plot(listx,np.array(list_zp_div_lambda_p),'-',ms = ms,color = 'purple')
#    plt.plot(listx,np.array(list_lambda_p)*1e-3,'--',ms = ms,color = 'lightseagreen')
#    plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=0.1,handletextpad=0.2, handlelength=length_marker)
    plt.tight_layout()
    os.chdir(path_save)
    
    plt.savefig( 'zp_optimum_for_decay_rate_Silver_resonance' + labelp + '.png',bbox_inches='tight',pad_inches = 0.01, format='png', dpi=dpi)  
#plt.yscale('log')



#


#

#



