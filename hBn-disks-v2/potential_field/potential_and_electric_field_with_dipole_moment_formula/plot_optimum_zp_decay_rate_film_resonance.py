
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


create_data = 0
load_data = 1

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
    from decay_rate_film_resonance import EELS_film_ana_f, EELS_film_pole_aprox_f
except ModuleNotFoundError:
    print(err)



try:
    sys.path.insert(1, path_constants)
    from hBn_PP import epsilon_x,hBn_lambda_p
except ModuleNotFoundError:
    print('epsilon_z.py no se encuentra en ' + path_constants)


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


b = -0.01

d_nano_film = 1

D_disk_nano = 100
d_thickness_disk_nano = 1
#title1 = r'$\kappa$ = %.2f$\om

int_v = 10
 
#title1 = r'$\kappa$ = %.2f$\omega_0$, $\kappa_r$ = %.2f$\kappa$, $E_0$=%i meV' %(kappa_factor_omega0, kappa_r_factor, energy0_pol)     
#title2 = r'$\hbar\mu$ = %.2feV, $\hbar\gamma$ = %.4feV' %(hbmu,hbgama) 
#title3 = r'$z_p$=%inm, px=%i, py=%i, pz=%i' %(zp*1e3,px,py,pz)
title1 = r'b = %i nm, v = c/%i, d = %i nm, hBN' %(b*1e3,int_v, d_nano_film)
#title2 = r'$\hbar\gamma_{in}$ = %i meV, $\epsilon_b$ = %i' %(hbgamma_DL*1e3,epsilon_b)
labelp = r'_dfilm%.1fnm_ddisk%.1fnm_v%i' %(d_nano_film,d_thickness_disk_nano,int_v)
title = title1 

x1 = 0.09260651629072682 
x2 = 0.10112781954887218

x3 = 0.17030075187969923
x4 = 0.19937343358395992


#### la parte real de lambda_p para d = 0.4 nm es positiva si #####

# 0.08684904004367106 intervalo 1   # 0.093 - 0.1 esta dentro de x1-x2 pero no funciona porque Im{Gself} es negativo ---> el decay rate es negativo 
# 0.13150932790273412

# 0.15239881738519911 intervalo 2   # 0.153 - 0.169 esta dentro de x3-x4 . funciona desde 0.141 eV hasta 0.162 eV (entre 0.163 y 0.169 Im{Gself} es negativo ---> el decay rate es negativo )
# 0.16968667074999771

# 0.20714368637372804 intervalo 3 # es mayor a x4 
# 0.8064559363534132


#########################

N = 50
listx = np.linspace(0.09,0.195,N)  
#if primer_intervalo == 1:
#    listx = np.linspace(x1*1.01,x2*0.99,N)
#else:
#    listx = np.linspace(0.171,0.195,N)
#
#listx = np.linspace(0.095 , 0.162, N )
   
listx = np.linspace(0.1525 , 0.162, N ) ## para d = 0.4 nm y primer intervalo 

listx = np.linspace(0.093 , 0.1625, N ) ## para d = 0.1 nm y primer intervalo 

#listx = np.linspace(0.093 , 0.1, N ) ## para d = 0.1 nm y primer intervalo 

listx = np.linspace(0.13 , 0.1625, N ) ## para d = 1 nm y primer intervalo (largo )
listx = np.linspace(0.09 , 0.14, 25) ## para d = 1 nm y primer intervalo (una parte, un zoom )

# con el cambio de signo en Gself y dip moment, ahora se puede usar el segundo intervalo (entre x3 y x4)
listx  = np.linspace(0.171 , 0.195, 100)

#%%

def lambda_p(energy0):
    
    E = energy0
    

    
#    d_micros = d_nano*1e-3
#    Rp = Silver_Rp(E,epsi1,epsi3)
    lambda_p_v = hBn_lambda_p(E,epsilon_Silica(energy0),epsilon_Silica(energy0))*d_nano_film


    return lambda_p_v


def lambda_p_v2(energy0): ## creo que esta no tiene mucho que ver 
    
    epsi_x = epsilon_x(energy0)
    epsi_HBN_par = epsi_x
    epsi_silica = epsilon_Silica(energy0)        
    
    omegac = energy0/aux
#    d_micros = d_nano*1e-3
    d_micro = d_nano_film*1e-3
    alfa_p = epsi_silica*2/(omegac*d_micro*(epsi_HBN_par-1))
    kp = alfa_p*omegac
      

    return (2*np.pi/kp)*1e3 ## en nano

def function_imag_ana(energy0): ## devuelve el zp optimo en nanometros
    omegac0 = energy0/aux
#    if primer_intervalo == 1:
#        listx = np.linspace(50,450,100)
#    else:
#        listx = np.linspace(200,600,100)
    N = 400
    N = 350
#    if d_nano_film == 1:    
##        if energy0 <= 0.187:
##            listx = np.linspace(300,700,N) # segundo intervalo
##        else:
##            listx = np.linspace(100,500,N) # segundo intervalo
#
#        if energy0 <= 0.187:
#            list_zp_nano = np.linspace(50,750,N)
#        else:
#            list_zp_nano = np.linspace(100,400,N)    
#    elif d_nano_film == 1:
#        list_zp_nano = np.linspace(90,290,N) ## para todo el intervalo 
#        list_zp_nano = np.linspace(150,250,N) ## para una parte del intervalo (un zoom)
#    else:
#        list_zp_nano = np.linspace(10,400,N)
#        list_zp_nano = np.linspace(10,200,N)
        
    if D_disk_nano == 50:
        
        list_zp_nano = np.linspace(0.001,3.5,N)
    else:
        list_zp_nano = np.linspace(1,10,N)#         
        list_zp_nano = np.linspace(2,11,N) # mas preciso          
        
    listy = []
    for zp_nano in list_zp_nano:
        zp = zp_nano*1e-3
        rta = EELS_film_pole_aprox_f(omegac0,epsilon_Silica, d_nano_film, d_thickness_disk_nano, D_disk_nano,int_v,b,zp)
        listy.append(rta)
#    print(energy0,v_sobre_c)
        
    lambda_p_value = lambda_p(energy0)    
    peaks, _ = find_peaks(listy, height=0)
    maxi = list_zp_nano[peaks]
    
    print('energy: ', energy0, 'eV', ' zp: ' , maxi, 'nm', r'$\lambda_p$: ', lambda_p_value, 'nm')
    
    
    if len(maxi ) > 1 :
        listy = np.array(listy)
        list_maxis_y = listy[peaks]
        
        maxi_ind = np.argmax(list_maxis_y)
        maxi =  list_zp_nano[peaks[maxi_ind]]
        print(maxi,)
    print('')


   
    return float(maxi)


    
#%%    
    



#%%
    
labelx = r'Plasmon energy $\hbar\omega$ (eV)'
labely = r'optimal $z_p$ [nm]'
    
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
    plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
    plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in", pad = pad)
#    plt.title(title,fontsize=int(tamtitle*0.9))

    return   

#%%
if create_data == 1:
    
    listy = []
    for x in listx:
        listy.append(function_imag_ana(x))   

    list_lambda_p = []
    for x in listx:
        list_lambda_p.append(np.real(lambda_p(x)))
        

    
    graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    plt.plot(listx[0:len(listy)],listy,'.-',ms = ms,color = 'purple',label = r'opt $z_p$')
#    plt.plot(listx,list_lambda_p,'.-',ms = ms,color = 'lightseagreen',label = r'$\lambda_p$')
    plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=0.1,handletextpad=0.2, handlelength=length_marker)
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'zp_optimum_for_decay_rate_hBN_disks_resonance' + labelp + '.png', format='png', dpi=dpi)  
    
    
    tabla = np.array([listx,listy,list_lambda_p])
    tabla = np.transpose(tabla)
    header1 = 'E [eV]     zp [nm]     Re(lambda_p) [nm]' + ', ' + title + ', ' + name_this_py
    np.savetxt('zp_optimum_for_decay_rate_hBN_disks_resonance' + labelp + '.txt', tabla, fmt='%1.11e', delimiter='\t', header = header1)

#%%
#
if load_data == 1:

    from scipy.signal import savgol_filter
    os.chdir(path_save)
    tabla = np.loadtxt('zp_optimum_for_decay_rate_hBN_disks_resonance' + labelp + '.txt', delimiter='\t', skiprows=1)
    tabla = np.transpose(tabla)
    [listx,listy,list_lambda_p] = tabla
    
    labely = 'Surface-dipole distance \n for power scattered \n improvement ($\mu$m)'
    labely = 'Optimal dipole-surface \n separation ($\mu$m)'

    list_lambda_p = []
    for x in listx:
        list_lambda_p.append(np.real(lambda_p(x)))


    
    list_zp_div_lambda_p = []
    for ind in range(len(listx)):
        rta = listy[ind]/list_lambda_p[ind]
        list_zp_div_lambda_p.append(rta)
    
    labely = 'Surface-dipole' + '\n' +  r'distance  $z_{\rm 0}$/$\lambda_{\rm p}$'
    labelx = r'Frequency $\omega/\omega_\parallel$'
    
    omega_D = 170.1*1e-3  ## eV sin el hbar porque se cancela con el listx 
    
    listx_2 = np.array(listx)/omega_D ## sin el hbar porque se cancela con el listx 
    if d_nano_film == 0.1:
        list_zp_div_lambda_p = savgol_filter(list_zp_div_lambda_p, 27, 3)
        lim = -1
        
    elif d_nano_film == 1:
        list_zp_div_lambda_p = savgol_filter(list_zp_div_lambda_p, 15, 3)
        lim = -1

    graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    plt.plot(listx_2[0:lim],np.array(list_zp_div_lambda_p[0:lim]),'-',ms = ms,color = 'purple')
#    plt.plot(listx,np.array(list_lambda_p)*1e-3,'--',ms = ms,color = 'lightseagreen')
#    plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=0.1,handletextpad=0.2, handlelength=length_marker)
    plt.tight_layout()
    os.chdir(path_save)
    
    plt.savefig( 'zp_optimum_for_decay_rate_hBN_disks_resonance' + labelp + '.png',bbox_inches='tight',pad_inches = 0.01, format='png', dpi=dpi)  

#
#    graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#    plt.plot(listx_2,np.array(listy),'-',ms = ms,color = 'purple')
#    plt.plot(listx_2,np.array(list_lambda_p),'--',ms = ms,color = 'lightseagreen')
##    plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=0.1,handletextpad=0.2, handlelength=length_marker)
#    plt.tight_layout()

#



