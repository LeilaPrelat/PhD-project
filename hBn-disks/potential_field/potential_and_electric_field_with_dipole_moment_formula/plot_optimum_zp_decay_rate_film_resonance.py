
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
load_data = 1
save_data = 0


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
    from decay_rate_film_resonance import EELS_film_ana_f,EELS_film_num_f
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

d_nano = 1
d_micro = d_nano*1e-3


int_v = 10

dir_dip_moment = 'all'


title4 =  r'b = %i nm, d = %i nm, v = c/%i, hBN disks, D = 120 nm' %(b*1e3,d_nano,int_v)
labelp = r'_d%inm_p%s' %(d_nano,dir_dip_moment)

x1 = 0.09260651629072682 
x2 = 0.10112781954887218

x3 = 0.17030075187969923
x4 = 0.19937343358395992

N = 60
listx = np.linspace(0.09,0.195,N)  
if primer_intervalo == 1:
    listx = np.linspace(x1*1.01,x2*0.99,N)
else:
    listx = np.linspace(0.171,0.195,N)

#%%
    
def function_imag_ana(energy0): ## devuelve el zp optimo en nanometros
    
    omegac0 = energy0/aux
    listx = np.linspace(0.1,800,500)
    
    if d_nano == 10:
        listx = np.linspace(10,800,500)
        
    listy = []
    for zp_nano in listx:
        zp = zp_nano*1e-3
        rta = EELS_film_ana_f(omegac0,epsi1,epsi3,d_nano,int_v,b,zp,dir_dip_moment)
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

labelx = r'Plasmon energy $\hbar\omega$ (eV)'
labely = r'optimal $z_p$ [nm]'
    
    
#%%

    
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
    if save_data == 1:
        plt.title(title,fontsize=int(tamtitle*0.9))

    return   

#%%

list_lambda_p = []
for x in listx:
    list_lambda_p.append(np.real(lambda_p(x)))
    
if save_data == 1:

    listy = []
    for x in listx:
        listy.append(function_imag_ana(x))
    
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
    
    info =  title4
    
    tabla = np.array([listx,listy])
    tabla = np.transpose(tabla)
    info = 'zp_optimum_for_decay_rate_resonance' + labelp
    header1 = 'E [eV]     zp [nm]' + ', ' + info + ', ' + name_this_py
    np.savetxt('zp_optimum_for_decay_rate' + labelp + '.txt', tabla, fmt='%1.11e', delimiter='\t', header = header1)


#%%

if load_data == 1:
    os.chdir(path_save)
    tabla = np.loadtxt('zp_optimum_for_decay_rate' + labelp + '.txt', delimiter='\t', skiprows=1)
    tabla = np.transpose(tabla)
    [listx,listy] = tabla
    
    labely = 'Surface-dipole distance \n for power scattered \n improvement ($\mu$m)'
    labely = 'Optimal dipole-surface \n separation ($\mu$m)'
        
    graph(title4,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    plt.plot(listx,np.array(listy)*1e-3,'-',ms = ms,color = 'purple')
    plt.plot(listx,np.array(list_lambda_p)*1e-3,'--',ms = ms,color = 'lightseagreen')
#    plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=0.1,handletextpad=0.2, handlelength=length_marker)
    plt.tight_layout()

    plt.savefig( 'zp_optimum_for_decay_rate' + labelp + '.png', format='png', dpi=dpi)  
    
    
    
    
    


