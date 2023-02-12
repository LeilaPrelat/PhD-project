
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

comparar_otro_hBN = 0

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_save = path_basic + '/' + 'Green_self_image'
if not os.path.exists(path_save):
    print('Creating folder to save graphs')
    os.mkdir(path_save)

err = 'fieldE_direct_numerical.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from green_self_image import green_self_ana_v1,green_self_ana_v2, green_self_pole_aprox_v1, green_self_pole_aprox_v2, green_self_num
except ModuleNotFoundError:
    print(err)

try:
    sys.path.insert(1, path_basic)
    from Silica_epsilon import epsilon_Silica
except ModuleNotFoundError:
    print(err)


try:
    sys.path.insert(1, path_basic)
    from constants import constantes
except ModuleNotFoundError:
    print('constants.py no se encuentra en ' + path_basic)

pi,hb,c,alfac,mu1,mu2 = constantes()
aux = c*hb
aux2 = 1e12/c

#%%

print('Definir parametros del problema')


zp_nano = 50


d_nano = 0.4

title = r'$z_{\rm p}$ = %i nm, d = %.2f nm' %(zp_nano,d_nano)

N = 100

    # z0 = 0.06*1e3
labelx = r'$\hbar\omega$ [eV]'   
labelp = 'vs_E_d%.2fnm' %(d_nano) 


x1 = 0.09260651629072682 
x2 = 0.10112781954887218


x3 = 0.17030075187969923
x4 = 0.19937343358395992

listx = np.linspace(0.087, 0.8, N)
listx = np.linspace(0.087, 0.2, N)
#listx = np.linspace(0.087, 0.195, N)
#listx = np.linspace(0.2, 0.8, N)


def function_num_xx_re(energy0):
    omegac0 = energy0/aux 
    zp = zp_nano*1e-3
    rtaself_x, rtaself_y, rtaself_z = green_self_num(omegac0,epsilon_Silica,d_nano,zp)
    # print(rta    
    return np.real(rtaself_x)

def function_num_xx_im(energy0):
    omegac0 = energy0/aux 
    zp = zp_nano*1e-3
    rtaself_x, rtaself_y, rtaself_z = green_self_num(omegac0,epsilon_Silica,d_nano,zp)
    # print(rta    
    return np.imag(rtaself_x)


def function_ana_xx_re_v1(energy0):
    omegac0 = energy0/aux 
    zp = zp_nano*1e-3
    rtaself_x, rtaself_y, rtaself_z  = green_self_ana_v1(omegac0,epsilon_Silica,d_nano,zp)
    
    return np.real(rtaself_x)


def function_ana_xx_re_v2(energy0):
    omegac0 = energy0/aux 
    zp = zp_nano*1e-3
    rtaself_x, rtaself_y, rtaself_z  = green_self_ana_v2(omegac0,epsilon_Silica,d_nano,zp)
    
    return np.real(rtaself_x) 



#def function_ana3_xx_re(energy0):
#    omegac0 = energy0/aux 
#    zp = zp_nano*1e-3
#    rtaself_x, rtaself_y, rtaself_z  = green_self_ana3(omegac0,epsilon_Silica,d_nano,zp)
#    
#    return np.real(rtaself_x) 

def function_ana_xx_im_v1(energy0):
    omegac0 = energy0/aux 
    zp = zp_nano*1e-3
    rtaself_x, rtaself_y, rtaself_z  = green_self_ana_v1(omegac0,epsilon_Silica,d_nano,zp)
    
    return np.imag(rtaself_x)



def function_ana_xx_im_v2(energy0):
    omegac0 = energy0/aux 
    zp = zp_nano*1e-3
    rtaself_x, rtaself_y, rtaself_z  = green_self_ana_v2(omegac0,epsilon_Silica,d_nano,zp)
    
    return np.imag(rtaself_x)



#def function_ana3_xx_im(energy0):
#    omegac0 = energy0/aux 
#    zp = zp_nano*1e-3
#    rtaself_x, rtaself_y, rtaself_z  = green_self_ana3(omegac0,epsilon_Silica,d_nano,zp)
#    
#    return np.imag(rtaself_x)

def function_pole_aprox_xx_re_v1(energy0):
    omegac0 = energy0/aux 
    zp = zp_nano*1e-3
    rtaself_x, rtaself_y, rtaself_z  = green_self_pole_aprox_v1(omegac0,epsilon_Silica,d_nano,zp)
    
    return np.real(rtaself_x)



def function_pole_aprox_xx_im_v1(energy0):
    omegac0 = energy0/aux 
    zp = zp_nano*1e-3

    rtaself_x, rtaself_y, rtaself_z  = green_self_pole_aprox_v1(omegac0,epsilon_Silica,d_nano,zp)
    
    return np.imag(rtaself_x)



def function_pole_aprox_xx_re_v2(energy0):
    omegac0 = energy0/aux 
    zp = zp_nano*1e-3
    rtaself_x, rtaself_y, rtaself_z  = green_self_pole_aprox_v2(omegac0,epsilon_Silica,d_nano,zp)
    
    return np.real(rtaself_x)



def function_pole_aprox_xx_im_v2(energy0):
    omegac0 = energy0/aux 
    zp = zp_nano*1e-3

    rtaself_x, rtaself_y, rtaself_z  = green_self_pole_aprox_v2(omegac0,epsilon_Silica,d_nano,zp)
    
    return np.imag(rtaself_x)



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

listy_re_ana_v1 = []
listy_re_ana_v2 = []
#listy_re_ana3 = []
listy_re_num = []
listy_re_pole_aprox_v1 = []
listy_re_pole_aprox_v2 = []

listy_im_ana_v1 = []
listy_im_ana_v2 = []
#listy_im_ana3 = []
listy_im_num = []
listy_im_pole_aprox_v1 = []
listy_im_pole_aprox_v2 = []

for value in listx: 

    y_re_ana_v1 = function_ana_xx_re_v1(value)       
    y_re_ana_v2 = function_ana_xx_re_v2(value)       
#    y_re_ana3 = function_ana3_xx_re(value)    
    y_re_num = function_num_xx_re(value)
    y_re_pole_aprox_v1 = function_pole_aprox_xx_re_v1(value)
    y_re_pole_aprox_v2 = function_pole_aprox_xx_re_v2(value)

    
    listy_re_ana_v1.append(y_re_ana_v1)
    listy_re_ana_v2.append(y_re_ana_v2)
#    listy_re_ana3.append(y_re_ana3)
    listy_re_num.append(y_re_num)
    listy_re_pole_aprox_v1.append(y_re_pole_aprox_v1)
    listy_re_pole_aprox_v2.append(y_re_pole_aprox_v2)


    y_im_ana_v1 = function_ana_xx_im_v1(value)       
    y_im_ana_v2 = function_ana_xx_im_v2(value)   
#    y_im_ana3 = function_ana3_xx_im(value)   
    y_im_num = function_num_xx_im(value)
    y_im_pole_aprox_v1 = function_pole_aprox_xx_im_v1(value)
    y_im_pole_aprox_v2 = function_pole_aprox_xx_im_v2(value)

    listy_im_ana_v1.append(y_im_ana_v1)
    listy_im_ana_v2.append(y_im_ana_v2)
#    listy_im_ana3.append(y_im_ana3)
    listy_im_num.append(y_im_num)
    listy_im_pole_aprox_v1.append(y_im_pole_aprox_v1)
    listy_im_pole_aprox_v2.append(y_im_pole_aprox_v2)    
    
    
#%%

if comparar_otro_hBN == 1:
    os.chdir(path_save)
    label1 = 'vs_E_d%.2fnm' %(d_nano) 
       
    tabla = np.loadtxt( 'Gself_hBN_num_' + label1 + '.txt', delimiter='\t', skiprows=1)
    tabla = np.transpose(tabla)
    [listx_num_S,listy_re_num_S,listy_im_num_S] = tabla 
    
    
    tabla = np.loadtxt( 'Gself_hBN_ana_' + label1 + '.txt', delimiter='\t', skiprows=1)
    tabla = np.transpose(tabla)
    [listx_ana_S,listy_re_ana2_S,listy_im_ana2_S] = tabla 
    
    
    tabla = np.loadtxt( 'Gself_hBN_pole_aprox_' + label1 + '.txt', delimiter='\t', skiprows=1)
    tabla = np.transpose(tabla)
    [listx_pole_S,listy_re_pole_S,listy_im_pole_S] = tabla 
    
    

label1 = r'$r_{\rm p} = k_\parallel/(k_\parallel - k_{\rm p})$'    
label2 = r'$r_{\rm p} = k_{\rm p}/(k_\parallel - k_{\rm p})$'   

graph(title,labelx,r'Re{G$_{self}$} ($\mu$m$^{-3}$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#plt.plot(listx,listy_re_ana3,'.',ms = ms,color = 'blue',label = 'PP analytical 3')
#plt.plot(listx,listy_re_ana,'.',ms = ms,color = 'purple',label = 'PP analytical')
plt.plot(listx,listy_re_ana_v1,'-',ms = ms,color = 'purple',label = 'PP ana')
plt.plot(listx,listy_re_num,'--',ms = ms,color = 'darkred',label = 'full numerical')
#plt.plot(listx,listy_re_pole_aprox_v1,'.',ms = ms+1,color = 'darkred',label = 'PP num')
#plt.plot(listx,listy_re_pole_aprox_v2,'--',ms = 3,color = 'darkred',label = 'PP num 2')
if comparar_otro_hBN == 1:
    plt.plot(listx_ana_S,listy_re_ana2_S,'--',ms = ms,color = 'blue',label = 'PP ana dip')
    plt.plot(listx_num_S,listy_re_num_S,'--',ms = ms,color = 'lightseagreen',label = 'full num dip')
    labelp = labelp + '_v2'
#    plt.plot(listx_pole_S,listy_re_pole_S,'--',ms = ms,color = 'darkred',label = 'PP num sphere')

    

plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0.1,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
#    plt.yscale('log')
os.chdir(path_save)
plt.savefig( 'Re_Gself' + labelp + '.png', format='png')   


graph(title,labelx,r'Im{G$_{self}$} ($\mu$m$^{-3}$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#plt.plot(listx,listy_re_ana3,'.',ms = ms,color = 'purple',label = 'PP analytical 3')
#plt.plot(listx,listy_im_ana,'.',ms = ms,color = 'purple',label = 'PP analytical')
plt.plot(listx,listy_im_ana_v2,'-',ms = ms,color = 'purple',label = 'PP ana')
#plt.plot(listx,listy_im_ana_v2,'--',ms = ms+1,color = 'purple',label = 'PP ana 2 ' +  label2)
plt.plot(listx,listy_im_num,'--',ms = ms,color = 'darkred',label = 'full numerical')
#plt.plot(listx,listy_im_pole_aprox_v1,'.',ms = ms+1,color = 'darkred',label = 'PP num')
#plt.plot(listx,listy_im_pole_aprox_v2,'--',ms = 3,color = 'darkred',label = 'PP num 2')
if comparar_otro_hBN == 1:
    plt.plot(listx_ana_S,listy_im_ana2_S,'--',ms = ms,color = 'blue',label = 'PP ana dip')
    plt.plot(listx_num_S,listy_im_num_S,'--',ms = ms,color = 'lightseagreen',label = 'full num dip')
#    plt.plot(listx_pole_S,listy_im_pole_S,'--',ms = ms,color = 'darkred',label = 'PP num sphere')

plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0.1,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
#    plt.yscale('log')
os.chdir(path_save)
plt.savefig( 'Im_Gself' + labelp + '.png', format='png')   


#%%
