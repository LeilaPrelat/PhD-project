
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

plot_vs_E = 0
plot_vs_c = 1

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/potential_field/potential_and_electric_field/potential_and_electric_field_simple_version','')
path_save = path_basic + '/' + 'EELS'
if not os.path.exists(path_save):
    print('Creating folder to save graphs')
    os.mkdir(path_save)

err = 'EELS_film.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from EELS_1dip import EELS_ana_f,EELS_num_f,EELS_ana_f_dir
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
cota = 75 #nanometros
epsi1,epsi2 = 1,1
hbmu,hbgama = 0.3,0.0001
zp = 0.05
b = -0.01


px = 1
py = 0
pz = 0

title2 = r'$\hbar\mu$ = %.2feV, $\hbar\gamma$ = %.4feV' %(hbmu,hbgama) 
title4 = r'$z_p$=%inm, px = %i, py = %i, pz = %i' %(zp*1e3,px,py,pz)
labelp = r'_px%i_py%i_pz%i' %(px,py,pz)


N = 200


if plot_vs_E == 1:

    int_v = 10    
    labelx = r'$\hbar\omega$ [meV]'
    
    title4 = title4 + ', v = c/%i' %(int_v)
    
    label1 = 'vs_E' + labelp
    #    listx = np.linspace(0.0001,2,N)
    listx = np.linspace(15,65,N)

else:
    E0 = 43
    labelx = r'v/c'
    
    title4 = title4 + ', $\hbar\omega$ = %imeV' %(E0)
    
    label1 = 'vs_c' + labelp
    #    listx = np.linspace(0.0001,2,N)
    listx = np.linspace(0.05,0.99,N)

title = title2 + '\n' + title4

#%%

def function_ana_re(energy0,int_v0):
    omegac0 = energy0*1e-3/aux 

    rta = EELS_ana_f(omegac0,epsi1,epsi2,hbmu,hbgama,px,py,pz,int_v0,b,zp)
    
    return np.real(rta)

def function_num_re(energy0,int_v0):
    omegac0 = energy0*1e-3/aux 

    rta = EELS_num_f(omegac0,epsi1,epsi2,hbmu,hbgama,px,py,pz,int_v0,b,zp)
    
    return np.real(rta)

def function_ana_re_dir(energy0,int_v0):
    omegac0 = energy0*1e-3/aux 

    rta = EELS_ana_f_dir(omegac0,epsi1,epsi2,px,py,pz,int_v0,b)
    
    return np.real(rta)

def function_ana_im(energy0,int_v0):
    omegac0 = energy0*1e-3/aux 

    rta = EELS_ana_f(omegac0,epsi1,epsi2,hbmu,hbgama,px,py,pz,int_v0,b,zp)
    
    return np.imag(rta)

def function_num_im(energy0,int_v0):
    omegac0 = energy0*1e-3/aux 

    rta = EELS_num_f(omegac0,epsi1,epsi2,hbmu,hbgama,px,py,pz,int_v0,b,zp)
    
    return np.imag(rta)


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

listy_re_ana = []
listy_re_ana_dir = []
listy_re_ana_tot = []


listy_re_num = []

listy_im_ana = []
listy_im_num = []


if plot_vs_E == 1:
    for value in listx: 
    
        y_re_ana = function_ana_re(value,int_v)        
        y_re_num = function_num_re(value,int_v)        
        
        y_re_ana_dir = function_ana_re_dir(value,int_v)        
        
        listy_re_ana.append(y_re_ana)
        listy_re_num.append(y_re_num)
    
        listy_re_ana_dir.append(y_re_ana_dir)
        listy_re_ana_tot.append(y_re_ana_dir + y_re_ana)
    
    
        y_im_ana = function_ana_im(value,int_v)        
        y_im_num = function_num_im(value,int_v)        
        
        listy_im_ana.append(y_im_ana)
        listy_im_num.append(y_im_num)

else:
    
    for value in listx: 
        
            value2 = 1/value
        
            y_re_ana = function_ana_re(E0,value2)        
            y_re_num = function_num_re(E0,value2)        

            y_re_ana_dir = function_ana_re_dir(E0,value2)    
            
            listy_re_ana.append(y_re_ana)
            listy_re_num.append(y_re_num)
            
            listy_re_ana_dir.append(y_re_ana_dir)
            listy_re_ana_tot.append(y_re_ana_dir + y_re_ana)
        
            y_im_ana = function_ana_im(E0,value2)        
            y_im_num = function_num_im(E0,value2)        
            
            listy_im_ana.append(y_im_ana)
            listy_im_num.append(y_im_num)
    
    
    
#%%



graph(title,labelx,'$\Gamma_{ind}$',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(listx,listy_re_ana,'.-',ms = ms,color = 'purple')
#plt.plot(listx,listy_re_num,'.-',ms = ms,color = 'lightseagreen',label ='full numerical' )
plt.tight_layout()
#if plot_vs_c == 1:
#    plt.yscale('log')
os.chdir(path_save)
plt.savefig( 'EELS_ind' + label1 + '.png', format='png')   



graph(title,labelx,'$\Gamma_{dir}$',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(listx,listy_re_ana_dir,'.-',ms = ms,color = 'purple')
#plt.plot(listx,listy_re_num,'.-',ms = ms,color = 'lightseagreen',label ='full numerical' )
plt.tight_layout()
#if plot_vs_c == 1:
#    plt.yscale('log')
os.chdir(path_save)
plt.savefig( 'EELS_dir' + label1 + '.png', format='png')   


graph(title,labelx,'$\Gamma_{tot}$',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(listx,listy_re_ana_tot,'.-',ms = ms,color = 'purple')
#plt.plot(listx,listy_re_num,'.-',ms = ms,color = 'lightseagreen',label ='full numerical' )
plt.tight_layout()
#if plot_vs_c == 1:
#    plt.yscale('log')
os.chdir(path_save)
plt.savefig( 'EELS_tot' + label1 + '.png', format='png')   


graph(title,labelx,'$\Gamma$',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(listx,listy_re_ana,'.-',ms = ms,color = 'purple',label = 'analytical')
plt.plot(listx,listy_re_num,'.-',ms = ms,color = 'lightseagreen',label ='numerical' )
plt.legend(loc = 'best',markerscale=1.5,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
#if plot_vs_c == 1:
#    plt.yscale('log')


graph(title,labelx,'$\Gamma$',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(listx,listy_im_ana,'.-',ms = ms,color = 'purple',label = 'analytical')
plt.plot(listx,listy_im_num,'.-',ms = ms,color = 'lightseagreen',label ='numerical' )
plt.legend(loc = 'best',markerscale=1.5,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
#if plot_vs_c == 1:
#    plt.yscale('log')

    