
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

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/potential_field/potential_and_electric_field_with_dipole_moment_formula','')
path_save = path_basic + '/' + 'Green_self_image'
if not os.path.exists(path_save):
    print('Creating folder to save graphs')
    os.mkdir(path_save)

err = 'fieldE_direct_numerical.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from green_self_image_integrand import green_self_num,green_self_pole_aprox
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

int_v = 10
#v = c/int_v
#omega = 0.7*1e12
cota = 75 #nanometros
epsi1,epsi3 = 1,1

zp = 0.05
energy0 = 0.175
omegac0 = energy0/aux 

d_nano = 1

title = r'$z_p$ = %inm, d = %i nm, $\hbar\omega$ = %.3f eV' %(zp*1e3,d_nano,energy0)

N = 100


    # z0 = 0.06*1e3
labelx = r'$\alpha_\parallel$'   
label1 = '_d%inm' %(d_nano) 

listx = np.linspace(0.01,400,N)


def function_num_xx_re(u):


    rtaself_x, rtaself_y, rtaself_z = green_self_num(omegac0,epsi1,epsi3,d_nano,zp,u)
    # print(rta    
    return np.real(rtaself_x)

def function_num_xx_im(u):


    rtaself_x, rtaself_y, rtaself_z = green_self_num(omegac0,epsi1,epsi3,d_nano,zp,u)
    # print(rta    
    return np.imag(rtaself_x)


def function_pole_aprox_xx_re(u):


    rtaself_x, rtaself_y, rtaself_z  = green_self_pole_aprox(omegac0,epsi1,epsi3,d_nano,zp,u)
    
    return np.real(rtaself_x)



def function_pole_aprox_xx_im(u):


    rtaself_x, rtaself_y, rtaself_z  = green_self_pole_aprox(omegac0,epsi1,epsi3,d_nano,zp,u)
    
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


listy_re_num = []
listy_re_pole_aprox = []


listy_im_num = []
listy_im_pole_aprox = []


for value in listx: 

   
    y_re_num = function_num_xx_re(value)
    y_re_pole_aprox = function_pole_aprox_xx_re(value)
    

    listy_re_num.append(y_re_num)
    listy_re_pole_aprox.append(y_re_pole_aprox)

   
    y_im_num = function_num_xx_im(value)
    y_im_pole_aprox = function_pole_aprox_xx_im(value)

    listy_im_num.append(y_im_num)
    listy_im_pole_aprox.append(y_im_pole_aprox)
    
    
    
#%%
graph(title,labelx,r'Re(G$_{self}$ inside)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#    plt.plot(listx,listy_re_ana,'.',ms = ms,color = 'purple',label = 'analytical')
plt.plot(listx,listy_re_num,'.-',ms = ms,color = 'lightseagreen',label = 'full numerical')
plt.plot(listx,listy_re_pole_aprox,'.-',ms = 3,color = 'darkred',label = 'PP numerical')
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
#    plt.yscale('log')
os.chdir(path_save)
plt.savefig( 'Re_Gself_integrand' + label1 + '.png', format='png')   


graph(title,labelx,r'Im(G$_{self}$ inside)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#    plt.plot(listx,listy_re_ana,'.',ms = ms,color = 'purple',label = 'analytical')
plt.plot(listx,listy_im_num,'.-',ms = ms,color = 'lightseagreen',label = 'full numerical')
plt.plot(listx,listy_im_pole_aprox,'.-',ms = 3,color = 'darkred',label = 'PP numerical')
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
#    plt.yscale('log')
os.chdir(path_save)
plt.savefig( 'Im_Gself_integrand' + label1 + '.png', format='png')   


#%%
