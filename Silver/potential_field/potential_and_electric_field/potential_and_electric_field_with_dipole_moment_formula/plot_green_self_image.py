
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

plot_vs_E = 1
plot_vs_zp = 0

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/potential_field/potential/potential_and_electric_field_with_dipole_moment_formula','')
path_save = path_basic + '/' + 'Green_self_image'
if not os.path.exists(path_save):
    print('Creating folder to save graphs')
    os.mkdir(path_save)

err = 'fieldE_direct_numerical.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from green_self_image import green_self_ana,green_self_pole_aprox, green_self_ana2, green_self_num,green_self_ana3
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
epsi1,epsi3 = 1,1
d_nano = 1

if plot_vs_E == 1:
    zp0_nano = 0.25*1e3
    listx = np.linspace(0.001,3.5,100)
    title = r'$z_p$ = %inm, d = %i nm' %(zp0_nano,d_nano)
    labelx = r'$\hbar\omega$ [eV]'   
    label1 = '_vs_E_d%inm' %(d_nano) 


elif plot_vs_zp == 1:
    E0 = 3
    listx = np.linspace(10,400,100)
    title = r'$\hbar\omega$ = %i eV, d = %i nm' %(E0,d_nano)
    labelx = r'$z_p$ [nm]'   
    label1 = '_vs_zp_d%inm' %(d_nano) 


N = 100

    # z0 = 0.06*1e3


def function_num_xx_re(energy0,zp_nano):
    zp = zp_nano*1e-3

    omegac0 = energy0/aux 

    rtaself_x, rtaself_y, rtaself_z = green_self_num(omegac0,epsi1,epsi3,d_nano,zp)
    # print(rta    
    return np.real(rtaself_x)

def function_num_xx_im(energy0,zp_nano):
    zp = zp_nano*1e-3

    omegac0 = energy0/aux  

    rtaself_x, rtaself_y, rtaself_z = green_self_num(omegac0,epsi1,epsi3,d_nano,zp)
    # print(rta    
    return np.imag(rtaself_x)


def function_ana_xx_re(energy0,zp_nano):
    zp = zp_nano*1e-3

    omegac0 = energy0/aux 

    rtaself_x, rtaself_y, rtaself_z  = green_self_ana(omegac0,epsi1,epsi3,d_nano,zp)
    
    return np.real(rtaself_x)


def function_ana2_xx_re(energy0,zp_nano):
    zp = zp_nano*1e-3

    omegac0 = energy0/aux 

    rtaself_x, rtaself_y, rtaself_z  = green_self_ana2(omegac0,epsi1,epsi3,d_nano,zp)
    
    return np.real(rtaself_x)




def function_ana3_xx_re(energy0,zp_nano):
    zp = zp_nano*1e-3

    omegac0 = energy0/aux 

    rtaself_x, rtaself_y, rtaself_z  = green_self_ana3(omegac0,epsi1,epsi3,d_nano,zp)
    
    return np.real(rtaself_x)



def function_ana_xx_im(energy0,zp_nano):
    zp = zp_nano*1e-3

    omegac0 = energy0/aux 

    rtaself_x, rtaself_y, rtaself_z  = green_self_ana(omegac0,epsi1,epsi3,d_nano,zp)
    
    return np.imag(rtaself_x)



def function_ana2_xx_im(energy0,zp_nano):
    zp = zp_nano*1e-3

    omegac0 = energy0/aux 

    rtaself_x, rtaself_y, rtaself_z  = green_self_ana2(omegac0,epsi1,epsi3,d_nano,zp)
    
    return np.imag(rtaself_x)


def function_ana3_xx_im(energy0,zp_nano):
    zp = zp_nano*1e-3

    omegac0 = energy0/aux 

    rtaself_x, rtaself_y, rtaself_z  = green_self_ana3(omegac0,epsi1,epsi3,d_nano,zp)
    
    return np.imag(rtaself_x)

def function_pole_aprox_xx_re(energy0,zp_nano):
    zp = zp_nano*1e-3

    omegac0 = energy0/aux 

    rtaself_x, rtaself_y, rtaself_z  = green_self_pole_aprox(omegac0,epsi1,epsi3,d_nano,zp)
    
    return np.real(rtaself_x)



def function_pole_aprox_xx_im(energy0,zp_nano):
    
    zp = zp_nano*1e-3
    omegac0 = energy0/aux 

    rtaself_x, rtaself_y, rtaself_z  = green_self_pole_aprox(omegac0,epsi1,epsi3,d_nano,zp)
    
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

listy_re_ana = []
listy_re_ana2 = []
listy_re_ana3 = []
listy_re_num = []
listy_re_pole_aprox = []


listy_im_ana = []
listy_im_ana2 = []
listy_im_ana3 = []
listy_im_num = []
listy_im_pole_aprox = []

if plot_vs_E == 1:

    for value in listx: 
    
        y_re_ana = function_ana_xx_re(value,zp0_nano)       
        y_re_ana2 = function_ana2_xx_re(value,zp0_nano)       
        y_re_ana3 = function_ana3_xx_re(value,zp0_nano)  
        y_re_num = function_num_xx_re(value,zp0_nano)
        y_re_pole_aprox = function_pole_aprox_xx_re(value,zp0_nano)
        
        listy_re_ana.append(y_re_ana)
        listy_re_ana2.append(y_re_ana2)
        listy_re_ana3.append(y_re_ana3)
        listy_re_num.append(y_re_num)
        listy_re_pole_aprox.append(y_re_pole_aprox)
    
        y_im_ana = function_ana_xx_im(value,zp0_nano)       
        y_im_ana2 = function_ana2_xx_im(value,zp0_nano)   
        y_im_ana3 = function_ana3_xx_im(value,zp0_nano)  
        y_im_num = function_num_xx_im(value,zp0_nano)
        y_im_pole_aprox = function_pole_aprox_xx_im(value,zp0_nano)
    
        listy_im_ana.append(y_im_ana)
        listy_im_ana2.append(y_im_ana2)
        listy_im_ana3.append(y_im_ana3)
        listy_im_num.append(y_im_num)
        listy_im_pole_aprox.append(y_im_pole_aprox)
        
    
elif plot_vs_zp == 1:
    

    for value in listx: 
    
        y_re_ana = function_ana_xx_re(E0,value)       
        y_re_ana2 = function_ana2_xx_re(E0,value)       
        y_re_ana3 = function_ana3_xx_re(E0,value) 
        y_re_num = function_num_xx_re(E0,value)
        y_re_pole_aprox = function_pole_aprox_xx_re(E0,value)
        
        listy_re_ana.append(y_re_ana)
        listy_re_ana2.append(y_re_ana2)
        listy_re_ana3.append(y_re_ana3)
        listy_re_num.append(y_re_num)
        listy_re_pole_aprox.append(y_re_pole_aprox)
    
        y_im_ana = function_ana_xx_im(E0,value)       
        y_im_ana2 = function_ana2_xx_im(E0,value)   
        y_im_ana3 = function_ana3_xx_im(E0,value)  
        y_im_num = function_num_xx_im(E0,value)
        y_im_pole_aprox = function_pole_aprox_xx_im(E0,value)
    
        listy_im_ana.append(y_im_ana)
        listy_im_ana2.append(y_im_ana2)
        listy_im_ana3.append(y_im_ana3)
        listy_im_num.append(y_im_num)
        listy_im_pole_aprox.append(y_im_pole_aprox)
        

#%%
graph(title,labelx,r'Re(G$_{self}$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#plt.plot(listx,listy_re_ana3,'.-',ms = ms,color = 'blue',label = 'PP analytical 3')
#plt.plot(listx,listy_re_ana,'.',ms = ms,color = 'purple',label = 'PP analytical')
plt.plot(listx,listy_re_ana2,'.',ms = ms+1,color = 'darkorange',label = 'PP analytical 2')
plt.plot(listx,listy_re_num,'.',ms = ms,color = 'lightseagreen',label = 'full numerical')
plt.plot(listx,listy_re_pole_aprox,'.-',ms = 3,color = 'darkred',label = 'PP numerical')
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0.1,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
#plt.yscale('log')
os.chdir(path_save)
plt.savefig( 'Re_Gself' + label1 + '.png', format='png')   


graph(title,labelx,r'Im(G$_{self}$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#plt.plot(listx,listy_im_ana3,'.-',ms = ms,color = 'blue',label = 'PP analytical 3')
#plt.plot(listx,listy_im_ana,'.',ms = ms,color = 'purple',label = 'PP analytical')
plt.plot(listx,listy_im_ana2,'.',ms = ms+1,color = 'darkorange',label = 'PP analytical 2')
plt.plot(listx,listy_im_num,'.',ms = ms,color = 'lightseagreen',label = 'full numerical')
plt.plot(listx,listy_im_pole_aprox,'.-',ms = 3,color = 'darkred',label = 'PP numerical')
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0.1,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
#plt.yscale('log')
os.chdir(path_save)
plt.savefig( 'Im_Gself' + label1 + '.png', format='png')   


#%%
