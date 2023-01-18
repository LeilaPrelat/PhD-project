
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
path_save = path_basic + '/' + 'Green_self_image_integrand'
if not os.path.exists(path_save):
    print('Creating folder to save graphs')
    os.mkdir(path_save)

err = 'fieldE_direct_numerical.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from green_self_image_integrand import green_self_pole_aprox_integrand, green_self_pole_aprox_integrand_v2, green_self_num_integrand
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


zp_nano = 0.05*1e3

energy0 = 0.01 # eV  np.linspace(0.087, 0.8, N) ### for higher freq de 

d_nano = 0.4

title = r'$\hbar\omega$ = %.2f eV, $z_p$=%inm, d = %.2f nm' %(energy0,zp_nano,d_nano)

N = 100

    # z0 = 0.06*1e3
labelx = r'$\alpha_\parallel$'   
label1 = 'integrand_vs_E_d%inm' %(d_nano) 


list_alpha_parallel = np.linspace(0.01, 200, N)
#listx = np.linspace(0.087, 0.195, N)
#listx = np.linspace(0.2, 0.8, N)

#%%
    

def function_num_xx(alpha_parallel):
    omegac0 = energy0/aux 
    zp = zp_nano*1e-3
    rtaself_x, rtaself_y, rtaself_z = green_self_num_integrand(omegac0,epsilon_Silica,d_nano,zp,alpha_parallel)
    # print(rta    
    return rtaself_x


def function_pole_aprox_xx(alpha_parallel):
    omegac0 = energy0/aux 
    zp = zp_nano*1e-3
    rtaself_x, rtaself_y, rtaself_z  = green_self_pole_aprox_integrand(omegac0,epsilon_Silica,d_nano,zp,alpha_parallel)
    
    return rtaself_x


def function_pole_aprox_v2_xx(alpha_parallel):
    omegac0 = energy0/aux 
    zp = zp_nano*1e-3
    rtaself_x, rtaself_y, rtaself_z  = green_self_pole_aprox_integrand_v2(omegac0,epsilon_Silica,d_nano,zp,alpha_parallel)
    
    return rtaself_x


#%%
    
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
listy_re_pole_aprox2 = []


listy_im_num = []
listy_im_pole_aprox = []
listy_im_pole_aprox2 = []

for value in list_alpha_parallel: 

    y_pole = function_pole_aprox_xx(value)       
    y_pole2 = function_pole_aprox_v2_xx(value)        
    y_num = function_num_xx(value)

    listy_re_num.append(np.real(y_num))
    listy_im_num.append(np.imag(y_num))
    
    listy_re_pole_aprox.append(np.real(y_pole))
    listy_im_pole_aprox.append(np.imag(y_pole))

    listy_re_pole_aprox2.append(np.real(y_pole2))
    listy_im_pole_aprox2.append(np.imag(y_pole2))

      
#%%
graph(title,labelx,r'Re(G$_{self}$) integrand',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(list_alpha_parallel,listy_re_pole_aprox2,'.-',ms = ms+1,color = 'darkorange',label = r'$r_{\rm p} = k_{\rm p}/(k_\parallel - k_{\rm p})$')
plt.plot(list_alpha_parallel,listy_re_pole_aprox,'.-',ms = 3,color = 'darkred',label = r'$r_{\rm p} = k_\parallel/(k_\parallel - k_{\rm p})$')
plt.plot(list_alpha_parallel,listy_re_num,'.',ms = ms,color = 'lightseagreen',label = 'rp fresnel')
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0.1,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
#    plt.yscale('log')
os.chdir(path_save)
plt.savefig( 'Re_Gself' + label1 + '.png', format='png')   


graph(title,labelx,r'Im(G$_{self}$) integrand',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(list_alpha_parallel,listy_im_pole_aprox2,'.-',ms = ms+1,color = 'darkorange',label =  r'$r_{\rm p} = k_{\rm p}/(k_\parallel - k_{\rm p})$')
plt.plot(list_alpha_parallel,listy_im_pole_aprox,'.-',ms = 3,color = 'darkred',label = r'$r_{\rm p} = k_\parallel/(k_\parallel - k_{\rm p})$')
plt.plot(list_alpha_parallel,listy_im_num,'.',ms = ms,color = 'lightseagreen',label = 'rp fresnel')
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0.1,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
#    plt.yscale('log')
os.chdir(path_save)
plt.savefig( 'Im_Gself' + label1 + '.png', format='png')   


#%%
