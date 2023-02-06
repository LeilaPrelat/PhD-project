
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
plot_num = 0
plot_color_map = 0

plot_vs_zp = 0
plot_vs_freq = 1

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/potential_field/infinite_dipoles','')
path_save = path_basic + '/' + 'decay_rate_theta_n'
if not os.path.exists(path_save):
    print('Creating folder to save graphs')
    os.mkdir(path_save)

err = 'decay_rate_theta_n.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from decay_rate_theta_n import decay_rate_theta_inf_dipoles_ana_res_div_gamma0_v3
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

epsi1, epsi3 = 1,1

print('Definir parametros del problema')


b = - 0.01

d_nano = 1

int_v = 10

Nmax = 5

labely = r'$\Gamma_n$'


if plot_vs_zp == 1:


    theta_degree = 60 ##
    theta = theta_degree*np.pi/180
    E =1.5  # eV

    omegac0 = E/(c*hb)
    lambda_SP = 2*np.pi/omegac0 

    a_min = np.real(lambda_SP)*Nmax/(int_v - 1)
    a_max = np.real(lambda_SP)*Nmax/(int_v + 1)


    a = 0.05
    a = np.mean([a_min,a_max])
    a_nm = a*1e3



    listx = np.linspace(0.001,12,200) ### nanometros 
    listx_2 = np.array(listx)*1e-3/np.real(lambda_SP)
    
    labelx = r'$z_p$/$\lambda_p$'  



    title2 = r'v = c/%i, b = %i nm, $\theta$ = %.2fº' %(int_v,b*1e3,theta_degree) 
    title3 = r'a = %i nm' %(a*1e3)
    title4 = r', $\hbar\omega$ = %.2f eV' %(E)


    labelp = r'_a%inm_E%ieV_theta%i' %(a*1e3,E,theta_degree)
    label1 = 'vs_zp'  + labelp


elif plot_vs_freq == 1:    
    

    
    theta_degree = 0
    theta = theta_degree*np.pi/180
    zp = 0.05*1e3


#    omegac0 = E/(c*hb)
#    lambda_SP = 2*np.pi/omegac0
#
#    a_min = np.real(lambda_SP)*Nmax/(int_v - 1)
#    a_max = np.real(lambda_SP)*Nmax/(int_v + 1)


    a = 500*1e-3
    
    a_nm = a*1e3



    listx = np.linspace(0.5,3,200)

    listx_3 = np.array(listx)/(c*hb)

    listx_2 = zp/(2*np.pi/listx_3)
    
    labelx = r'$z_p$/$\lambda_p$'  



    title2 = r'v = c/%i, b = %i nm, $\theta$ = %.2fº' %(int_v,b*1e3,theta_degree) 
    title3 = r'a = %i nm' %(a*1e3)
    title4 = r', $z_p$ = %i nm' %(zp)


    labelp = r'_a%inm_zp%inm_theta%i' %(a*1e3,zp,theta_degree)
    label1 = 'vs_zp'  + labelp

    
    
    
title = title2 + '\n' +  title3  + title4

#%%


def function_real_ana(zp_nano,energy0,Nmax):
                
    omegac0 = energy0/(c*hb)  
    zp = zp_nano*1e-3
         
    rta = decay_rate_theta_inf_dipoles_ana_res_div_gamma0_v3(omegac0,epsi1,epsi3,d_nano,int_v,zp,a,b,n)

    return rta


tamfig = (4.5,3.5)
tamlegend = 12
tamletra = 12
tamtitle = 11
tamnum = 10
labelpady = -1.5
labelpadx = 2
pad = 0
mk = 2
lw = 1.5
ms = 4
hp = 0.3
length_marker = 1


#%%
    


def graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad):
    plt.figure(figsize=tamfig)
    plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
    plt.ylabel(labely,fontsize=tamletra,labelpad =labelpadx)
    plt.tick_params(labelsize = tamnum, pad = pad)
    plt.title(title,fontsize=int(tamtitle*0.9))

    return   

#%%

    
graph(title,labelx,labely ,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
for n in [0,1,2,3,4]:
    
    list_y_re = []

    if plot_vs_zp == 1:
        for x in listx:
            list_y_re.append(function_real_ana(x,E,n))
        
    elif plot_vs_freq == 1:
        for x in listx:
            list_y_re.append(function_real_ana(zp,x,n))
        
    maxi = np.max(list_y_re)
    print(n,maxi)
    list_y_re = np.array(list_y_re)/maxi
    
    plt.plot(listx_2,list_y_re,'.-',ms = lw, label = 'n = %i'%(n))
    
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
plt.grid(1)
plt.tight_layout()
#plt.yscale('log')
os.chdir(path_save)
plt.savefig('decay_rate_res_' + label1 + '.png', format='png')



#%%

