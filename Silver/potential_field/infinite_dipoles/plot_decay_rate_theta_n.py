
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

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/potential_field/infinite_dipoles/potential_and_electric_field_with_dipole_moment_formula/decay_rate_second_try','')
path_save = path_basic + '/' + 'decay_rate_theta_n'
if not os.path.exists(path_save):
    print('Creating folder to save graphs')
    os.mkdir(path_save)

err = 'decay_rate_theta_n.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from decay_rate_theta_n import decay_rate_theta_inf_dipoles_ana
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


epsi1, epsi3 = 1,1

zp = 0.05
b = - 0.01



#Ndipoles = 50
#list_xD = np.linspace(-0.01,0.01,Ndipoles)


d_nano = 1

int_v = 20
int_v = 10


a = 5
a_nm = a*1e3



energy0_pol = 43
omega0 = energy0_pol*1e-3/hb 
#R = 10 # 10 nm en unidades de micrometros
kappa_factor_omega0 = 0.1
kappa_r_factor= 0.5

E = 124.55*1e-3
omegac0 = E/(c*hb)


title2 = r'v = c/%i, b = %i nm, $z_p$ = %i nm, d = %i nm' %(int_v,b*1e3,zp*1e3,d_nano) 
title3 = r'$\hbar\omega_0$ = %i meV, a = %i nm' %(energy0_pol,a*1e3)



labelp = r'_a%inm' %(a*1e3)

labelx = r'$\theta$ [radians]'  
labely = r'$\Gamma_n$'
title = title2 + '\n' +  title3  + ', E = %.2f meV' %(E*1e3)
label1 = 'vs_theta'  + labelp

listx = np.linspace(-np.pi,np.pi,200)

#listx = np.linspace(0.05,np.pi,200)
#listy = np.linspace(-cota_x,cota_x,100)
#listx = listy

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
    
def function_real_ana(theta,Nmax):
                
    rta = decay_rate_theta_inf_dipoles_ana(omegac0,epsi1,epsi3,d_nano,int_v,zp,a,b,n,omega0,kappa_factor_omega0,kappa_r_factor,theta)
    # print(rta)
    return rta

graph(title,labelx,labely ,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)

for n in [0,1,2]:
    list_y_re = function_real_ana(listx,n)
    maxi = np.max(list_y_re)
    print(n,maxi)
    list_y_re = np.array(list_y_re)/maxi
    plt.plot(listx,list_y_re,'.-',ms = lw, label = 'n = %i'%(n))
    

plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
plt.grid(1)
plt.tight_layout()
#plt.yscale('log')
os.chdir(path_save)
plt.savefig('decay_rate' + label1 + '.png', format='png')

