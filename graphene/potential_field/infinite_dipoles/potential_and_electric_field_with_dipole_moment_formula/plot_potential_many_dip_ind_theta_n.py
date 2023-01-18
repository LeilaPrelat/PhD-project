
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
#import seaborn as sns

plot_num = 0
plot_color_map = 0

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/potential_field/infinite_dipoles','')
path_save = path_basic + '/' + 'many_dipoles_phi_ind3D'
if not os.path.exists(path_save):
    print('Creating folder to save graphs')
    os.mkdir(path_save)

err = 'potential_many_dipoles_ind.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from potential_theta_n_many_dipoles_ind import P_theta_many_dipoles_ana_n
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


epsi1, epsi2 = 1,1
hbmu, hbgama = 0.3,0.0001

zp = 0.05
b = - 0.01

z0 = zp 
x0 = 500*1e-3 ## en micrones
y0 = 500*1e-3 ## en micrones

#Ndipoles = 50
#list_xD = np.linspace(-0.01,0.01,Ndipoles)




int_v = 20
int_v = 10


a = 0.05
a_nm = a*1e3



energy0_pol = 43
omega0 = energy0_pol*1e-3/hb 
#R = 10 # 10 nm en unidades de micrometros
kappa_factor_omega0 = 0.1
kappa_r_factor= 0.5

E = 43*1e-3
omegac0 = E/(c*hb)


title2 = r'v = c/%i, b = %i nm, $z_p$ = z =%i nm, x = %i nm' %(int_v,b*1e3,zp*1e3,x0*1e3) 
title3 = r'y = %i nm, $\hbar\omega_0$ = %i meV, a = %i nm' %(y0*1e3,energy0_pol,a*1e3)



labelp = r'_a%inm' %(a*1e3)

labelx = r'$\theta$ [radians]'  
labely = r'$P_n$'
title = title2 + '\n' +  title3  + ', E = %.2f meV' %(E*1e3)
label1 = 'vs_theta'  + labelp

listx = np.linspace(-2*np.pi,2*np.pi,200)


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
                
    rta = P_theta_many_dipoles_ana_n(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,zp,x0,y0,z0,a,b,Nmax,omega0,kappa_factor_omega0,kappa_r_factor,theta)
    # print(rta)
    return rta.real

graph(title,labelx,labely + ' normalized [a.u]' ,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)

for n in [-15,1,50]:
    list_y_re = function_real_ana(listx,n)
    maxi = np.max(list_y_re)
    print(n,maxi)
    list_y_re = np.array(list_y_re)/np.abs(maxi)
    plt.plot(listx,list_y_re,'.-',ms = lw, label = 'n = %i'%(n))
    

plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
plt.grid(1)
plt.tight_layout()
os.chdir(path_save)
plt.savefig('disp_relation_for_potential' + label1 + '.png', format='png')


