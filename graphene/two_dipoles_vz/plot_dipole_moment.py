
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
path_constants =  path_basic.replace('/two_dipoles_vz','')
path_save = path_basic + '/' + 'dipole_moment'
if not os.path.exists(path_save):
    print('Creating folder to save graphs')
    os.mkdir(path_save)

err = 'fieldE_direct_numerical.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from dipole_moment_z import dipole_moment_z_ana, dipole_moment_z_num
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
epsi1,epsi2 = 1,1
hbmu,hbgama = 0.3,0.0001
zp = 0.05
b = -0.01



#omega0THz = 65
#omega0 = omega0THz*1e12 
energy0_pol1 = 35
omega01 = energy0_pol1*1e-3/hb 

energy0_pol2 = 55
omega02 = energy0_pol2*1e-3/hb 



kappa_factor_omega01 = 0.1
kappa_r_factor1 = 0.5

kappa_factor_omega02 = 0.1
kappa_r_factor2 = 0.5


x1 = -0.05
x2 = 0.05

z1 = 0
z2 = 0
 
title1 = r'$\hbar\omega_{01}$=%i meV, $\hbar\omega_{02}$=%i meV' %(energy0_pol1, energy0_pol2)     
#title2 = r'$\kappa_2$ = %.2f$\omega_{02}$, $\kappa_{r2}$ = %.2f$\kappa$, $\hbar\omega_{02}$=%i meV' %(kappa_factor_omega02, kappa_r_factor2, energy0_pol2)     
#title2 = r'$\hbar\mu$ = %.2feV, $\hbar\gamma$ = %.4feV' %(hbmu,hbgama) 
#title3 = r'$z_p$=%inm, px=%i, py=%i, pz=%i' %(zp*1e3,px,py,pz)
title4 = r'v = c/%i, $x_1$ = %i nm, $x_2$ = %i nm' %(int_v,x1*1e3,x2*1e3)


N = 100

    # z0 = 0.06*1e3
labelx = r'$\hbar\omega$ [meV]'   
label1 = '_vs_E_x1_%i_x2_%i' %(x1*1e3,x2*1e3) 
listx = np.linspace(20,70,N)

title = title1   + '\n' + title4 

def function_num(energy0):
    omegac0 = energy0*1e-3/aux 

    px_f  = dipole_moment_z_num(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,zp,x1,z1,x2,z2,omega01,kappa_factor_omega01,kappa_r_factor1,omega02,kappa_factor_omega02,kappa_r_factor2)
    # print(rta)
    return px_f

def function_ana(energy0):
    omegac0 = energy0*1e-3/aux 

    px_f  = dipole_moment_z_ana(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,zp,x1,z1,x2,z2,omega01,kappa_factor_omega01,kappa_r_factor1,omega02,kappa_factor_omega02,kappa_r_factor2)
       
    return px_f

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
listy_re_num = []

listy_im_ana = []
listy_im_num = []


for value in listx: 

#    y_re_ana = function_ana(value)       
    y_ana = function_ana(value)  
    
    y_num = function_num(value)
 
#    listy_re_ana.append(y_re_ana)
    
    listy_re_ana.append(np.real(y_ana))
    listy_re_num.append(np.real(y_num))

    listy_im_ana.append(np.imag(y_ana))
    listy_im_num.append(np.imag(y_num))


Emax = listx[np.argmax(listy_re_num)]
print(Emax)

#%%
graph(title,labelx,r'Re$\{p_z\}$/e',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#    plt.plot(listx,listy_re_ana,'.',ms = ms,color = 'purple',label = 'analytical')
plt.plot(listx,listy_re_ana,'.-',ms = ms,color = 'purple',label = 'analytical')
plt.plot(listx,listy_re_num,'.-',ms = ms,color = 'lightseagreen',label = 'numerical')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
#plt.yscale('log')
os.chdir(path_save)
plt.savefig( 'p_re' + label1 + '.png', format='png')   



graph(title,labelx,r'Im$\{p_z\}$/e',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#    plt.plot(listx,listy_re_ana,'.',ms = ms,color = 'purple',label = 'analytical')
plt.plot(listx,listy_im_ana,'.-',ms = ms,color = 'purple',label = 'analytical')
plt.plot(listx,listy_im_num,'.-',ms = ms,color = 'lightseagreen',label = 'numerical')
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
#plt.yscale('log')
os.chdir(path_save)
plt.savefig( 'p_im' + label1 + '.png', format='png')  

#%%
