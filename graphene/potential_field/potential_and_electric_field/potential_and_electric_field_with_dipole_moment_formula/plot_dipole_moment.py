
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
path_constants =  path_basic.replace('/potential_field/potential/potential_and_electric_field_with_dipole_moment_formula','')
path_save = path_basic + '/' + 'dipole_moment'
if not os.path.exists(path_save):
    print('Creating folder to save graphs')
    os.mkdir(path_save)

err = 'fieldE_direct_numerical.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from dipole_moment import dipole_moment_anav2, dipole_moment_num, dipole_moment_pole_aprox
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
energy0_pol = 80
omega0 = energy0_pol*1e-3/hb 
#R = 10 # 10 nm en unidades de micrometros
kappa_factor_omega0 = 0.1
kappa_r_factor= 0.5
 
title1 = r'v = c/%i, $\kappa$ = %.2f$\omega_0$, $\kappa_r$ = %.2f$\kappa$, $E_0$=%i meV' %(int_v, kappa_factor_omega0, kappa_r_factor, energy0_pol)     
#title2 = r'$\hbar\mu$ = %.2feV, $\hbar\gamma$ = %.4feV' %(hbmu,hbgama) 
#title3 = r'$z_p$=%inm, px=%i, py=%i, pz=%i' %(zp*1e3,px,py,pz)
title4 = r'$z_p$=%i nm, b = %i nm' %(zp*1e3,b*1e3)
labelp = r'_E0_%i' %(energy0_pol)

N = 100

    # z0 = 0.06*1e3
labelx = r'$\hbar\omega$ [meV]'   
label1 = 'vs_E' + labelp
listx = np.linspace(1,50,N)

title = title1  + '\n'  + title4 

def function_num(energy0):
    omegac0 = energy0*1e-3/aux 

    px_f,py_f,pz_f = dipole_moment_num(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,b,zp,omega0,kappa_factor_omega0,kappa_r_factor)
    # print(rta)
    rta = np.abs(px_f)**2 + np.abs(py_f)**2 + np.abs(pz_f)**2 

    return np.sqrt(rta)


#def function_ana(energy0):
#    omegac0 = energy0*1e-3/aux 
#
#    px_f,py_f,pz_f  = dipole_moment_ana(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,px,py,pz,b,zp,omega0,kappa_factor_omega0,kappa_r_factor)
#    
#    rta = np.abs(px_f)**2 + np.abs(py_f)**2 + np.abs(pz_f)**2 
#    
#    return np.sqrt(rta)


def function_anav2(energy0):
    omegac0 = energy0*1e-3/aux 

    px_f,py_f,pz_f  = dipole_moment_anav2(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,b,zp,omega0,kappa_factor_omega0,kappa_r_factor)
    
    rta = np.abs(px_f)**2 + np.abs(py_f)**2 + np.abs(pz_f)**2 
    
    return np.sqrt(rta)


def function_pole_approx(energy0):
    omegac0 = energy0*1e-3/aux 

    px_f,py_f,pz_f  = dipole_moment_pole_aprox(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,b,zp,omega0,kappa_factor_omega0,kappa_r_factor)
    
    rta = np.abs(px_f)**2 + np.abs(py_f)**2 + np.abs(pz_f)**2 
    
    return np.sqrt(rta)


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
listy_re_anav2 = []
listy_re_num = []
listy_re_pole = []

for value in listx: 

#    y_re_ana = function_ana(value)       
    y_re_anav2 = function_anav2(value)  
    
    y_re_num = function_num(value)
    y_re_pole = function_pole_approx(value)
    
#    listy_re_ana.append(y_re_ana)
    
    listy_re_anav2.append(y_re_anav2)

    listy_re_num.append(y_re_num)
    listy_re_pole.append(y_re_pole)

Emax = listx[np.argmax(listy_re_num)]
print(Emax)

#%%
graph(title,labelx,r'$|p|$/e',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#    plt.plot(listx,listy_re_ana,'.',ms = ms,color = 'purple',label = 'analytical')
plt.plot(listx,listy_re_anav2,'.',ms = ms,color = 'purple',label = 'PP analytical')
plt.plot(listx,listy_re_num,'.',ms = ms,color = 'lightseagreen',label = 'full numerical')
plt.plot(listx,listy_re_pole,'.-',ms = 3,color = 'darkred',label = 'PP numerical')
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
plt.yscale('log')
os.chdir(path_save)
plt.savefig( 'p_tot' + label1 + '.png', format='png')   

#%%
