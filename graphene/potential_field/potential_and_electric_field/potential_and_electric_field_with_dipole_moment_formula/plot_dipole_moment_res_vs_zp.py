
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
    from dipole_moment import dipole_moment_anav2_resonance, dipole_moment_num_res, dipole_moment_pole_aprox_res
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

epsi1,epsi2 = 1,1
hbmu,hbgama = 0.3,0.0001
b = -0.01


energy0 = 45 #meV
 

title4 = r'v = c/%i, $\hbar\omega$ = %i meV, b = %i nm' %(int_v, energy0,b*1e3)
labelp = r'_res' 

N = 100

    # z0 = 0.06*1e3
labelx = r'$z_0$ [nm]'   
label1 = 'vs_zp' + labelp
listx = np.array(np.linspace(200,1000,N))*1e-3
listx_2 = listx*1e3
title =  title4 

def function_num(zp):
    omegac0 = energy0*1e-3/aux 

    px_f,py_f,pz_f = dipole_moment_num_res(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,b,zp)
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


def function_anav2(zp):
    omegac0 = energy0*1e-3/aux 

    px_f,py_f,pz_f  = dipole_moment_anav2_resonance(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,b,zp)
    
    rta = np.abs(px_f)**2 + np.abs(py_f)**2 + np.abs(pz_f)**2 
    
    return np.sqrt(rta)


def function_pole_approx(zp):
    omegac0 = energy0*1e-3/aux 

    px_f,py_f,pz_f  = dipole_moment_pole_aprox_res(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,b,zp)
    
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
graph(title,labelx,r'$|\mathbf{p}|^2/(e/(2\pi v))^2$ [$\mu$m$^2$]',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#    plt.plot(listx,listy_re_ana,'.',ms = ms,color = 'purple',label = 'analytical')
plt.plot(listx_2,listy_re_anav2,'.',ms = ms,color = 'purple',label = 'PP analytical')
plt.plot(listx_2,listy_re_num,'.',ms = ms,color = 'lightseagreen',label = 'full numerical')
plt.plot(listx_2,listy_re_pole,'.-',ms = 3,color = 'darkred',label = 'PP numerical')
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
plt.yscale('log')
os.chdir(path_save)
plt.savefig( 'p_tot' + label1 + '.png', format='png')   

#%%
