
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

paper = 1

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/two_dipoles_vz','')
path_save = path_basic + '/' + 'EELS_1'
if not os.path.exists(path_save):
    print('Creating folder to save graphs')
    os.mkdir(path_save)

err = 'fieldE_direct_numerical.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from EELS_1 import EELS_1_ana, EELS_1_num
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


#omega0THz = 65
#omega0 = omega0THz*1e12 
energy0_pol1 = 45
omega01 = energy0_pol1*1e-3/hb 

energy0_pol2 = 45
omega02 = energy0_pol2*1e-3/hb 



kappa_factor_omega01 = 0.1
kappa_r_factor1 = 0.5

kappa_factor_omega02 = 0.1
kappa_r_factor2 = 0.5


x1 = -0.05

z1 = 0
z2 = 0
 
title1 = r'$\kappa_1$ = %.2f$\omega_{01}$, $\kappa_{r1}$ = %.2f$\kappa$, $\hbar\omega_{01}$=%i meV' %(kappa_factor_omega01, kappa_r_factor1, energy0_pol1)     
title2 = r'$\kappa_2$ = %.2f$\omega_{02}$, $\kappa_{r2}$ = %.2f$\kappa$, $\hbar\omega_{02}$=%i meV, v = c/%i' %(kappa_factor_omega02, kappa_r_factor2, energy0_pol2,int_v)     
#title2 = r'$\hbar\mu$ = %.2feV, $\hbar\gamma$ = %.4feV' %(hbmu,hbgama) 
#title3 = r'$z_p$=%inm, px=%i, py=%i, pz=%i' %(zp*1e3,px,py,pz)
title4 = r'$x_1$ = %i nm, zp = %i nm' %(x1*1e3,zp*1e3)


N = 100

    # z0 = 0.06*1e3
labelx = r'$\hbar\omega$ [meV]'   
labely = r'$|x_1 - x_2|$ [nm]'   

label1 = '_vs_E_3D_x1_%i' %(x1*1e3) 
listx = np.linspace(39,52,N)

listy = np.linspace(0,250,N)


listy_final = []
for x2 in listy:
    listy_final.append(np.abs(x2-x1*1e3))

title = title1  + '\n' + title2 + '\n' + title4 

def function_num3D(energy0,x2_nano):
    omegac0 = energy0*1e-3/aux 
    
    x2 = x2_nano*1e-3

    px_f  = EELS_1_num(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,zp,x1,z1,x2,z2,omega01,kappa_factor_omega01,kappa_r_factor1,omega02,kappa_factor_omega02,kappa_r_factor2)
    # print(rta)
    return np.real(px_f)

def function_ana3D(energy0,x2_nano):
    omegac0 = energy0*1e-3/aux 
    
    x2 = x2_nano*1e-3

    px_f  = EELS_1_ana(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,zp,x1,z1,x2,z2,omega01,kappa_factor_omega01,kappa_r_factor1,omega02,kappa_factor_omega02,kappa_r_factor2)
    
    
    print(px_f)
    return np.real(px_f)

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

labelz =r'$\Gamma_{dip}$ [meV]'


X, Y = np.meshgrid(listx, listy)
#f_num = np.vectorize(function_num3D)
#Z_num = f_num(X, Y)

f_ana = np.vectorize(function_ana3D)
Z_ana = f_ana(X, Y)

#%%
limits = [np.min(listx) , np.max(listx), np.min(listy_final) , np.max(listy_final)]
graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)

# im = plt.imshow(Z, extent = limits,  cmap='RdBu', interpolation='bilinear') 
im = plt.imshow(Z_ana, extent = limits, cmap=plt.cm.hot, aspect=1/20) 
cbar = plt.colorbar(im, extend='both', fraction=0.046, pad=0.04)
cbar.ax.tick_params(labelsize = tamnum)
if paper == 0:
    cbar.set_label(labelz,fontsize=tamlegend,labelpad = 1)
plt.tight_layout(pad = 0.1)
os.chdir(path_save)
plt.savefig( 'EELS' + label1 + '.png', format='png')   

#%%


