
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: leila

graficar campo total approximacion analitica
"""

import numpy as np
import sys
import os 
import matplotlib.pyplot as plt
import seaborn as sns

#%%

save_graphs = 1 #guardar los graficos 2D del campo

plot_vs_omegaTHz = 1 # graficar vs omega THz

sns.set()    

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/EELS','')
path_save = path_basic + '/' + 'EELS'

if save_graphs==1:    
    if not os.path.exists(path_save):
        print('Creating folder to save graphs')
        os.mkdir(path_save)

err2 = 'EELS_ana.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from EELS_ana import EELS_f
except ModuleNotFoundError:
    print(err2)

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

int_v = 400
v = c/int_v
b = -0.01 #electron position in z < 0 (arriba del plano)
#omega = 0.7*1e12
cota = 2
epsi1,epsi2 = 1,1
hbmu,hbgama = 0.3,0.0001
zp = 0.05
z0 = 0 # z tiene que ser 0 para que la parte directa tenga solucion analitica
zD = 0
xD = 0.05
omega0THz = 60
omega0 = omega0THz*1e12 
R = 10 # 10 nm en unidades de micrometros
kappa_factor_omega0 = 0.1
kappa_r_factor= 0.5
E0 = omega0THz*hb*1e12

list_OmegaTHz = np.linspace(0.01,2.01,51)
px, py, pz = 1,0,0
px, py, pz = 1,0,0
    
list_OmegaTHz = np.linspace(5,101,201)
list_E = np.array(list_OmegaTHz)*1e12*hb
label1 = '_vsOmega_px%i_py%i_pz%i' %(px,py,pz)
title1 = r'$\epsilon_1$ = %i, v = c/%i $\mu$m/s, b = %inm,xD = %inm, zD = %inm' %(epsi1,int_v,b*1e3,xD*1e3,zD*1e3) 
labelx = r'E [eV]'
title2A = r'$\epsilon_2$ = %i, $\hbar\mu$ = %.2feV, $\hbar\gamma$ = %.4feV, $z_p$=%inm, z = %inm' %(epsi2,hbmu,hbgama,zp*1e3,z0*1e3) 
title3 = r'$\kappa$ = %.2f$\omega_0$, $\kappa_r$ = %.2f$\kappa$, $E_0$ =%.2feV, px = %i, py = %i, pz = %i' %(kappa_factor_omega0,kappa_r_factor,E0,px,py,pz)   

title = title1 + '\n'  + title2A + '\n' + title3   

#%%

tamfig = (4.5,3.5)
tamlegend = 10
tamletra = 11
tamtitle = 10
tamnum = 9
labelpady = -1.5
labelpadx = 0.8
pad = 0
mk = 2
ms = 4
hp = 0.3

def graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad):
    plt.figure(figsize=tamfig)
    plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
    plt.ylabel(labely,fontsize=tamletra,labelpad =labelpadx)
    plt.tick_params(labelsize = tamnum, pad = pad)
    plt.title(title,fontsize=int(tamtitle*0.9))

    return   

color_rgb1 = (1.0, 0.01, 0.24)
color_rgb2 = (0.8, 0.58, 0.46)
color_rgb3 = (0.55, 0.71, 0)
color_rgb4 = (0.6, 0.73, 0.89)

#%%

list_v_re = []
list_v_im = []

for omegaTHz in list_OmegaTHz: 
    
    omegaTHz = np.round(omegaTHz,8)
    omegac = omegaTHz*aux2

    rta = EELS_f(omegac,epsi1,epsi2,hbmu,hbgama,z0,xD,zD,b,zp,int_v,omega0,kappa_factor_omega0,kappa_r_factor,px,py,pz)

    list_v_re.append(rta.real)
    list_v_im.append(rta.imag)


#%%
   
print('Graficar el External field tensor')
    
#graph(title,labelx,r'Re$(E_{x})$',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#plt.plot(list_E,list_v_re,'.',ms = ms,color = 'lightseagreen')# plt.yscale('log')
##plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
#if save_graphs==1:
#    plt.tight_layout()
#    os.chdir(path_save)
#    plt.savefig( 'ReEELS' + label1 + '.png', format='png')

graph(title,labelx,r'EELS',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(list_E,list_v_im,'.',ms = ms,color = 'lightseagreen')
# plt.yscale('log')
#plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
if save_graphs==1:
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'ImEELS_2' + label1+ '.png', format='png')

#%%