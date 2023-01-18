
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
path_save = path_basic + '/' + 'green_self_tensor'
if not os.path.exists(path_save):
    print('Creating folder to save graphs')
    os.mkdir(path_save)

err = 'green_self_tensor.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from green_self_tensor import green_tensor_self_ana, green_tensor_self_num
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
cota = 75 #nanometros
epsi1,epsi2 = 1,1
hbmu,hbgama = 0.3,0.0001
zp = 0.05
b = -0.01

y0,x0 = 0,0 

xD = 0.005
yD = 0.005
zD = 0

px = 1
py = 0
pz = 0
  
title2 = r'$\hbar\mu$ = %.2feV, $\hbar\gamma$ = %.4feV' %(hbmu,hbgama) 
title3 = r'$z_p$=%inm, px=%i, py=%i, pz=%i' %(zp*1e3,px,py,pz)
labelp = r'_px%i_py%i_pz%i' %(px,py,pz)

N = 100
# z0 = 0.06*1e3
labelx = 'E [meV]'   
label1 = 'vs_E' + labelp
listx = np.linspace(1,50,N)
    
title = title2 + '\n' + title3 

def function_num(energy0):
    omegac0 = energy0*1e-3/aux 

    rta = green_tensor_self_num(omegac0,epsi1,epsi2,hbmu,hbgama,px,py,pz,zp)
    
    return rta


def function_ana(energy0):
    omegac0 = energy0*1e-3/aux 

    rta  = green_tensor_self_ana(omegac0,epsi1,epsi2,hbmu,hbgama,px,py,pz,zp)
    
    return rta

#%%
    
tamfig = (4.5,3.5)
tamlegend = 10
tamletra = 11
tamtitle = 10
tamnum = 9
labelpady = -1.5
labelpadx = 2
pad = 0
mk = 2
ms = 5
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

    y_ana = function_ana(value)       
    y_num = function_num(value)
    
    listy_re_ana.append(np.real(y_ana))
    listy_re_num.append(np.real(y_num))

    listy_im_ana.append(np.imag(y_ana))
    listy_im_num.append(np.imag(y_num))


graph(title,labelx,r'Re $G_{self}$ $[\mu m]^3$',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(listx,listy_re_ana,'.',ms = ms,color = 'purple',label = 'pole aprox')
plt.plot(listx,listy_re_num,'.-',ms = 3,color = 'lightseagreen',label = 'sin aprox')
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
#    plt.yscale('log')
os.chdir(path_save)
plt.savefig( 'Re_Gself' + label1 + '.png', format='png')   


graph(title,labelx,r'Im $G_{self}$ $[\mu m]^3$',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(listx,listy_im_ana,'.',ms = ms,color = 'purple',label = 'pole aprox')
plt.plot(listx,listy_im_num,'.-',ms = 3,color = 'lightseagreen',label = 'sin aprox')
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
#    plt.yscale('log')
os.chdir(path_save)
plt.savefig( 'Im_Gself' + label1 + '.png', format='png')   

#%%
