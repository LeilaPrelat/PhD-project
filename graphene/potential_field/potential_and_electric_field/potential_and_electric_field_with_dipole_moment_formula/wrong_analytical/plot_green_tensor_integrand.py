
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
#sns.set()

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/potential_field/potential/potential_and_electric_field_with_dipole_moment_formula','')
path_save = path_basic + '/' + 'green_self_tensor_integrand'
if not os.path.exists(path_save):
    print('Creating folder to save graphs')
    os.mkdir(path_save)

err = 'green_self_tensor.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from green_self_tensor_integrand import green_tensor_self_num
except ModuleNotFoundError:
    print(err)

try:
    sys.path.insert(1, path_constants)
    from graphene_sigma import sigma_DL
except ModuleNotFoundError:
    print('graphene_sigma.py no se encuentra en ' + path_basic)

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
labely = r'$\alpha_\parallel$'
label1 = 'vs_E' + labelp
    
title = title2 + '\n' + title3 

def function_num_re(energy0,u):
    omegac0 = energy0*1e-3/aux 

    rta = green_tensor_self_num(omegac0,epsi1,epsi2,hbmu,hbgama,px,py,pz,zp,u)
    
    return np.real(rta)

def function_num_im(energy0,u):
    omegac0 = energy0*1e-3/aux 

    rta = green_tensor_self_num(omegac0,epsi1,epsi2,hbmu,hbgama,px,py,pz,zp,u)
    
    return np.imag(rta)

def function_num_abs(energy0,u):
    omegac0 = energy0*1e-3/aux 

    rta = green_tensor_self_num(omegac0,epsi1,epsi2,hbmu,hbgama,px,py,pz,zp,u)
    
    return np.abs(rta)

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

list_x = np.linspace(1,150,N)
list_y = np.linspace(1,150,N)
energy0 = 43
omega0THz = energy0*1e-12/hb
omegac0 = energy0/(c*hb)

cond = 4*np.pi*alfac*sigma_DL(energy0*1e-3,hbmu,hbgama)

Rp = 2*epsi1/(epsi1 + epsi2)
alfa_p = 1j*(epsi1 + epsi2)/(cond)
k_p = alfa_p*omegac0
print(np.exp(-alfa_p*omegac0*2*zp))

list_y1D = []
for u in list_x:
#    energy0 = 2*alfac*hbmu*u
    
    rta = function_num_abs(energy0,u)
    list_y1D.append(rta)

#%%
    
ejey = np.linspace(np.min(list_y1D),np.max(list_y1D),N)

plt.figure(figsize=tamfig)      
plt.title(title,fontsize=tamtitle)
plt.plot(list_x,list_y1D,'-',ms = 2)
plt.xlabel(r'$\alpha_{\parallel}$',fontsize=tamletra, labelpad = 0)
plt.ylabel(r'F($\alpha_{\parallel}$,k)',fontsize=tamletra, labelpad = 0)
#plt.plot(alfa_p*np.ones(N),ejey,'-')
plt.legend(loc='best',markerscale=2,fontsize=int(tamlegend*0.8))
plt.tick_params(labelsize = tamnum,pad = pad)
plt.grid(1)
plt.tight_layout()
os.chdir(path_save)
plt.savefig('integrandF.png')

#%%


X, Y = np.meshgrid(list_x, list_y)
limits = [min(list_x) , max(list_x), min(list_y) , max(list_y)]

f1_real = np.vectorize(function_num_re)
Z1_real = f1_real(X, Y)

graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
im = plt.imshow(Z1_real, extent = limits, cmap=plt.cm.hot, interpolation='bilinear')
cbar = plt.colorbar(im)
cbar.ax.tick_params(labelsize = tamnum)
# cbar.set_label('Re(Ex)',fontsize=tamlegend,labelpad = 1)

# plt.plot(listx,listy_re_ana,'.',ms = ms,color = 'purple',label = 'pole aprox')
# plt.plot(listx,listy_re_num,'.-',ms = 3,color = 'lightseagreen',label = 'sin aprox')
# plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
# plt.tight_layout()
#    plt.yscale('log')
os.chdir(path_save)
plt.savefig( 'Re_Gself' + label1 + '.png', format='png')   


graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
f1_imag = np.vectorize(function_num_im)
Z1_imag = f1_imag(X, Y)
im = plt.imshow(Z1_imag, extent = limits, cmap=plt.cm.hot, interpolation='bilinear')
cbar = plt.colorbar(im)
cbar.ax.tick_params(labelsize = tamnum)
os.chdir(path_save)
plt.savefig( 'Im_Gself' + label1 + '.png', format='png')   

#%%
