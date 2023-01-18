
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
path_constants =  path_basic.replace('/potential_field/infinite_dipoles/num_vs_ana','')
path_save = path_basic + '/' + 'term1_phi_integrand'
if not os.path.exists(path_save):
    print('Creating folder to save graphs')
    os.mkdir(path_save)

err = 'green_self_tensor.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from term_integrand_phi_many_dipole import phi_many_dipoles_integrand_term1
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

x =  0.00
y =  0.00
int_v = 10
px = 1
a = 0.5 #micrones 
Ndip = 50

#   minimum frequency to have plasmons
n = Ndip/2 - 1
omega_p = 2*n*np.pi*c/(a*(int_v - 1))
energy_p = omega_p*hb # energy in electron volts

energy0 = 43
omega0THz = energy0*1e-12/hb
omegac0 = energy0/(c*hb)
cond = 4*np.pi*alfac*sigma_DL(energy0*1e-3,hbmu,hbgama)


print(energy_p < energy0)

Rp = 2*epsi1/(epsi1 + epsi2)
alfa_p = 1j*(epsi1 + epsi2)/(cond)
print('we can neglect the imaginary part of alpha_p', np.imag(alfa_p)/np.real(alfa_p)<1e-2)

alfa_x = int_v - 2*pi*n/(omegac0*a)
print(np.sqrt(alfa_p**2 - alfa_x**2))

title2 = r'E = %imeV, $\hbar\mu$ = %.2feV, $\hbar\gamma$ = %.4feV, v = c/%i' %(energy0,hbmu,hbgama,int_v) 
title3 = r'px = %i, a = %inm, N = %i, $z_p$=%inm, x=%inm, y=%inm' %(px,a*1e3, Ndip*2+1, zp*1e3,x*1e3,y*1e3)


# z0 = 0.06*1e3
labelx = 'E [meV]'   
labely = r'$\alpha_\parallel$'
label1 = 'vs_E' 
    
title = title2 + '\n' + title3 

def function_num_re(energy0,alpha_y): # energy entra en unidades de meV
    omegac0 = energy0*1e-3/aux 

    rta = phi_many_dipoles_integrand_term1(omegac0,epsi1,epsi2,hbmu,hbgama,x,y,a,Ndip,zp,int_v,px,alpha_y)
    
    return np.real(rta)

def function_num_im(energy0,alpha_y): # energy entra en unidades de meV
    omegac0 = energy0*1e-3/aux 

    rta = phi_many_dipoles_integrand_term1(omegac0,epsi1,epsi2,hbmu,hbgama,x,y,a,Ndip,zp,int_v,px,alpha_y)
    
    return np.imag(rta)

def function_num_abs(energy0,alpha_y): # energy entra en unidades de meV
    omegac0 = energy0*1e-3/aux 

    rta = phi_many_dipoles_integrand_term1(omegac0,epsi1,epsi2,hbmu,hbgama,x,y,a,Ndip,zp,int_v,px,alpha_y)
    
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
Narray = 100
list_x = np.linspace(0.1,150,Narray)
list_y = np.linspace(0.1,150,Narray)


list_y1D = []
for u in list_x:
#    energy0 = 2*alfac*hbmu*u
    
    rta = function_num_abs(energy0,u)
    list_y1D.append(rta)

#%%
    
ejey = np.linspace(np.min(list_y1D),np.max(list_y1D),Narray)

plt.figure(figsize=tamfig)      
plt.title(title,fontsize=tamtitle)
plt.plot(list_x,list_y1D,'-',ms = 2)
plt.xlabel(r'$\alpha_{y}=k_y/k$',fontsize=tamletra, labelpad = 0)
plt.ylabel(r'|M$_1$($\alpha_{y}$,k)|',fontsize=tamletra, labelpad = 0)
#plt.plot(alfa_p*np.ones(N),ejey,'-')
plt.legend(loc='best',markerscale=2,fontsize=int(tamlegend*0.8))
plt.tick_params(labelsize = tamnum,pad = pad)
plt.grid(1)
plt.tight_layout()
os.chdir(path_save)
plt.savefig('integrandM1_E0%i.png' %(energy0))

#%%



#%%
