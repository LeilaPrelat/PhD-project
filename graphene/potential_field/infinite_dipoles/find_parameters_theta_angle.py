
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

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/potential_field/infinite_dipoles','')
path_save = path_basic + '/' + 'theta_ang'
if not os.path.exists(path_save):
    print('Creating folder to save graphs')
    os.mkdir(path_save)

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

#Ndipoles = 50
#list_xD = np.linspace(-0.01,0.01,Ndipoles)

def omega_n_THz(int_v,a_nm,Nmax):  ## omega puede valer esto o menos para que se generen plasmones
    
    a = a_nm*1e-3
    
    aux = alfac*int_v*hbmu/(hb)
    rta = aux + 0.5*np.sqrt((2*aux)**2 + 16*alfac*c*hbmu*Nmax*np.pi/(a*hb))
    
    return rta



def theta(int_v,a_nm,Nmax):
    

    
    omega_n = omega_n_THz(int_v,a_nm,Nmax)
    omegac = omega_n/c
#    print('omega/c:', omegac)
    
 #   omegac = omegac - 0.5 ## ver si funciona con una freq menor 

    E = omegac*aux    
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    alfa_p = 1j*(epsi1 + epsi2)/(cond)
    kp = alfa_p*omegac    
    
    lambdda_p = 2*pi/kp
    
#    
#    den1 = omegac*int_v/(2*np.pi) - 1/lambdda_p
#    den2 = omegac*int_v/(2*np.pi) + 1/lambdda_p
#    
#    a_min = Nmax/den1
#    
#    a_max = Nmax/den2
#    
#    a = np.mean([a_min,a_max])
#    print('a:', a)
    a = a_nm*1e-3
    
    theta0 = np.arccos(lambdda_p*(omegac*int_v/(2*np.pi) - Nmax/a))
    
    return theta0 
    
#%%
    
tamfig = (4.5,3.5)
tamlegend = 11
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

Nnum = 100
labelx = 'v/c'
labely = 'a [nm]' 


list_n = [1,2,3,4,5]

for Nmax in list_n : 
    
    label1 = '_n%i' %(Nmax) 

    title = r'$E_F$ = %.2f eV, n = %i' %(hbmu,Nmax)
    
    def theta_n(int_v_div,a_nm):
        int_v = 1/int_v_div
        theta0 = theta(int_v,a_nm,Nmax)
        
        if np.isnan(theta0) == False:
            
            return np.real(theta0)
    

    listx = np.linspace(0.001,1,Nnum)
    listy = np.linspace(10,5000,Nnum)
    
    
    N = 100
    X, Y = np.meshgrid(listx, listy)
    limits = [min(listx) , max(listx), min(listy) , max(listy)]
    
    long_x = max(listx) - min(listx) 
    long_y = max(listy) - min(listy) 

    
    f_real_ana = np.vectorize(theta_n)
    Z_real_ana = f_real_ana(X, Y)
    
    
    graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    vmin, vmax = np.min(Z_real_ana), np.max(Z_real_ana)
    im = plt.imshow(Z_real_ana, extent = limits, cmap=plt.cm.hot, aspect=long_x/long_y)
    cbar = plt.colorbar(im)
    cbar.ax.tick_params(labelsize = tamnum)
    cbar.set_label(r'$\theta_n$')
    os.chdir(path_save)
    plt.tight_layout()
    plt.savefig('theta_n' + label1 + '.png', format='png')

            
#%%         
        
    