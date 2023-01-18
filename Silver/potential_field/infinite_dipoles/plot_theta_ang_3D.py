
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
    from Silver_PP import Silver_lambda_p, Silver_Rp
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

epsi1,epsi3 = 1,1

#Ndipoles = 50
#list_xD = np.linspace(-0.01,0.01,Ndipoles)

d_nano = 1
int_v = 10
Nmax = 0
#def omega_n_THz(int_v,a_nm,Nmax):  ## omega puede valer esto o menos para que se generen plasmones
#    
#    a = a_nm*1e-3
#    
#    aux = alfac*int_v*hbmu/(hb)
#    rta = aux + 0.5*np.sqrt((2*aux)**2 + 16*alfac*c*hbmu*Nmax*np.pi/(a*hb))
#    
#    return rta*1.5

title = r'v = c/%i, n = %i, d = %i nm' %(int_v,Nmax,d_nano) 
labelx = r'Plasmon energy, $\hbar\omega$ (eV)'  
labely = r'Lattice period, $a$ ($\mu$m)'
def theta(E,a_micros):
    

#    
#    omega_n = omega_n_THz(int_v,a_nm,Nmax)
#    omegac = omega_n/c
#    print('omega/c:', omegac)
    
 #   omegac = omegac - 0.5 ## ver si funciona con una freq menor 

    omegac = E/(hb*c)
    d_micros = d_nano*1e-3
    Rp = Silver_Rp(E,epsi1,epsi3)
    lambda_p_v = Silver_lambda_p(E,epsi1,epsi3)*d_micros
    kp = 2*np.pi/lambda_p_v
    
#    lambdda_p = 2*pi/kp
    
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
    a = a_micros
    
    theta0 = np.arccos((omegac*int_v + 2*np.pi*Nmax/a)/np.real(kp))
#    print(theta0)
    return theta0*180/np.pi


    
#%%
    
tamfig = [2.5, 2]
tamletra = 7
tamtitle  = 8
tamnum = 6
tamlegend = 6
labelpady = 2
labelpadx = 3
pad = 2.5
mk = 1
ms = 1
hp = 0.5
length_marker = 1.5
dpi = 500

def graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad):
    plt.figure(figsize=tamfig)
    plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
    plt.ylabel(labely,fontsize=tamletra,labelpad =labelpadx)
    plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in", pad = pad)
#    plt.title(title,fontsize=int(tamtitle*0.9))

    return   

#%% 
N2 = 2000
if Nmax == 1 or Nmax == 2:
    listx = np.linspace(1.8,1.98,N2)
elif Nmax == 0:
    listx = np.linspace(2,3.5,N2)
listy = np.array(np.linspace(0.2,0.6,N2))

X, Y = np.meshgrid(listx, listy)
f_num = np.vectorize(theta)
Z_num = f_num(X, Y)

#%% 

limits = [np.min(listx) , np.max(listx), np.min(listy) , np.max(listy)]
graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#for x in [x1,x2,x3,x4]:
#    plt.plot(ejey_aux1,x*np.ones(10),'--',color = 'grey' )
im = plt.imshow(Z_num, extent = limits, cmap=plt.cm.hot, aspect='auto', interpolation = 'none',origin = 'lower') 
#plt.plot(maxis,0.18*np.ones(len(maxis)),'o',color = 'green' )
cbar = plt.colorbar(im, fraction=0.046, pad=0.04, orientation = 'vertical')
cbar.ax.tick_params(labelsize = tamnum, width=0.5, length = 2,pad = 1.5)
cbar.set_label(r'$\theta$ ' +  '(' + u"\N{DEGREE SIGN}" + ')',fontsize=tamlegend,labelpad = 2)
plt.tight_layout()
os.chdir(path_save)
plt.savefig( 'theta_N%i' %(Nmax) + '.png', format='png',bbox_inches='tight',pad_inches = 0,dpi = dpi)   

#%%         
        












    