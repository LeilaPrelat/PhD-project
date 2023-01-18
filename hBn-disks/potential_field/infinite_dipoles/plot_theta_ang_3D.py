
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
    from hBn_PP import hBn_lambda_p, hBn_Rp, epsilon_x, epsilon_z, hBn_lambda_p_Gself_image, hBn_Rp_Gself_image
except ModuleNotFoundError:
    print('graphene_sigma.py no se encuentra en ' + path_constants)

    
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

d_nano = 10
#Ndipoles = 50
#list_xD = np.linspace(-0.01,0.01,Ndipoles)


int_v = 10
Nmax = 1
#def omega_n_THz(int_v,a_nm,Nmax):  ## omega puede valer esto o menos para que se generen plasmones
#    
#    a = a_nm*1e-3
#    
#    aux = alfac*int_v*hbmu/(hb)
#    rta = aux + 0.5*np.sqrt((2*aux)**2 + 16*alfac*c*hbmu*Nmax*np.pi/(a*hb))
#    
#    return rta*1.5

title = r'v = c/%i, n = %i, d = %i nm' %(int_v,Nmax,d_nano) 
labelx = r'$\hbar\omega$ [eV]'  
labely = r'$a$ [nm]'
def theta(E,a_nm):
    

#    
#    omega_n = omega_n_THz(int_v,a_nm,Nmax)
#    omegac = omega_n/c
#    print('omega/c:', omegac)
    
 #   omegac = omegac - 0.5 ## ver si funciona con una freq menor 


    omegac = E/(hb*c)
    d_micros = d_nano*1e-3
    Rp = hBn_Rp_Gself_image(E,epsi1,epsi3)
    lambda_p_v = hBn_lambda_p_Gself_image(E,epsi1,epsi3)*d_micros
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
    a = a_nm*1e-3
    
    theta0 = np.arccos((omegac*int_v + 2*np.pi*Nmax/a)/np.real(kp))
#    print(theta0)
    return theta0*180/np.pi


    
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
N2 = 400


x1 = 0.09260651629072682 
x2 = 0.10112781954887218
x3 = 0.17030075187969923
x4 = 0.19937343358395992

listx = np.linspace(0.09,0.195,N2)
if Nmax == 0:
    listx = np.linspace(x1+1e-3,x4-1e-3,N2)
else:
    listx = np.linspace(x1+1e-3,0.85,N2)


listy = np.linspace(200,6000,N2)

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
cbar = plt.colorbar(im, extend='both', fraction=0.046, pad=0.04)
cbar.ax.tick_params(labelsize = tamnum)
cbar.set_label(r'$\theta$ [degrees]',fontsize=tamlegend,labelpad = 1)
plt.tight_layout()
os.chdir(path_save)
plt.savefig( 'theta_N%i' %(Nmax) + '.png', format='png')   

#%%         
        












    