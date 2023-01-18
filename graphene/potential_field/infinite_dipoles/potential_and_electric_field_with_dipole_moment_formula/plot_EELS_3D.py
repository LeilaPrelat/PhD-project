
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
from scipy import special

#%%
paper = 1

plot_vs_b = 1
plot_vs_v = 0
plot_vs_y = 0

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/potential_field/infinite_dipoles','')
path_save = path_basic + '/' + 'EELS_3D'
if not os.path.exists(path_save):
    print('Creating folder to save graphs')
    os.mkdir(path_save)


try:
    sys.path.insert(1, path_constants)
    from dipole_moment import dipole_moment_anav2
except ModuleNotFoundError:
    print('dipole_moment.py no se encuentra en ' + path_constants)
    
try:
    sys.path.insert(1, path_constants)
    from graphene_sigma import sigma_DL
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

e_charge = 4.8032*1e-10

#%%

print('Definir parametros del problema')



#v = c/int_v
#omega = 0.7*1e12
cota = 75 #nanometros
epsi1,epsi2 = 1,1
hbmu,hbgama = 0.3,0.0001
zp = 0.05

energy0_pol = 43
omega0 = energy0_pol*1e-3/hb 
#R = 10 # 10 nm en unidades de micrometros
kappa_factor_omega0 = 0.1
kappa_r_factor= 0.5

x0,z0 = 0, zp 

#Ndipoles = 50
#list_xD = np.linspace(-0.01,0.01,Ndipoles)
a = 0.05




N = 150
title1 = r'$\kappa$ = %.2f$\omega_0$, $\kappa_r$ = %.2f$\kappa$, $E_0$=%i meV' %(kappa_factor_omega0, kappa_r_factor, energy0_pol)    
#title2 = r'$\hbar\mu$ = %.2feV, $\hbar\gamma$ = %.4feV' %(hbmu,hbgama) 
title2 = r'$z_p$=%inm, a = %inm' %(zp*1e3,a*1e3)

labelp = r'_a%inm' %(a*1e3)


def EELS(omegac,epsi1,epsi2,hbmu,hbgama,y,zp,int_v,a,b):
    E = omegac*aux
    omega = omegac*c
    cond = 4*pi*alfac*sigma_DL(E,hbmu,hbgama)
    # Rp = 2*epsi1/(epsi1 + epsi2)
    # alfa_p = 1j*(epsi1 + epsi2)/(cond*cte1)
    # kp = alfa_p*k1     

    alfa_p = 1j*(epsi1 + epsi2)/cond
    kp = alfa_p*omegac
    
    Rp = 2*epsi1/(epsi1 + epsi2)   
    
    px,py,pz  = dipole_moment_anav2(omegac,epsi1,epsi2,hbmu,hbgama,int_v,b,zp,omega0,kappa_factor_omega0,kappa_r_factor)    
    
    cte = Rp/(4*(np.pi**3)*a*E)
    
    zero_cte = omegac*int_v
    
    term_den = np.sqrt(kp**2 - zero_cte**2)
    
    term = 1j*kp*zero_cte*(px*zero_cte/term_den + py + 1j*pz*kp/term_den)*np.exp(1j*term_den*y)*np.exp(-kp*(2*zp-b))
    
    return np.real(term*cte)



def EELS_dir(omegac,epsi1,epsi2,hbmu,hbgama,y,zp,int_v,a,b):
    E = omegac*aux
    omega = omegac*c

    # Rp = 2*epsi1/(epsi1 + epsi2)
    # alfa_p = 1j*(epsi1 + epsi2)/(cond*cte1)
    # kp = alfa_p*k1     
    px,py,pz  = dipole_moment_anav2(omegac,epsi1,epsi2,hbmu,hbgama,int_v,b,zp,omega0,kappa_factor_omega0,kappa_r_factor)
    
    cte = 1/(2*(np.pi**3)*a*E)
    arg = omegac*int_v*np.abs(b)
    K0 =  special.kv(0,arg)    
    K1 =  special.kv(1,arg)    
    zero_cte = 3*omegac*int_v
    
    term = 1j*(zero_cte**2)*( -1j*px*K0 - pz*np.sign(b)*K1)
    
#    print(np.real(K0), np.imag(K1))
#    
    return np.real(term*cte)

if plot_vs_v == 1:
    b = -0.01    
    y0 = 0
    def EELS_3D_ind(int_v_value,energy0):
        omegac0 = energy0*1e-3/aux     
        return EELS(omegac0,epsi1,epsi2,hbmu,hbgama,y0,zp,1/int_v_value,a,b)

    def EELS_3D_dir(int_v_value,energy0):
        omegac0 = energy0*1e-3/aux     
        return EELS_dir(omegac0,epsi1,epsi2,hbmu,hbgama,y0,zp,1/int_v_value,a,b)

    listx = np.linspace(0.01,0.2,N)
    labelx = 'v/c'  
    title3 = 'b = %inm , y = %inm' %(b*1e3,y0*1e3)
    title = title2 + '\n' + title3
    labelp = labelp + '_vs_v'
    listy = np.linspace(15,75,N)    
    
elif plot_vs_b == 1:    
    int_v = 10
    y0 =0
    def EELS_3D_ind(b_value_nm,energy0):
        omegac0 = energy0*1e-3/aux 
        b_value = b_value_nm*1e-3
        return EELS(omegac0,epsi1,epsi2,hbmu,hbgama,y0,zp,int_v,a,b_value)


    def EELS_3D_dir(b_value_nm,energy0):
        omegac0 = energy0*1e-3/aux 
        b_value = b_value_nm*1e-3
        return EELS_dir(omegac0,epsi1,epsi2,hbmu,hbgama,y0,zp,int_v,a,b_value)
   
    listx1 = np.linspace(-0.05*1e-3,0.00001,N)
    listx = np.array(listx1)*1e3
    labelx = 'b [nm]'      
    title3 = 'v = c/%i, y = %inm' %(int_v,y0*1e3)
    title = title2 + '\n' + title3 
    labelp = labelp + '_vs_b'
    
    listy = np.linspace(15,55,N)

elif plot_vs_y == 1:    
    int_v = 10
    b = -0.01     
    def EELS_3D_ind(y_value,energy0):
        omegac0 = energy0*1e-3/aux     
        return EELS(omegac0,epsi1,epsi2,hbmu,hbgama,y_value,zp,int_v,a,b)

    def EELS_3D_dir(y_value,energy0):
        omegac0 = energy0*1e-3/aux     
        return EELS_dir(omegac0,epsi1,epsi2,hbmu,hbgama,y_value,zp,int_v,a,b)
    
    listx1 = np.linspace(0,500,N)
    listx = np.array(listx1)
    labelx = 'y [nm]'      
    title = title2 + '\n' +  title3  + ', v = c/%i' %(int_v)
    labelp = labelp + '_vs_y'


labely = r'$\hbar\omega$ [meV]'   


#%%
    
tamfig = (4.5,3.5)
tamlegend = 12
tamletra = 12
tamtitle = 11
tamnum = 10
labelpady = -0.5
labelpadx = -0.5
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
    if paper == 0:
        plt.title(title,fontsize=int(tamtitle*0.9))

    return   

#%%
    
labelz1 = r'EELS direct $/(Le)$'
labelz2 = r'EELS induce $/(Le)$'

X, Y = np.meshgrid(listx, listy, sparse=True)
f_dir = np.vectorize(EELS_3D_dir)
Z_dir = f_dir(X, Y)

f_ind = np.vectorize(EELS_3D_ind)
Z_ind = f_ind(X, Y)

#%%    

limits = [np.min(listx) , np.max(listx), np.min(listy) , np.max(listy)]
graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)

# im = plt.imshow(Z, extent = limits,  cmap='RdBu', interpolation='bilinear') 
if plot_vs_y == 1:
    im = plt.imshow(Z_dir, extent = limits, cmap=plt.cm.hot, aspect=5) 
elif plot_vs_b:
    im = plt.imshow(Z_dir, extent = limits, cmap=plt.cm.hot, aspect=1/1000) 
else:
    im = plt.imshow(Z_dir, extent = limits, cmap=plt.cm.hot, aspect=1/300) 
cbar = plt.colorbar(im, extend='both', fraction=0.046, pad=0.04)
cbar.ax.tick_params(labelsize = tamnum)
if paper == 0:
    cbar.set_label(labelz1,fontsize=tamlegend,labelpad = 1)
plt.tight_layout()
os.chdir(path_save)
plt.savefig( 'EELS_3D_dir' + labelp + '.png', format='png')   

#%%


limits = [np.min(listx) , np.max(listx), np.min(listy) , np.max(listy)]
graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)

# im = plt.imshow(Z, extent = limits,  cmap='RdBu', interpolation='bilinear') 
if plot_vs_y == 1:
    im = plt.imshow(Z_ind, extent = limits, cmap=plt.cm.hot, aspect=5) 
elif plot_vs_b:
    im = plt.imshow(Z_ind, extent = limits, cmap=plt.cm.hot, aspect=1/1000) 
else:
    im = plt.imshow(Z_ind, extent = limits, cmap=plt.cm.hot, aspect=1/300) 
cbar = plt.colorbar(im, extend='both', fraction=0.046, pad=0.04)
cbar.ax.tick_params(labelsize = tamnum)
if paper == 0:
    cbar.set_label(labelz2,fontsize=tamlegend,labelpad = 1)
plt.tight_layout()
os.chdir(path_save)
plt.savefig( 'EELS_3D_ind' + labelp + '.png', format='png')   