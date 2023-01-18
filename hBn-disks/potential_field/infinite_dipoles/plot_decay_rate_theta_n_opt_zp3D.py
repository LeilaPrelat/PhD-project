
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
from scipy.interpolate import interp1d

#import math

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/potential_field/infinite_dipoles','')
path_save = path_basic + '/' + 'decay_rate_theta_n'
if not os.path.exists(path_save):
    print('Creating folder to save graphs')
    os.mkdir(path_save)

err = 'decay_rate_theta_n.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from decay_rate_theta_n import decay_rate_theta_inf_dipoles_ana
except ModuleNotFoundError:
    print(err)

try:
    sys.path.insert(1, path_constants)
    from hBn_PP import hBn_lambda_p, hBn_Rp, epsilon_x, epsilon_z, hBn_lambda_p_Gself_image, hBn_Rp_Gself_image
except ModuleNotFoundError:
    print('hBn_PP.py no se encuentra en ' + path_constants)


try:
    sys.path.insert(1, path_constants)
    from constants import constantes
except ModuleNotFoundError:
    print('constants.py no se encuentra en ' + path_constants)

pi,hb,c,alfac,mu1,mu2 = constantes()
aux = c*hb
aux2 = 1e12/c

#%%

epsi1, epsi3 = 1,1

print('Definir parametros del problema')


b = - 0.01


int_v = 20
int_v = 10

Nmax = 1

#
if Nmax == 0:
    theta_degree = 0
else:
    theta_degree = 0
theta = theta_degree*np.pi/180

d_nano = 10 
tabla = np.loadtxt('zp_optimum_for_decay_rate_d%inm_pall.txt' %(d_nano), delimiter='\t', skiprows=1)
tabla = np.transpose(tabla)
[listx,listy] = tabla
f = interp1d(listx, listy)  ## energy and z_p 

N = 1e3
listx_2 = np.linspace(np.min(listx),np.max(listx),N)
if Nmax == 0:
    listx_3 = listx_2
elif Nmax == 1:
    listx_3 = np.linspace(np.min(listx),0.205,N)
elif Nmax == 2:
    listx_3 = np.linspace(150,500,N)
    
    
listy_2 = np.array(f(listx_2)) ## zp 

list_a_nm = np.linspace(0.1,1000,len(listx))
list_a_nm = np.linspace(200,600,N)
#    listx_2 = zp/(2*np.pi/listx_3)


labelx = r'Plasmon energy, $\hbar\omega$ (eV)'  
labely = r'Lattice period, $a$ (nm)'

title2 = r'v = c/%i, b = %i nm' %(int_v,b*1e3) 
title4 = r'$z_p$ = $z^{opt}_p$, n = %i, hBN disks, D = 120 nm' %(Nmax)


labelp = r'_zp_opt_theta%i_Nmax%i_v%i' %(105,Nmax,int_v)
label1 = ''  + labelp

    
energy0_pol = 0.18
omega0 = energy0_pol/hb 
#R = 10 # 10 nm en unidades de micrometros
kappa_factor_omega0 = 0.1
kappa_r_factor= 0.5   

title1 = r'$\kappa$ = %.2f$\omega_o$, $\kappa_r$ = %.2f$\kappa$, $\hbar\omega_o$ = %.2f eV' %(kappa_factor_omega0, kappa_r_factor, energy0_pol)    
title = title1 + '\n' + title2 + '\n'   + title4

def theta_f(int_v,a_nm,Nmax,E):
    
#    def omega_n_THz(int_v,a_nm,Nmax):  ## omega puede valer esto o menos para que se generen plasmones
#        
#        a = a_nm*1e-3
#        
#        aux = alfac*int_v*hbmu/(hb)
#        rta = aux + 0.5*np.sqrt((2*aux)**2 + 16*alfac*c*hbmu*Nmax*np.pi/(a*hb))
#        
#        return rta*1.5
    
#    omega_n = omega_n_THz(int_v,a_nm,Nmax)
#    omegac = omega_n/c
#    print('omega/c:', omegac)
    
 #   omegac = omegac - 0.5 ## ver si funciona con una freq menor 

    omegac = E/(hb*c)
    d_micros = d_nano*1e-3
    Rp = hBn_Rp_Gself_image(E,epsi1,epsi3)
    lambda_p_v = hBn_lambda_p_Gself_image(E,epsi1,epsi3)*d_micros
    kp = 2*np.pi/lambda_p_v
    alfa_p = kp/omegac
    

    kp = alfa_p*omegac    
    
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
    
    return theta0*180/np.pi

def function_real_ana(zp_nano,energy0_eV,a_nm):
                
    omegac0 = energy0_eV/(c*hb)  
    zp = zp_nano*1e-3
    a = a_nm*1e-3
    theta = theta_f(int_v,a_nm,Nmax,energy0_eV)
    
    theta_degree = 105
    theta = theta_degree*np.pi/180
    
    rta = decay_rate_theta_inf_dipoles_ana(omegac0,epsi1,epsi3,d_nano,int_v,zp,a,b,Nmax,omega0,kappa_factor_omega0,kappa_r_factor,theta)
#    print(rta)

    print(rta)

    rta2 = np.log10(rta)
    
    return rta2

#    return rta2


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
lw = 1.5
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
corte = int(N/2)
if Nmax == 0:
    list_z = []
    for a_nm in list_a_nm :
        list_z_1 = []
        for ind in range(corte,len(listx_3)):## decay rate negativo
            
            energy0_eV = listx_3[ind]
            zp_nano = listy_2[ind]
            

            
            list_z_1.append(function_real_ana(zp_nano,energy0_eV,a_nm))
                
        list_z.append(list_z_1)

else:
#    listx2 =  np.linspace(100,350,len(listy))
    list_z = []
    for a_nm in list_a_nm :
        list_z_1 = []
        for ind in range(len(listx_3)):
            energy0_eV = listx_3[ind]
            zp_nano = listy_2[ind]
    
            list_z_1.append(function_real_ana(zp_nano,energy0_eV,a_nm))
        list_z.append(list_z_1)


#%%

limits = [np.min(listx_3) , np.max(listx_3), np.min(list_a_nm) , np.max(list_a_nm)]
graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
im = plt.imshow(list_z, extent = limits, cmap=plt.cm.hot, aspect='auto', interpolation = 'none',origin = 'lower')  
cbar = plt.colorbar(im, extend='both', fraction=0.046, pad=0.04)
cbar.set_label('log($\Gamma_{SPE}$)',fontsize=tamlegend,labelpad = 1)
cbar.ax.tick_params(labelsize = tamnum)
plt.tight_layout()
os.chdir(path_save)
plt.savefig( 'SPE_decay_rate' + labelp + '.png', format='png')   

#%%
    


