
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
    from decay_rate_theta_n import decay_rate_theta_inf_dipoles_ana_res
except ModuleNotFoundError:
    print(err)

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

epsi1, epsi3 = 1,1

print('Definir parametros del problema')


b = - 0.01

d_nano = 1

int_v = 20
int_v = 10

Nmax = 1

#
if Nmax == 0:
    theta_degree = 0
else:
    theta_degree = 0
theta_degree = 30
theta = theta_degree*np.pi/180


tabla = np.loadtxt('zp_optimum_for_decay_rate_resonance_d%inm_v%i.txt' %(d_nano,int_v), delimiter='\t', skiprows=1)
tabla = np.transpose(tabla)
[listx,listy,listz] = tabla
f = interp1d(listx, listy)

N = 1e3
listx_2 = np.linspace(np.min(listx),np.max(listx),N)
if Nmax == 0:
    listx_3 = listx_2
elif Nmax == 1:
    listx_3 = np.linspace(0.5,4.5,N)
elif Nmax == 2:
    listx_3 =  np.linspace(0.5,4.5,N)
    
    
listy_2 = np.array(f(listx_2))

list_a_nm = np.linspace(0.1,1000,len(listx))
list_a_nm = np.array(np.linspace(200,600,N))*1e-3
#    listx_2 = zp/(2*np.pi/listx_3)

labelx = r'Plasmon energy, $\hbar\omega$ (eV)'  
labely = r'Lattice period, $a$ ($\mu$m)'

title2 = r'v = c/%i, b = %i nm' %(int_v,b*1e3) 
title4 = r'$z_p$ = $z^{opt}_p$, n = %i' %(Nmax)


labelp = r'_zp_opt_Nmax%i_v%i_theta%i' %(Nmax,int_v,theta_degree)
label1 = ''  + labelp

      
title =  title2 + '\n'   + title4

def theta_f(int_v,a,Nmax,E):
    
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
    Rp = Silver_Rp(E,epsi1,epsi3)
    lambda_p_v = Silver_lambda_p(E,epsi1,epsi3)*d_micros
    kp = 2*np.pi/lambda_p_v
    alfa_p = kp/omegac
    

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
#    a = a_nm*1e-3
    
    theta0 = np.arccos((omegac*int_v + 2*np.pi*Nmax/a)/np.real(kp))
    
    return theta0*180/np.pi

def function_real_ana(zp_nano,energy0,a):
                
    omegac0 = energy0/(c*hb)  
    zp = zp_nano*1e-3
#    a = a_nm*1e-3
#    theta = theta_f(int_v,a,Nmax,energy0)

    rta = decay_rate_theta_inf_dipoles_ana_res(omegac0,epsi1,epsi3,d_nano,int_v,zp,a,b,Nmax,theta)
    print(rta)

    rta2 = np.log10(rta)
    
    return rta2

#    return rta2


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
    plt.tick_params(labelsize = tamnum, pad = pad)
    plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in", pad = pad)

    return   

#%%
corte = int(N/2)
if Nmax == 0:
    list_z = []
    for a_nm in list_a_nm :
        list_z_1 = []
        for ind in range(corte,len(listx_3)):## decay rate negativo
            
            energy0_meV = listx_3[ind]
            zp_nano = listy_2[ind]
            

            
            list_z_1.append(function_real_ana(zp_nano,energy0_meV,a_nm))
                
        list_z.append(list_z_1)

else:
#    listx2 =  np.linspace(100,350,len(listy))
    list_z = []
    for a_nm in list_a_nm :
        list_z_1 = []
        for ind in range(len(listx_3)):
            energy0_meV = listx_3[ind]
            zp_nano = listy_2[ind]
    
            list_z_1.append(function_real_ana(zp_nano,energy0_meV,a_nm))
        list_z.append(list_z_1)


#%%

if Nmax == 0:
    limits = [np.min(listx_3[corte:-1]) , np.max(listx_3[corte:-1]), np.min(list_a_nm[corte:-1]) , np.max(list_a_nm[corte:-1])]
else:
    limits = [np.min(listx_3) , np.max(listx_3), np.min(list_a_nm) , np.max(list_a_nm)]
graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
im = plt.imshow(list_z, extent = limits, cmap=plt.cm.hot, aspect='auto', interpolation = 'none',origin = 'lower')  
cbar = plt.colorbar(im, fraction=0.046, pad=0.04, orientation = 'vertical')
cbar.set_label('log($\Gamma_{SPE}$) (arb. units)',fontsize=tamlegend,labelpad = 1)
cbar.ax.tick_params(labelsize = tamnum, width=0.5, length = 2,pad = 1.5)
plt.tight_layout()
os.chdir(path_save)
plt.savefig( 'SPE_decay_rate_res_' + labelp + '.png', format='png',dpi = dpi,bbox_inches='tight',pad_inches = 0.01)   

#%%
    


