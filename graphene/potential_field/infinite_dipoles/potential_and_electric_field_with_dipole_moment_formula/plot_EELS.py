
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
import seaborn as sns
sns.set()

#%%

plot_vs_y = 0
plot_vs_E = 1

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/potential_field/infinite_dipoles/potential_and_electric_field_with_dipole_moment_formula','')
path_save = path_basic + '/' + 'EELS'
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
b = -0.01

x0,z0 = 0, zp 

energy0_pol = 43
omega0 = energy0_pol*1e-3/hb 
#R = 10 # 10 nm en unidades de micrometros
kappa_factor_omega0 = 0.1
kappa_r_factor= 0.5

#Ndipoles = 50
#list_xD = np.linspace(-0.01,0.01,Ndipoles)
a = 0.05



int_v = 10

title1 = r'$\kappa$ = %.2f$\omega_0$, $\kappa_r$ = %.2f$\kappa$, $E_0$=%i meV' %(kappa_factor_omega0, kappa_r_factor, energy0_pol)    
title2 = r'v = c/%i, b = %inm' %(int_v,b*1e3) 
title3 = r'$z_p$=%inm, a = %inm' %(zp*1e3,a*1e3)

labelp = r'_a%inm' %(a*1e3)

N = 500
if plot_vs_y == 1: 
    E0 = 36 # meV 
    labelx = 'y [nm]'  
    title = title1 + '\n' + title2 + ', ' +  title3  + ', E = %.2f meV' %(E0)
    label1 = 'vs_y_E0%i' %(E0) + labelp
    listx = np.linspace(50,6000,N)
    
    listx = np.linspace(-50,40000,N)
else:
    y0 = 0
#    y0 = 0
    # z0 = 0.06*1e3
    labelx = '$\hbar\omega$ [meV]'   
    title = title1 + '\n' + title2 + ', ' + title3  + ', y = %inm' %(y0)
    label1 = 'vs_E_y0%.2fnm' %(y0) + labelp
    listx = np.linspace(11,60,N)

#%%
    
    
def EELS(omegac,epsi1,epsi2,hbmu,hbgama,y,zp,int_v,a,b):
    E = omegac*aux
    omega = omegac*c
    cond = 4*pi*alfac*sigma_DL(E,hbmu,hbgama)
    # Rp = 2*epsi1/(epsi1 + epsi2)
    # alfa_p = 1j*(epsi1 + epsi2)/(cond*cte1)
    # kp = alfa_p*k1     


    px,py,pz  = dipole_moment_anav2(omegac,epsi1,epsi2,hbmu,hbgama,int_v,b,zp,omega0,kappa_factor_omega0,kappa_r_factor)
    
    alfa_p = 1j*(epsi1 + epsi2)/cond
    kp = alfa_p*omegac
    
    Rp = 2*epsi1/(epsi1 + epsi2)   
    
    
    
    cte = Rp/(4*(np.pi**3)*a*E)
    
    zero_cte = omegac*int_v
    
    term_den = np.sqrt(kp**2 - zero_cte**2)
    
    term = -1j*kp*zero_cte*(-px*zero_cte/term_den + py + 1j*pz*kp/term_den)*np.exp(1j*term_den*y)*np.exp(-kp*(2*zp-b))
    
    return term*cte



def EELS_dir(omegac,epsi1,epsi2,hbmu,hbgama,y,zp,int_v,a,b):
    E = omegac*aux
    omega = omegac*c
    cond = 4*pi*alfac*sigma_DL(E,hbmu,hbgama)
    # Rp = 2*epsi1/(epsi1 + epsi2)
    # alfa_p = 1j*(epsi1 + epsi2)/(cond*cte1)
    # kp = alfa_p*k1     

    alfa_p = 1j*(epsi1 + epsi2)/cond
    kp = alfa_p*omegac
    
    Rp = 2*epsi1/(epsi1 + epsi2)   
    
    arg = omegac*int_v*np.abs(b)
    K0 =  special.kv(0,arg)    
    K1 =  special.kv(1,arg)    
    
    cte = Rp/(2*(np.pi**3)*a*E)
    
    zero_cte = omegac*int_v


    px,py,pz  = dipole_moment_anav2(omegac,epsi1,epsi2,hbmu,hbgama,int_v,b,zp,omega0,kappa_factor_omega0,kappa_r_factor)
    
    term_den = np.sqrt(kp**2 - zero_cte**2)
    
    term = 1j*(zero_cte**2)*( 1j*px*K0 - pz*np.sign(b)*K1)
    
    return term*cte

def function_re(y_nano,energy0):
    omegac0 = energy0*1e-3/aux 
    y = y_nano*1e-3 

    rta = EELS(omegac0,epsi1,epsi2,hbmu,hbgama,y,zp,int_v,a,b)
    # print(rta)
    return np.real(rta)

def function_im(y_nano,energy0):
    omegac0 = energy0*1e-3/aux 
    y = y_nano*1e-3 

    rta = EELS(omegac0,epsi1,epsi2,hbmu,hbgama,y,zp,int_v,a,b)
    # print(rta)
    return np.imag(rta)


def function_re_dir(y_nano,energy0):
    omegac0 = energy0*1e-3/aux 
    y = y_nano*1e-3 

    rta = EELS_dir(omegac0,epsi1,epsi2,hbmu,hbgama,y,zp,int_v,a,b)
    # print(rta)
    return np.real(rta)

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
labely1 = r'EELS$_{ind}$ $/(Le)$'
labely2 = r'EELS$_{dir}$ $/(Le)$'
if plot_vs_y == 1: 
    
    listy_re_ana = []    
    listy_im_ana = [] 
    for value in listx: 

        y_re_ana = function_re(value,E0)    
        listy_re_ana.append(y_re_ana)
        
        y_im_ana = function_im(value,E0)    
        listy_im_ana.append(y_im_ana)
        

if plot_vs_E == 1: 
    
    listy_re_ana = []    
    listy_im_ana = [] 

    listy_re_ana_dir = []   
    for value in listx: 

        y_re_ana = function_re(y0,value)     
        listy_re_ana.append(y_re_ana)

        y_im_ana = function_im(y0,value)    
        listy_im_ana.append(y_im_ana)
        
        y_re_ana_dir = function_re_dir(y0,value)     
        listy_re_ana_dir.append(y_re_ana_dir)

        
#%%       
graph(title,labelx,labely1,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(listx,listy_re_ana,'.-',ms = ms,color = 'purple')
#    plt.plot(listx,listy_re_anaG1_v2,'.',ms = ms,color = 'pink',label = 'analytical v2')
#    plt.plot(listx,listy_re_num,'.',ms = 4,color = 'lightseagreen',label = 'full numerical')
#    plt.plot(listx,listy_re_pole_approx,'.-',ms = 2,color = 'darkred',label = 'PP numerical')
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
os.chdir(path_save)
plt.savefig( 'EELS_ind' + label1 + '.png', format='png')   
   

graph(title,labelx,labely2,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(listx,listy_re_ana_dir,'.-',ms = ms,color = 'darkred')
#    plt.plot(listx,listy_re_anaG1_v2,'.',ms = ms,color = 'pink',label = 'analytical v2')
#    plt.plot(listx,listy_re_num,'.',ms = 4,color = 'lightseagreen',label = 'full numerical')
#    plt.plot(listx,listy_re_pole_approx,'.-',ms = 2,color = 'darkred',label = 'PP numerical')
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
os.chdir(path_save)
plt.savefig( 'EELS_dir' + label1 + '.png', format='png')   

#%%

#graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#plt.plot(listx,listy_im_ana,'.-',ms = ms,color = 'purple')
##    plt.plot(listx,listy_re_anaG1_v2,'.',ms = ms,color = 'pink',label = 'analytical v2')
##    plt.plot(listx,listy_re_num,'.',ms = 4,color = 'lightseagreen',label = 'full numerical')
##    plt.plot(listx,listy_re_pole_approx,'.-',ms = 2,color = 'darkred',label = 'PP numerical')
#plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
#plt.tight_layout()
#os.chdir(path_save)
#plt.savefig( 'EELS_ind_im' + label1 + '.png', format='png')   
#   
