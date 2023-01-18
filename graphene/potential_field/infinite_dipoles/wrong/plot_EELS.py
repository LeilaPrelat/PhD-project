
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

plot_vs_y = 1
plot_vs_E = 0

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/potential_field/infinite_dipoles','')
path_save = path_basic + '/' + 'EELS_1'
if not os.path.exists(path_save):
    print('Creating folder to save graphs')
    os.mkdir(path_save)

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

#%%

print('Definir parametros del problema')

def EELS(omegac,epsi1,epsi2,hbmu,hbgama,y,zp,int_v,px,py,pz,a):
    E = omegac*aux
    omega = omegac*c
    cond = 4*pi*alfac*sigma_DL(E,hbmu,hbgama)
    # Rp = 2*epsi1/(epsi1 + epsi2)
    # alfa_p = 1j*(epsi1 + epsi2)/(cond*cte1)
    # kp = alfa_p*k1     

    alfa_p = 1j*(epsi1 + epsi2)/cond
    kp = alfa_p*omegac
    
    Rp = 2*epsi1/(epsi1 + epsi2)   
    
    term_den = np.sqrt(kp**2 - (omegac*int_v)**2)
    
    cte = 1j*Rp*omegac*int_v/(2*(np.pi**2)*a*omega)
    
    term = kp*(px*omegac*int_v/term_den + py + 1j*pz*kp/term_den)*np.exp(1j*term_den*y)*np.exp(-kp*(2*zp-b))
    
    return np.real(term*cte)

#v = c/int_v
#omega = 0.7*1e12
cota = 75 #nanometros
epsi1,epsi2 = 1,1
hbmu,hbgama = 0.3,0.0001
zp = 0.05
b = -0.01

x0,z0 = 0, zp 

#Ndipoles = 50
#list_xD = np.linspace(-0.01,0.01,Ndipoles)
a = 0.05

px = 1
py = 0
pz = 0

int_v = 10

title2 = r'$\hbar\mu$ = %.2feV, $\hbar\gamma$ = %.4feV, v = c/%i, b = %inm' %(hbmu,hbgama,int_v,b*1e3) 
title3 = r'$z_p$=%inm, px=%i, py=%i, pz=%i, a = %inm' %(zp*1e3,px,py,pz,a*1e3)

labelp = r'_px%i_py%i_pz%i_a%inm' %(px,py,pz,a*1e3)

N = 500
if plot_vs_y == 1: 
    E0 = 60 # meV 
    labelx = 'y [nm]'  
    title = title2 + '\n' +  title3  + ', E = %.2f meV' %(E0)
    label1 = 'vs_y_E0%i' %(E0) + labelp
    listx = np.linspace(50,6000,N)
    
    listx = np.linspace(-50,80000,N)
else:
    y0 = 400
#    y0 = 0
    # z0 = 0.06*1e3
    labelx = 'E [meV]'   
    title = title2 + '\n' + title3  + ', y = %inm' %(y0)
    label1 = 'vs_E_y0%.2fnm' %(y0) + labelp
    listx = np.linspace(20,60,N)

#%%

def function(y_nano,energy0):
    omegac0 = energy0*1e-3/aux 
    y = y_nano*1e-3 

    rta = EELS(omegac0,epsi1,epsi2,hbmu,hbgama,y,zp,int_v,px,py,pz,a)
    # print(rta)
    return rta


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

if plot_vs_y == 1: 
    
    listy_re_ana = []    

    for value in listx: 

        y_re_ana = function(value,E0)    
        listy_re_ana.append(y_re_ana)

    graph(title,labelx,'EELS/(e/$\hbar$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    plt.plot(listx,listy_re_ana,'.-',ms = ms,color = 'purple')
#    plt.plot(listx,listy_re_anaG1_v2,'.',ms = ms,color = 'pink',label = 'analytical v2')
#    plt.plot(listx,listy_re_num,'.',ms = 4,color = 'lightseagreen',label ='full numerical')
#    plt.plot(listx,listy_re_pole_approx,'.-',ms = 2,color = 'darkred',label =  'PP numerical')
    plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'EELS' + label1 + '.png', format='png')   

if plot_vs_E == 1: 
    
    listy_re_ana = []    

    for value in listx: 

        y_re_ana = function(y0,value)     
        listy_re_ana.append(y_re_ana)

#        
    graph(title,labelx,'EELS/(e/$\hbar$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    plt.plot(listx,listy_re_ana,'.-',ms = ms,color = 'purple')
#    plt.plot(listx,listy_re_anaG1_v2,'.',ms = ms,color = 'pink',label = 'analytical v2')
#    plt.plot(listx,listy_re_num,'.',ms = 4,color = 'lightseagreen',label = 'full numerical')
#    plt.plot(listx,listy_re_pole_approx,'.-',ms = 2,color = 'darkred',label = 'PP numerical')
    plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'EELS' + label1 + '.png', format='png')   
   
#%%
