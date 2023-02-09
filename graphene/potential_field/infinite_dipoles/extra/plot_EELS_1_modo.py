
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

plot_vs_E = 0

plot_vs_v = 0

plot_vs_a = 0

plot_vs_zp = 1

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/potential_field/infinite_dipoles','')
path_save = path_basic + '/' + 'EELS_1_modo'
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

e_charge = 4.8032*1e-10

#%%

print('Definir parametros del problema')

def EELS(omegac,epsi1,epsi2,hbmu,hbgama,y,zp,int_v,px,py,pz,a,b,n):

    kx = omegac*int_v - 2*np.pi*n/a 
    
    E = omegac*aux

    cond = 4*pi*alfac*sigma_DL(E,hbmu,hbgama)
    # Rp = 2*epsi1/(epsi1 + epsi2)
    # alfa_p = 1j*(epsi1 + epsi2)/(cond*cte1)
    # kp = alfa_p*k1     

    alfa_p = 1j*(epsi1 + epsi2)/cond
    kp = alfa_p*omegac
    
    Rp = 2*epsi1/(epsi1 + epsi2)   
        
    cte = Rp/(2*(np.pi**2)*E)
    
    term_den = np.sqrt(kp**2 - kx**2)
    
    term = 1j*kp*kx*(px*kx/term_den + py + 1j*pz*kp/term_den)*np.exp(1j*term_den*y)*np.exp(-kp*(2*zp-b))
    
    return term*cte

#%%

#v = c/int_v
#omega = 0.7*1e12
cota = 75 #nanometros
epsi1,epsi2 = 1,1
hbmu,hbgama = 0.3,0.0001

b = -0.01
#Ndipoles = 50
#list_xD = np.linspace(-0.01,0.01,Ndipoles)

px = 1
py = 0
pz = 0

y = 0

title3 = r'y = %i nm, px=%i, py=%i, pz=%i, b = %i nm' %(y*1e3,px,py,pz,b*1e3)

labelp = r'_px%i_py%i_pz%i' %(px,py,pz)

N = 500
if plot_vs_a == 1:
    
    E0 = 36 # meV 
    int_v0 = 10
    zp0 = 0.05*1e3

    labelx = 'a [nm]'
    
    title2 = r'E = %.2f meV, v = c/%i, $z_p$ = %i nm' %(E0,int_v0,zp0)
    
    title = title2 + '\n' +  title3  
    
    label1 = 'vs_a_E0%i' %(E0) + labelp
    listx = np.linspace(50,500,N)        ## en nanometros 
    

elif plot_vs_v == 1: 
    
    E0 = 36 # meV 
    a0 = 0.05*1e3   
    zp0 = 0.05*1e3
    
    
    labelx = 'v/c'  
    
    title2 = r'E = %.2f meV, a = %i nm, $z_p$ = %i nm' %(E0,a0,zp0)
    
    title = title2 + '\n' +  title3 
    
    
    label1 = 'vs_v_E0%i' %(E0) + labelp

    listx = np.linspace(0.01,0.99,N)
    
elif plot_vs_zp == 1: 
    
    E0 = 36 # meV 
    int_v0 = 10
    a0 = 0.05*1e3

    labelx = 'zp [nm]'
    
    title2 = r'E = %.2f meV, v = c/%i, a = %i nm' %(E0,int_v0,a0)
    
    title =  title2  + '\n' + title3  
    
    label1 = 'vs_zp_E0%i' %(E0) 
    listx = np.linspace(50,800,N)        ## en nanometros 
    
elif plot_vs_E == 1:
    
    a0 = 0.05*1e3    
    int_v0 = 10
    zp0 = 0.05*1e3

    
    labelx = '$\hbar\omega$ [meV]'   
    
    title2 = r'a = %i nm, v = c/%i, $z_p$ = %i nm' %(a0,int_v0,zp0)
    
    title = title2 + '\n' +  title3  
    
    label1 = 'vs_E' + labelp
    
    
    listx = np.linspace(11,60,N)

#%%

def function_re(a_nano,zp_nano,int_v,energy0,n):
    omegac0 = energy0*1e-3/aux 
    a = a_nano*1e-3 
    zp = zp_nano*1e-3

    rta = EELS(omegac0,epsi1,epsi2,hbmu,hbgama,y,zp,int_v,px,py,pz,a,b,n)
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

labely1 = r'EELS$_{ind}$ $/e$'

if plot_vs_a == 1:
    list_n = [1,2,3]
else:
    list_n = [0,1,2,3]
listy_re_ana_tot = []    



if plot_vs_a == 1: 
    
    for n in list_n:

        listy_re_ana = []    
        listy_im_ana = [] 

        for value in listx: 
    
            y_re_ana = function_re(value,zp0,int_v0,E0,n)    
            listy_re_ana.append(y_re_ana)
            
        listy_re_ana_tot.append(listy_re_ana)



elif plot_vs_v == 1: 
    
    for n in list_n:

        listy_re_ana = []    
        listy_im_ana = [] 

        for value in listx: 
    
            y_re_ana = function_re(a0,zp0,value,E0,n)    
            listy_re_ana.append(y_re_ana)
        
        
        listy_re_ana_tot.append(listy_re_ana)


elif plot_vs_zp == 1: 
    
    for n in list_n:

        listy_re_ana = []    
        listy_im_ana = [] 

        for value in listx: 
    
            y_re_ana = function_re(a0,value,int_v0,E0,n)    
            listy_re_ana.append(y_re_ana)
        
        
        listy_re_ana_tot.append(listy_re_ana)


elif plot_vs_E == 1: 

    for n in list_n:
        
        listy_re_ana = []    
        listy_im_ana = [] 
    
 
        for value in listx: 
    
            y_re_ana = function_re(a0,zp0,int_v0,value,n)     
            listy_re_ana.append(y_re_ana)
    

        listy_re_ana_tot.append(listy_re_ana)

        
#%%       
graph(title,labelx,labely1,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)

k = 0
for modo in list_n:
    plt.plot(listx,listy_re_ana_tot[k],'.-',ms = ms, label = 'n = %i' %(modo))
    k = k + 1
    
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
os.chdir(path_save)
plt.savefig( 'EELS_ind' + label1 + '.png', format='png')   
   
#%%
