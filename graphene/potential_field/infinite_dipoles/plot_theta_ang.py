
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

plot_vs_v = 0
plot_vs_a = 1

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


N = 1500
if plot_vs_v == 1:
    list_n = [1,2,3]
else:
    list_n = [1,5,10,15]


def omega_n_THz(int_v,a_nm,Nmax):  ## omega puede valer esto o menos para que se generen plasmones
    
    a = a_nm*1e-3
    
    aux = alfac*int_v*hbmu/(hb)
    rta = aux + 0.5*np.sqrt((2*aux)**2 + 16*alfac*c*hbmu*Nmax*np.pi/(a*hb))
    
    return rta*1.5



def theta(int_v,a_nm,Nmax):
    

    
    omega_n = omega_n_THz(int_v,a_nm,Nmax)
    omegac = omega_n/c
    print('omega/c:', omegac)
    
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
    
    theta0 = np.arccos((omegac*int_v + 2*np.pi*Nmax/a)/np.real(kp))
    
    return theta0*180/np.pi

if plot_vs_v == 1 :
    
    
#    
#    omega_n = omega_n_THz(int_v,a_nm,Nmax)*1e12
#    lambda0_nm =  2*pi*v/omega_n
    
    a_nm = 500
#    
#    v_c_min_tot = []
#    v_c_max_tot = []
#    for Nmax in list_n: 
#        v_c_min = (1 + lambda0_nm*Nmax/a_nm)**(-1)
#        v_c_max = (-1 + lambda0_nm*Nmax/a_nm)**(-1)
#        
#        v_c_min_tot.append(v_c_min)
#        v_c_max_tot.append(v_c_max)
    
    
    labelx = 'v/c' 
    label1 = 'vs_v' 
    
    title = 'a = %i nm, $\mu$ = %.2f eV' %(a_nm,hbmu)
    
#    if np.max(v_c_max_tot) <= 1:
#        listx = np.linspace(np.min(v_c_min_tot),np.max(v_c_max_tot),N)
#    else:
    listx = np.linspace(0.001,1,N)
    
elif plot_vs_a == 1:
    
    int_v0 = 10
#    E0 = 40
#    lambda0_nm =  2*np.pi/k_SP(E0,hbmu,hbgama)
#    lambda0_nm = np.real(lambda0_nm)
#    
#    a_min_tot = []
#    a_max_tot = []
#    for Nmax in list_n: 
#        a_min = lambda0_nm*Nmax/(int_v0 - 1)
#        a_max = lambda0_nm*Nmax/(int_v0 + 1)
#        
#        a_min_tot.append(a_min)
#        a_max_tot.append(a_max)
    
    labelx = 'a [nm]' 
    label1 = 'vs_a' 
    
    title = 'v = c/%i, $\mu$ = %.2f eV' %(int_v0,hbmu)
    
    listx = np.linspace(1, 5000,N)
    
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



listy_im_ana_tot = []

list_freqTHz_tot = []
if plot_vs_v == 1: 

    for n in list_n:


        listy_im_ana = []
        
        list_freq = []
        
        for value in listx: 
    
            y_im_ana = theta(value,a_nm,n)        
            listy_im_ana.append(y_im_ana)
            
            freq = omega_n_THz(value,a_nm,n)
            list_freq.append(freq)
            
            
        listy_im_ana_tot.append(listy_im_ana)
        list_freqTHz_tot.append(list_freq)

elif plot_vs_a == 1:       

    for n in list_n:


        listy_im_ana = []    
         
        list_freq = []
        for value in listx: 
            
    
            y_im_ana = theta(1/int_v0,value,n)        
            listy_im_ana.append(y_im_ana)

            freq = omega_n_THz(value,a_nm,n)
            list_freq.append(freq)
            

        listy_im_ana_tot.append(listy_im_ana)
        list_freqTHz_tot.append(list_freq)
            
#%%         
    
colors = ['darkred','steelblue','coral','yellowgreen']    

graph(title,labelx,r'$\theta_n$ [degrees]',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
k = 0
for modo in list_n:
    plt.plot(listx,listy_im_ana_tot[k],'.-',color = colors[k],ms = ms, label = 'n = %i' %(modo))
    k = k + 1

if plot_vs_a == 1:
    plt.legend(loc = 'upper right',markerscale=1.5,fontsize=tamlegend,frameon=0.2,handletextpad=0.2, handlelength=length_marker)
else:
    plt.legend(loc = 'best',markerscale=1.5,fontsize=tamlegend,frameon=0.2,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
#if plot_vs_c == 1:
#    plt.yscale('log')
os.chdir(path_save)
plt.savefig( 'phi' + label1 + '.png', format='png')   
                       
#%%         
        












    