
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
path_save = path_basic + '/' + 'omega_n'
if not os.path.exists(path_save):
    print('Creating folder to save graphs')
    os.mkdir(path_save)

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


N = 200
if plot_vs_v == 1:
    list_n = [0,1,2]
else:
    list_n = [1,2,3]
#def k_SP(E_mev,hbmu,hbgama):
#    num = 2/alfac 
#    omegac = E_mev*1e-3/aux
#    cte = E_mev*1e-3 + 1j*hbgama/(4*hbmu)
#    
#    return omegac*num*cte


def omega_n_THz(int_v,a_nm,Nmax):
    
    a = a_nm*1e-3
    
    aux = alfac*int_v*hbmu/(hb)
    rta = aux + 0.5*np.sqrt((2*aux)**2 + 16*alfac*c*hbmu*Nmax*np.pi/(a*hb))
    
    return rta*1e-12 


def omega_min_THz(int_v,a_nm,Nmax):
    a = a_nm*1e-3

    A = (epsi1 + epsi2)*hb/(4*hbmu*alfac*c)
    
    v = c/int_v
    
    aux = 1/(A*v)
    rta = aux*0.5  + 0.5*np.sqrt(aux**2 + 8*np.pi*n/(A*a) )
    return rta*1e-12   
    
    
if plot_vs_v == 1 :
    

#    lambda0_nm =  2*np.pi/k_SP(E0,hbmu,hbgama)
#    lambda0_nm = np.real(lambda0_nm)
    a_nm = 50
    
    
    labelx = 'v/c' 
    label1 = 'vs_v' 
    
    title = 'a = %i nm, $\mu$ = %.2f eV' %(a_nm,hbmu)
    

    listx = np.linspace(0.01,1,N)
    

    
elif plot_vs_a == 1:
    
    int_v0 = 10

#    lambda0_nm =  2*np.pi/k_SP(E0,hbmu,hbgama)
#    lambda0_nm = np.real(lambda0_nm)

    
    labelx = 'a [nm]' 
    label1 = 'vs_a' 
    
    title = 'v = c/%i, $\mu$ = %.2f eV' %(int_v0,hbmu)
    
    listx = np.linspace(50, 1000,N)
    
    
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



listy_im_ana_tot = []
listy_omega_min_tot = []


listy_energy_tot = []
listy_energy_min_tot = []


if plot_vs_v == 1: 

    for n in list_n:


        listy_im_ana = []
        listy_omega_min = []
        
        listy_energy = []
        listy_energy_min = []
        
        
        for value in listx: 
    
            y_im_ana = omega_n_THz(1/value,a_nm,n)        
            listy_im_ana.append(y_im_ana)
            listy_energy.append(y_im_ana*(1e12*hb)*1e3)
        
            omega_min = omega_min_THz(1/value,a_nm,n)
            listy_omega_min.append(omega_min)
            listy_energy_min.append(omega_min*(1e12*hb)*1e3)
            
        
        listy_im_ana_tot.append(listy_im_ana)
        listy_energy_tot.append(listy_energy)


        listy_omega_min_tot.append(listy_omega_min)
        listy_energy_min_tot.append(listy_energy_min)




elif plot_vs_a == 1:       

    for n in list_n:


        listy_im_ana = []
        listy_omega_min = []
        
        listy_energy = []
        listy_energy_min = []
        
        for value in listx: 
            
    
            y_im_ana = omega_n_THz(int_v0,value,n)        
            listy_im_ana.append(y_im_ana)
            listy_energy.append(y_im_ana*(1e12*hb)*1e3)
            
            omega_min = omega_min_THz(int_v0,value,n)
            listy_omega_min.append(omega_min)
            listy_energy_min.append(omega_min*(1e12*hb)*1e3)
            
        
        listy_im_ana_tot.append(listy_im_ana)
        listy_energy_tot.append(listy_energy)


        listy_omega_min_tot.append(listy_omega_min)
        listy_energy_min_tot.append(listy_energy_min)
            
#%%         
    
colors = ['darkred','steelblue','coral','yellowgreen']    

graph(title,labelx,r'$\omega_n$ [THz]',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
k = 0
for modo in list_n:
    plt.plot(listx,listy_im_ana_tot[k],'.-',color = colors[modo],ms = ms, label = 'n = %i' %(modo))
    k = k + 1

if plot_vs_a == 1:
    plt.legend(loc = 'upper right',markerscale=1.5,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
else:
    plt.legend(loc = 'best',markerscale=1.5,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
#if plot_vs_c == 1:
#    plt.yscale('log')
os.chdir(path_save)
plt.savefig( 'omega_n' + label1 + '.png', format='png')   


            
#%%         
    
graph(title,labelx,r'$\hbar\omega^{min}_n$ [meV]',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
k = 0
for modo in list_n:
    plt.plot(listx,listy_energy_tot[k],'.-',color = colors[modo],ms = ms, label = 'n = %i' %(modo))
    plt.plot(listx,listy_energy_min_tot[k],'--',color = colors[modo],ms = ms)
    k = k + 1

if plot_vs_a == 1:
    plt.legend(loc = 'upper right',markerscale=1.5,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
else:
    plt.legend(loc = 'best',markerscale=1.5,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
#if plot_vs_c == 1:
#    plt.yscale('log')
os.chdir(path_save)
plt.savefig( 'energy_n' + label1 + '.png', format='png')        
        
#%%         
        
    