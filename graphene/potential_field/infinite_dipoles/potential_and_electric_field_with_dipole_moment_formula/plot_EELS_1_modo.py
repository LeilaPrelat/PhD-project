
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
from scipy import special
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
path_constants =  path_basic.replace('/potential_field/infinite_dipoles/potential_and_electric_field_with_dipole_moment_formula','')
path_save = path_basic + '/' + 'EELS_1_modo'
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

def EELS(omegac,epsi1,epsi2,hbmu,hbgama,y,zp,int_v,a,b,n): ###en unidad de meV

    kx = omegac*int_v + 2*np.pi*n/a 
    
    E = omegac*aux

    cond = 4*pi*alfac*sigma_DL(E,hbmu,hbgama)
    # Rp = 2*epsi1/(epsi1 + epsi2)
    # alfa_p = 1j*(epsi1 + epsi2)/(cond*cte1)
    # kp = alfa_p*k1     

    px,py,pz  = dipole_moment_anav2(omegac,epsi1,epsi2,hbmu,hbgama,int_v,b,zp,omega0,kappa_factor_omega0,kappa_r_factor)
    
    cte_p = alfac*int_v
    
    px,py,pz = px*cte_p, py*cte_p, pz*cte_p
    
    
    alfa_p = 1j*(epsi1 + epsi2)/cond
    kp = alfa_p*omegac
    
    Rp = 2*epsi1/(epsi1 + epsi2)   
        
    cte = Rp/(4*(np.pi**3))
    
    term_den = np.sqrt(kp**2 - kx**2)
    
    term = 1j*kp*kx*(px*kx/term_den + py + 1j*pz*kp/term_den)*np.exp(1j*term_den*y)*np.exp(-kp*(2*zp-b))
    
    EmeV = E*1e3
    
    return term*cte*EmeV

def EELS_dir(omegac,epsi1,epsi2,hbmu,hbgama,y,zp,int_v,a,b,n): ###en unidad de meV

    kx = omegac*int_v + 2*np.pi*n/a 
    
    E = omegac*aux

    
    px,py,pz  = dipole_moment_anav2(omegac,epsi1,epsi2,hbmu,hbgama,int_v,b,zp,omega0,kappa_factor_omega0,kappa_r_factor)
    
    cte_p = alfac*int_v
    
    px,py,pz = px*cte_p, py*cte_p, pz*cte_p
    
    arg = kx*np.abs(b)
    K0 = special.kn(0,arg)
    K1 = special.kn(1,arg)
    

    cte = 1/(2*(np.pi**3))
    
    term = 1j*kx*kx*(-1j*px*K0 - pz*np.sign(b)*K1)
    
    EmeV = E*1e3
    
    return term*cte*EmeV

#%%


energy0_pol = 43
omega0 = energy0_pol*1e-3/hb 
#R = 10 # 10 nm en unidades de micrometros
kappa_factor_omega0 = 0.1
kappa_r_factor= 0.5

#v = c/int_v
#omega = 0.7*1e12
cota = 75 #nanometros
epsi1,epsi2 = 1,1
hbmu,hbgama = 0.3,0.0001

b = -0.01
#Ndipoles = 50
#list_xD = np.linspace(-0.01,0.01,Ndipoles)

y = 0

title1 = r'$\kappa$ = %.2f$\omega_0$, $\kappa_r$ = %.2f$\kappa$, $E_0$=%i meV' %(kappa_factor_omega0, kappa_r_factor, energy0_pol) 
title3 = r'y = %i nm, b = %inm' %(y*1e3,b*1e3)


N = 400
if plot_vs_a == 1:
    
    E0 = 36 # meV 
    int_v0 = 10
    zp0 = 0.05*1e3

    labelx = 'a [nm]'
    
    title2 = r'E = %i meV, v = c/%i, $z_p$ = %i nm' %(E0,int_v0,zp0)
    
    title = title1 + '\n' + title2 + '\n' +  title3  
    
    label1 = 'vs_a_E0%i' %(E0)
    listx = np.linspace(50,500,N)        ## en nanometros 
    

elif plot_vs_v == 1: 
    
    E0 = 36 # meV 
    a0 = 0.05*1e3   
    zp0 = 0.05*1e3
    
    
    labelx = 'v/c'  
    
    title2 = r'E = %i meV, a = %i nm, $z_p$ = %i nm' %(E0,a0,zp0)
    
    title = title1 + '\n' + title2 + '\n' +  title3 
    
    
    label1 = 'vs_v_E0%i' %(E0) 

    listx = np.linspace(0.01,0.99,N)
    
elif plot_vs_zp == 1: 
    
    E0 = 36 # meV 
    int_v0 = 10
    a0 = 0.05*1e3

    labelx = r'$z_p$ [nm]'
    
    title2 = r'E = %i meV, v = c/%i, a = %i nm' %(E0,int_v0,a0)
    
    title = title1 + '\n' + title2 + '\n' +  title3  
    
    label1 = 'vs_zp_E0%i' %(E0) 
    listx = np.linspace(50,800,N)        ## en nanometros 
    
elif plot_vs_E == 1:
    
    a0 = 0.05*1e3    
    int_v0 = 10
    zp0 = 0.2*1e3

    
    labelx = '$\hbar\omega$ [meV]'   
    
    title2 = r'a = %i nm, v = c/%i, $z_p$ = %i nm' %(a0,int_v0,zp0)
    
    title = title1 + '\n' + title2 + '\n' +  title3  
    
    label1 = 'vs_E_zp%inm' %(zp0) 
    
    
    listx = np.linspace(41,47,N)

#%%

def function_re(a_nano,zp_nano,int_v,energy0,n):
    omegac0 = energy0*1e-3/aux 
    a = a_nano*1e-3 
    zp = zp_nano*1e-3

    rta = EELS(omegac0,epsi1,epsi2,hbmu,hbgama,y,zp,int_v,a,b,n)
    # print(rta)
    return np.real(rta)

def function_re_dir(a_nano,zp_nano,int_v,energy0,n):
    omegac0 = energy0*1e-3/aux 
    a = a_nano*1e-3 
    zp = zp_nano*1e-3

    rta = EELS_dir(omegac0,epsi1,epsi2,hbmu,hbgama,y,zp,int_v,a,b,n)
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

labely1 = r'$\Gamma_{ind}$ [meV]'
labely2 = r'$\Gamma_{dir}$ [meV]'

colors = ['darkred','steelblue','coral','yellowgreen']

#%%
#
if plot_vs_a == 1:
    list_n = [1,2,3]
else:
    list_n = [0,1,2,3]
  
#list_n = [0,1,2,3]
  
listy_re_ana_tot = []    

listy_re_ana_tot_dir = []    

if plot_vs_a == 1: 
    
    for n in list_n:

        listy_re_ana = []    
        
        listy_re_ana_dir = [] 

        for value in listx: 
    
            y_re_ana = function_re(value,zp0,int_v0,E0,n)    
            listy_re_ana.append(y_re_ana)
            
            y_re_ana_dir = function_re_dir(value,zp0,int_v0,E0,n)    
            listy_re_ana_dir.append(y_re_ana_dir)
            
            
        listy_re_ana_tot.append(listy_re_ana)
        
        listy_re_ana_tot_dir.append(listy_re_ana_dir)



elif plot_vs_v == 1: 
    
    for n in list_n:

        listy_re_ana = []    
        
        listy_re_ana_dir = [] 

        for value in listx: 
    
            y_re_ana = function_re(a0,zp0,value,E0,n)    
            listy_re_ana.append(y_re_ana)
        
            y_re_ana_dir = function_re_dir(a0,zp0,value,E0,n)    
            listy_re_ana_dir.append(y_re_ana_dir)
            
            
        listy_re_ana_tot.append(listy_re_ana)
        
        listy_re_ana_tot_dir.append(listy_re_ana_dir)


elif plot_vs_zp == 1: 
    
    for n in list_n:

        listy_re_ana = []    
        
        listy_re_ana_dir = [] 

        for value in listx: 
    
            y_re_ana = function_re(a0,value,int_v0,E0,n)    
            listy_re_ana.append(y_re_ana)
            
            y_re_ana_dir = function_re_dir(a0,value,int_v0,E0,n)    
            listy_re_ana_dir.append(y_re_ana_dir)
            
            
        listy_re_ana_tot.append(listy_re_ana)
        
        listy_re_ana_tot_dir.append(listy_re_ana_dir)



elif plot_vs_E == 1: 

    for n in list_n:
        
        listy_re_ana = [] 
        
        listy_re_ana_dir = [] 
    
        for value in listx: 
    
            y_re_ana = function_re(a0,zp0,int_v0,value,n)     
            listy_re_ana.append(y_re_ana)
    
            y_re_ana_dir = function_re_dir(a0,zp0,int_v0,value,n)    
            listy_re_ana_dir.append(y_re_ana_dir)
            
            
        listy_re_ana_tot.append(listy_re_ana)
        
        listy_re_ana_tot_dir.append(listy_re_ana_dir)

        
#%%       
graph(title,labelx,labely1,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
k = 0
for modo in list_n:
    if plot_vs_zp == 1 and modo == 0:
        plt.plot(listx,np.array(listy_re_ana_tot[k])*10,'.-',color = colors[modo],ms = ms, label = 'n = %i (x 10)' %(modo))
    else:
        plt.plot(listx,listy_re_ana_tot[k],'.-',color = colors[modo],ms = ms, label = 'n = %i' %(modo))
    k = k + 1
    
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
os.chdir(path_save)
plt.savefig( 'EELS_ind' + label1 + '.png', format='png')   
   


graph(title,labelx,labely2,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
k = 0
for modo in list_n:
    
    if plot_vs_zp == 1 and modo == 0:
        plt.plot(listx,np.array(listy_re_ana_tot_dir[k])*10,'.-',color = colors[modo],ms = ms, label = 'n = %i (x 10)' %(modo))
    else:
        plt.plot(listx,listy_re_ana_tot_dir[k],'.-',color = colors[modo],ms = ms, label = 'n = %i' %(modo))
    k = k + 1
    
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
os.chdir(path_save)
plt.savefig( 'EELS_dir' + label1 + '.png', format='png')   
 

#%%


  