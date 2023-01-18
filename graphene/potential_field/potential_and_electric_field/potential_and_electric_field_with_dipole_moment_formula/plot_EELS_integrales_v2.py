
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

plot_vs_E = 1
plot_vs_c = 0

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/potential_field/potential_and_electric_field/potential_and_electric_field_with_dipole_moment_formula','')
path_save = path_basic + '/' + 'EELS_integrales_v2'
if not os.path.exists(path_save):
    print('Creating folder to save graphs')
    os.mkdir(path_save)

err = 'EELS_film.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from EELS_1dip_integrales_v2 import EELS_num_INT1, EELS_ana_INT1, EELS_num_INT2, EELS_ana_INT2,EELS_ana2_INT1, EELS_num_INT3, EELS_ana_INT3, EELS_ana_f_dir ,EELS_parallel_f_dir
except ModuleNotFoundError:
    print(err)

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

#v = c/int_v
#omega = 0.7*1e12

#v = c/int_v
#omega = 0.7*1e12


cota = 75 #nanometros
epsi1,epsi2 = 1,1
hbmu,hbgama = 0.3,0.0001
zp = 0.05
b = -0.01

energy0_pol = 44
omega0 = energy0_pol*1e-3/hb 
#R = 10 # 10 nm en unidades de micrometros
kappa_factor_omega0 = 0.1
kappa_r_factor= 0.5

L = 1000*1e-3
 
title1 = r'$\kappa$ = %.2f$\omega_0$, $\kappa_r$ = %.2f$\kappa$, $\hbar\omega_0$=%i meV' %(kappa_factor_omega0, kappa_r_factor, energy0_pol)     
#title2 = r'$\hbar\mu$ = %.2feV, $\hbar\gamma$ = %.4feV' %(hbmu,hbgama) 


title2 = r'$\hbar\mu$ = %.2feV, $\hbar\gamma$ = %.4f eV' %(hbmu,hbgama) 
title4 = r'$z_p$=%i nm, L = %i nm' %(zp*1e3, L*1e3)





N = 400
N = 150


if plot_vs_E == 1:

    int_v = 10    
    labelx = r'$\hbar\omega$ [meV]'
    
    title4 = title4 + ', v = c/%i' %(int_v)
    
    label1 = 'vs_E'
    #    listx = np.linspace(0.0001,2,N)
    listx = np.linspace(39,50,N)

else:
    E0 = 45
    labelx = r'v/c'
    
    title4 = title4 + ', $\hbar\omega$ = %i meV' %(E0)
    
    label1 = 'vs_c' 
    #    listx = np.linspace(0.0001,2,N)
    listx = np.linspace(0.01,0.4,N)

title = title1 + '\n' + title4

#%%

def function1_ana_re(energy0,int_v0):
    omegac0 = energy0*1e-3/aux 

    rta = EELS_ana_INT1(omegac0,epsi1,epsi2,hbmu,hbgama,int_v0,b,zp,omega0,kappa_factor_omega0,kappa_r_factor)
    
    return np.real(rta)


def function1_num_re(energy0,int_v0):
    omegac0 = energy0*1e-3/aux 

    rta = EELS_num_INT1(omegac0,epsi1,epsi2,hbmu,hbgama,int_v0,b,zp,omega0,kappa_factor_omega0,kappa_r_factor)
    
    return np.real(rta)

def function1_ana_re_v2(energy0,int_v0):
    omegac0 = energy0*1e-3/aux 

    rta = EELS_ana2_INT1(omegac0,epsi1,epsi2,hbmu,hbgama,int_v0,b,zp,omega0,kappa_factor_omega0,kappa_r_factor)
    
    return np.real(rta)


#%%

def function2_ana_re(energy0,int_v0):
    omegac0 = energy0*1e-3/aux 

    rta = EELS_ana_INT2(omegac0,epsi1,epsi2,hbmu,hbgama,int_v0,b,zp,omega0,kappa_factor_omega0,kappa_r_factor)
    
    return np.real(rta)


def function2_num_re(energy0,int_v0):
    omegac0 = energy0*1e-3/aux 

    rta = EELS_num_INT2(omegac0,epsi1,epsi2,hbmu,hbgama,int_v0,b,zp,omega0,kappa_factor_omega0,kappa_r_factor)
    
    return np.real(rta)

#%%
    
def function3_ana_re(energy0,int_v0):
    omegac0 = energy0*1e-3/aux 

    rta = EELS_ana_INT3(omegac0,epsi1,epsi2,hbmu,hbgama,int_v0,b,zp,omega0,kappa_factor_omega0,kappa_r_factor)
    
    return np.real(rta)


def function3_num_re(energy0,int_v0):
    omegac0 = energy0*1e-3/aux 

    rta = EELS_num_INT3(omegac0,epsi1,epsi2,hbmu,hbgama,int_v0,b,zp,omega0,kappa_factor_omega0,kappa_r_factor)
    
    return np.real(rta)

#%%
    
def function_tot_ana_re(energy0,int_v0):
    omegac0 = energy0*1e-3/aux 

    rta1 = EELS_ana2_INT1(omegac0,epsi1,epsi2,hbmu,hbgama,int_v0,b,zp,omega0,kappa_factor_omega0,kappa_r_factor)
    rta2 = EELS_ana_INT2(omegac0,epsi1,epsi2,hbmu,hbgama,int_v0,b,zp,omega0,kappa_factor_omega0,kappa_r_factor)
    rta3 = EELS_ana_INT3(omegac0,epsi1,epsi2,hbmu,hbgama,int_v0,b,zp,omega0,kappa_factor_omega0,kappa_r_factor)
    
    rta_dir = EELS_ana_f_dir(omegac0,epsi1,epsi2,hbmu,hbgama,int_v0,b,zp,omega0,kappa_factor_omega0,kappa_r_factor)
    
    rta_parallel = EELS_parallel_f_dir(omegac0,epsi1,epsi2,hbmu,hbgama,int_v0,b,zp,L)
    
    tot = rta1 +  rta2 + rta3 +  rta_dir + rta_parallel
    
    print(tot)
    
    return np.real(tot)


def function_tot_num_re(energy0,int_v0):
    omegac0 = energy0*1e-3/aux 

    rta1 = EELS_num_INT1(omegac0,epsi1,epsi2,hbmu,hbgama,int_v0,b,zp,omega0,kappa_factor_omega0,kappa_r_factor)
    rta2 = EELS_num_INT2(omegac0,epsi1,epsi2,hbmu,hbgama,int_v0,b,zp,omega0,kappa_factor_omega0,kappa_r_factor)
    rta3 = EELS_num_INT3(omegac0,epsi1,epsi2,hbmu,hbgama,int_v0,b,zp,omega0,kappa_factor_omega0,kappa_r_factor)
    
    rta_dir = EELS_ana_f_dir(omegac0,epsi1,epsi2,hbmu,hbgama,int_v0,b,zp,omega0,kappa_factor_omega0,kappa_r_factor)

    rta_parallel = EELS_parallel_f_dir(omegac0,epsi1,epsi2,hbmu,hbgama,int_v0,b,zp,L)
    
    tot = rta1 +  rta2 + rta3 +  rta_dir + rta_parallel
    
    return np.real(tot)

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

listy_re_ana_1 = []
listy_re_num_1 = []
listy_re_ana_1_v2 = []

listy_re_ana_2 = []
listy_re_num_2 = []

listy_re_ana_3 = []
listy_re_num_3 = []

listy_re_ana_tot = []
listy_re_num_tot = []


k = 0 
if plot_vs_E == 1:
    for value in listx: 
    
        y_re_ana_1 = function1_ana_re(value,int_v)  
        y_re_ana_1_v2 = function1_ana_re_v2(value,int_v)          
        y_re_num_1 = function1_num_re(value,int_v)          
        listy_re_num_1.append(y_re_num_1)
        listy_re_ana_1.append(y_re_ana_1)
        listy_re_ana_1_v2.append(y_re_ana_1_v2)
    

        y_re_ana_2 = function2_ana_re(value,int_v)  
        y_re_num_2 = function2_num_re(value,int_v)            
        listy_re_num_2.append(y_re_num_2)
        listy_re_ana_2.append(y_re_ana_2)

        y_re_ana_3 = function3_ana_re(value,int_v)  
        y_re_num_3 = function3_num_re(value,int_v)          
        listy_re_num_3.append(y_re_num_3)
        listy_re_ana_3.append(y_re_ana_3)
        
        
        y_re_tot_ana =  function_tot_ana_re(value,int_v)  
        y_re_tot_num = function_tot_num_re(value,int_v)  
        listy_re_ana_tot.append(y_re_tot_ana)
        listy_re_num_tot.append(y_re_tot_num)        
    
    
        print(k)
        k = k + 1


else:
    
    for value in listx: 
        
        value2 = 1/value
    
        y_re_ana_1 = function1_ana_re(E0,value2)  
        y_re_num_1 = function1_num_re(E0,value2)          
        y_re_ana_1_v2 = function1_ana_re_v2(E0,value2)    
        listy_re_num_1.append(y_re_num_1)
        listy_re_ana_1.append(y_re_ana_1)
        listy_re_ana_1_v2.append(y_re_ana_1_v2)   

        y_re_ana_2 = function2_ana_re(E0,value2)  
        y_re_num_2 = function2_num_re(E0,value2)          
        listy_re_num_2.append(y_re_num_2)
        listy_re_ana_2.append(y_re_ana_2)

        y_re_ana_3 = function3_ana_re(E0,value2)  
        y_re_num_3 = function3_num_re(E0,value2)          
        listy_re_num_3.append(y_re_num_3)
        listy_re_ana_3.append(y_re_ana_3)


        y_re_tot_ana =  function_tot_ana_re(E0,value2)  
        y_re_tot_num = function_tot_num_re(E0,value2)  
        listy_re_ana_tot.append(y_re_tot_ana)
        listy_re_num_tot.append(y_re_tot_num)     
    
        print(k)
        k = k + 1
    
#%%
##EELS en unidades de segundos --> agregar el hbar

graph(title,labelx,'Int 1 [eV$^{-1}$]',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(listx,np.array(listy_re_ana_1)/hb,'.-',ms = ms,color = 'purple',label ='analytical' )
plt.plot(listx,np.array(listy_re_ana_1_v2)/hb,'.-',ms = ms,color = 'darkred',label ='analytical 2' )
plt.plot(listx,np.array(listy_re_num_1)/hb,'.-',ms = ms,color = 'lightseagreen',label ='full numerical' )
plt.tight_layout()
plt.legend(loc = 'best',markerscale=1.5,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
#if plot_vs_c == 1:
#    plt.yscale('log')
os.chdir(path_save)
plt.savefig( 'EELS_int1' + label1 + '.png', format='png')   

graph(title,labelx,'Int 2 [eV$^{-1}$]',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(listx,np.array(listy_re_ana_2)/hb,'.-',ms = ms,color = 'purple',label ='analytical' )
plt.plot(listx,np.array(listy_re_num_2)/hb,'.-',ms = ms,color = 'lightseagreen',label ='full numerical' )
plt.tight_layout()
plt.legend(loc = 'best',markerscale=1.5,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
#if plot_vs_c == 1:
#    plt.yscale('log')
os.chdir(path_save)
plt.savefig( 'EELS_int2' + label1 + '.png', format='png')   

graph(title,labelx,'Int 3 [eV$^{-1}$]',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(listx,np.array(listy_re_ana_3)/hb,'.-',ms = ms,color = 'purple',label ='analytical' )
plt.plot(listx,np.array(listy_re_num_3)/hb,'.-',ms = ms,color = 'lightseagreen',label ='full numerical' )
plt.tight_layout()
plt.legend(loc = 'best',markerscale=1.5,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
#if plot_vs_c == 1:
#    plt.yscale('log')
os.chdir(path_save)
plt.savefig( 'EELS_int3' + label1 + '.png', format='png')   


graph(title,labelx,'EELS [eV$^{-1}$]',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(listx,np.array(listy_re_ana_tot)/hb,'.-',ms = ms,color = 'purple',label ='analytical' )
plt.plot(listx,np.array(listy_re_num_tot)/hb,'.-',ms = ms,color = 'lightseagreen',label ='full numerical' )
plt.tight_layout()
plt.legend(loc = 'best',markerscale=1.5,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
#if plot_vs_c == 1:
#    plt.yscale('log')
os.chdir(path_save)
plt.savefig( 'EELS_tot' + label1 + '.png', format='png')   

#%%

#
#    