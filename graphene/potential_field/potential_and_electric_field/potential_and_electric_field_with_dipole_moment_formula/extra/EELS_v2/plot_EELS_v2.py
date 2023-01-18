
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
path_save = path_basic + '/' + 'EELS_v2'
if not os.path.exists(path_save):
    print('Creating folder to save graphs')
    os.mkdir(path_save)

err = 'EELS_film.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from EELS_1dip_v2 import EELS_ana_f_ind1,EELS_PP_f_ind,EELS_ana_f_dir
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
 
title1 = r'$\kappa$ = %.2f$\omega_0$, $\kappa_r$ = %.2f$\kappa$, $\hbar\omega_0$=%i meV' %(kappa_factor_omega0, kappa_r_factor, energy0_pol)     
#title2 = r'$\hbar\mu$ = %.2feV, $\hbar\gamma$ = %.4feV' %(hbmu,hbgama) 


title2 = r'$\hbar\mu$ = %.2feV, $\hbar\gamma$ = %.4f eV' %(hbmu,hbgama) 
title4 = r'$z_p$=%i nm' %(zp*1e3)



N = 400
N = 150


if plot_vs_E == 1:

    int_v = 10    
    labelx = r'$\hbar\omega$ [meV]'
    
    title4 = title4 + ', v = c/%i' %(int_v)
    
    label1 = 'vs_E'
    #    listx = np.linspace(0.0001,2,N)
    listx = np.linspace(42,54,N)

else:
    E0 = 500
    labelx = r'v/c'
    
    title4 = title4 + ', $\hbar\omega$ = %i meV' %(E0)
    
    label1 = 'vs_c' 
    #    listx = np.linspace(0.0001,2,N)
    listx = np.linspace(0.01,0.99,N)

title = title1 + '\n' + title4

#%%

def function_ana_re(energy0,int_v0):
    omegac0 = energy0*1e-3/aux 

    rta = EELS_ana_f_ind1(omegac0,epsi1,epsi2,hbmu,hbgama,int_v0,b,zp,omega0,kappa_factor_omega0,kappa_r_factor)
    
    return np.real(rta)


def function_num_re(energy0,int_v0):
    omegac0 = energy0*1e-3/aux 

    rta = EELS_PP_f_ind(omegac0,epsi1,epsi2,hbmu,hbgama,int_v0,b,zp,omega0,kappa_factor_omega0,kappa_r_factor)
    
    return np.real(rta)



def function_ana_re_dir(energy0,int_v0):
    omegac0 = energy0*1e-3/aux 

    rta = EELS_ana_f_dir(omegac0,epsi1,epsi2,hbmu,hbgama,int_v0,b,zp,omega0,kappa_factor_omega0,kappa_r_factor)
    
    return np.real(rta)


#%%

def function_ana_im(energy0,int_v0):
    omegac0 = energy0*1e-3/aux 

    rta = EELS_ana_f_ind1(omegac0,epsi1,epsi2,hbmu,hbgama,int_v0,b,zp,omega0,kappa_factor_omega0,kappa_r_factor)
    
    return np.imag(rta)

def function_num_im(energy0,int_v0):
    omegac0 = energy0*1e-3/aux 

    rta = EELS_PP_f_ind(omegac0,epsi1,epsi2,hbmu,hbgama,int_v0,b,zp,omega0,kappa_factor_omega0,kappa_r_factor)
    
    return np.imag(rta)



def function_ana_im_dir(energy0,int_v0):
    omegac0 = energy0*1e-3/aux 

    rta = EELS_ana_f_dir(omegac0,epsi1,epsi2,hbmu,hbgama,int_v0,b,zp,omega0,kappa_factor_omega0,kappa_r_factor)
    
    return np.imag(rta)


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

listy_re_ana_ind = []
listy_re_num_ind = []



listy_re_ana_tot = []


k = 0 
if plot_vs_E == 1:
    for value in listx: 
    
        y_re_ana_ind = function_ana_re(value,int_v)  
        y_re_ana_dir = function_ana_re_dir(value,int_v)  
        y_re_num = function_num_re(value,int_v)        
        
    
        
        listy_re_num_ind.append(y_re_num)
        listy_re_ana_ind.append(y_re_ana_ind)
    
        
        listy_re_ana_tot.append(y_re_ana_dir + y_re_ana_ind)
    
        print(k)
        k = k + 1
#        y_im_ana = function_ana_im(value,int_v)        
#        y_im_num = function_num_im(value,int_v)        
#        
#        listy_im_ana.append(y_im_ana)
#        listy_im_num.append(y_im_num)

else:
    
    for value in listx: 
        
            value2 = 1/value
        
            y_re_ana_ind = function_ana_re(int_v,value2)  
            y_re_ana_dir = function_ana_re_dir(int_v,value2)  
            y_re_num = function_num_re(int_v,value2)        
            
        
            
            listy_re_num_ind.append(y_re_num)
            listy_re_ana_ind.append(y_re_ana_ind)       
            
            listy_re_ana_tot.append(y_re_ana_dir + y_re_ana_ind)
        
            print(k)
            k = k + 1
        
#            y_im_ana = function_ana_im(E0,value2)        
#            y_im_num = function_num_im(E0,value2)        
#            
#            listy_im_ana.append(y_im_ana)
#            listy_im_num.append(y_im_num)
#    
    
#%%

x1 = listx[22]
x2 = listx[23]

xmedio = np.mean([x1,x2])
mini = np.min(listy_re_ana_tot)
maxi = np.max(listy_re_ana_tot)
eje_extra = np.linspace(mini,maxi,10)

graph(title,labelx,'$\Gamma_{ind}$',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(listx,listy_re_ana_ind,'.-',ms = ms,color = 'purple')
plt.plot(listx,listy_re_num_ind,'.-',ms = ms,color = 'lightseagreen',label ='full numerical' )
plt.tight_layout()
plt.legend(loc = 'best',markerscale=1.5,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
#if plot_vs_c == 1:
#    plt.yscale('log')
os.chdir(path_save)
plt.savefig( 'EELS_ind' + label1 + '.png', format='png')   

graph(title,labelx,'$\Gamma = \Gamma_{dir} + \Gamma_{ind}$',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.legend(loc = 'best',markerscale=1.5,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.plot(listx,listy_re_ana_tot,'.-',ms = ms,color = 'purple')
#plt.plot(listx,listy_re_num,'.-',ms = ms,color = 'lightseagreen',label ='full numerical' )
plt.tight_layout()
#if plot_vs_c == 1:
#    plt.yscale('log')
os.chdir(path_save)
plt.savefig( 'EELS_tot' + label1 + '.png', format='png')   


#%%






#if plot_vs_E == 1:
#    
#    
#    hspace = 0.2
#    wspace = 0
#    fig, axs = plt.subplots(2,1, sharex=False, facecolor='w', figsize = (5.2,3.9))
#    plt.subplots_adjust(hspace =hspace,wspace = wspace)
#    
#    fig.suptitle(title,fontsize=int(tamtitle*0.9))
#    axs[0].plot(listx[0:int(N/2)],listy_re_ana[0:int(N/2)],'.-',ms = ms,color = 'purple')
#    axs[1].plot(listx[int(N/2):N],listy_re_ana[int(N/2):N],'.-',ms = ms,color = 'purple')
#    axs[1].minorticks_on()
#    axs[1].tick_params(labelsize = tamnum,pad = pad)
#    axs[0].tick_params(labelsize = tamnum,pad = pad)
#    axs[1].set_ylabel('$\Gamma_{ind}$',fontsize = tamletra,labelpad=labelpady) 
#    axs[0].set_ylabel('$\Gamma_{ind}$',fontsize = tamletra,labelpad=labelpady) 
#    axs[1].set_xlabel(labelx,fontsize = tamletra,labelpad=labelpadx)
#    #axs[1].title.set_text(title)
##    fig.tight_layout()
#    fig.savefig( 'EELS' + label1 + '.png', format='png')   


#
#graph(title,labelx,'$\Gamma_{ind}$',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#plt.plot(listx,listy_re_ana,'.-',ms = ms,color = 'purple',label = 'analytical')
#plt.plot(listx,listy_re_num,'.-',ms = ms,color = 'lightseagreen',label ='numerical' )
#plt.legend(loc = 'best',markerscale=1.5,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
#plt.tight_layout()
##if plot_vs_c == 1:
##    plt.yscale('log')
#
#
#graph(title,labelx,'$\Gamma_{ind}$',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#plt.plot(listx,listy_im_ana,'.-',ms = ms,color = 'purple',label = 'analytical')
#plt.plot(listx,listy_im_num,'.-',ms = ms,color = 'lightseagreen',label ='numerical' )
#plt.legend(loc = 'best',markerscale=1.5,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
#plt.tight_layout()
##if plot_vs_c == 1:
##    plt.yscale('log')
#
#    