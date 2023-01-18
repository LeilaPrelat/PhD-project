
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

plot_vs_R = 1
plot_vs_E = 0

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/potential_field/potential_and_electric_field/potential_and_electric_field_with_dipole_moment_formula','')
path_save = path_basic + '/' + 'potential_final_version'
if not os.path.exists(path_save):
    print('Creating folder to save graphs')
    os.mkdir(path_save)

err = 'fieldE_direct_numerical.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from potential_final_version import potential_final_ana, potential_final_pole_aprox, potential_final_num
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
int_v = 10
#v = c/int_v
#omega = 0.7*1e12
cota = 75 #nanometros
epsi1,epsi2 = 1,1
hbmu,hbgama = 0.3,0.0001
zp = 0.05
b = -0.01

phi = 0
z = zp



#omega0THz = 65
#omega0 = omega0THz*1e12 
#energy0_pol = omega0*hb
energy0_pol = 25
omega0 = energy0_pol*1e-3/hb 
#R = 10 # 10 nm en unidades de micrometros
kappa_factor_omega0 = 0.1
kappa_r_factor= 0.5
 
title1 = r'v = c/%i, $\kappa$ = %.2f$\omega_0$, $\kappa_r$ = %.2f$\kappa$, $\hbar\omega_0$=%i meV, $z_p$=%i nm' %(int_v, kappa_factor_omega0, kappa_r_factor, energy0_pol,zp*1e3)     
#title2 = r'$\hbar\mu$ = %.2feV, $\hbar\gamma$ = %.4feV' %(hbmu,hbgama) 
#title3 = r'$z_p$=%inm, px=%i, py=%i, pz=%i' %(zp*1e3,px,py,pz)
title4 = r'b = %i nm, $\phi$=%i, z = %i nm' %(b*1e3,phi,z*1e3)
labelp = r'_E0_%i' %(energy0_pol)

N = 400
if plot_vs_R == 1: 
    E0 = 14.36 # meV
    E0 = 35 # meV
    labelx = 'R [nm]'  
    title4 = title4 + ', ' + '$\hbar\omega$ = %i meV' %(E0)
    label1 = 'vs_R' + labelp
    listx = np.linspace(100,4000,N)
else:
    R0 = 1000 #nanos
    # z0 = 0.06*1e3
    labelx = r'$\hbar\omega$ [meV]'   
    title4 = title4 + ', ' + 'R = %i nm' %(R0)
    label1 = 'vs_E' + labelp
#    listx = np.linspace(0.0001,2,N)
    listx = np.linspace(20,55,N)
  
title = title1  + '\n'  + title4 

#%%

def function_real_ana(R_nano,energy0):
    omegac0 = energy0*1e-3/aux 
    R = R_nano*1e-3 

    rta = potential_final_ana(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,b,zp,phi,R,z,omega0,kappa_factor_omega0,kappa_r_factor)
    
    return rta.real

def function_imag_ana(R_nano,energy0):
    omegac0 = energy0*1e-3/aux 
    R = R_nano*1e-3 
    
    rta = potential_final_ana(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,b,zp,phi,R,z,omega0,kappa_factor_omega0,kappa_r_factor)
    return rta.imag

#%%

def function_real_num(R_nano,energy0):
    omegac0 = energy0*1e-3/aux 
    R = R_nano*1e-3

    rta = potential_final_num(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,b,zp,phi,R,z,omega0,kappa_factor_omega0,kappa_r_factor)
    
    return rta.real

def function_imag_num(R_nano,energy0):
    omegac0 = energy0*1e-3/aux 
    R = R_nano*1e-3
    
    rta = potential_final_num(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,b,zp,phi,R,z,omega0,kappa_factor_omega0,kappa_r_factor)
    return rta.imag


#%%


def function_real_pole_aprox(R_nano,energy0):
    omegac0 = energy0*1e-3/aux 
    R = R_nano*1e-3

    rta = potential_final_pole_aprox(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,b,zp,phi,R,z,omega0,kappa_factor_omega0,kappa_r_factor)
    
    return rta.real

def function_imag_pole_aprox(R_nano,energy0):
    omegac0 = energy0*1e-3/aux 
    R = R_nano*1e-3
    
    rta = potential_final_pole_aprox(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,b,zp,phi,R,z,omega0,kappa_factor_omega0,kappa_r_factor)
    return rta.imag


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

if plot_vs_R == 1: 
    listy_re_ana = []
    listy_im_ana = []

    listy_re_num = []
    listy_im_num = []
   
    listy_re_pole_aprox = []
    listy_im_pole_aprox = []
   
    
    
    for value in listx: 

        y_re_ana = function_real_ana(value,E0)
        y_im_ana = function_imag_ana(value,E0)        


        y_re_num = function_real_num(value,E0)
        y_im_num = function_imag_num(value,E0)  
        
        
        y_re_pole_aprox = function_real_pole_aprox(value,E0)
        y_im_pole_aprox = function_imag_pole_aprox(value,E0)  
        
        
        listy_re_ana.append(y_re_ana)
        listy_im_ana.append(y_im_ana)
        
        listy_re_num.append(y_re_num)
        listy_im_num.append(y_im_num)   
        
        listy_re_pole_aprox.append(y_re_pole_aprox)
        listy_im_pole_aprox.append(y_im_pole_aprox)   
        
    
#    graph(title,labelx,'Re($\phi$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#    plt.plot(listx,listy_re_ana,'.',ms = ms,color = 'purple',label = 'PP analytical')
#    plt.plot(listx,listy_re_pole_aprox,'.',ms = ms,color = 'darkred',label = 'PP numerical')
#    plt.plot(listx,listy_re_num,'.-',ms = ms,color = 'lightseagreen',label = 'full numerical')
#    plt.legend(loc = 'lower left',markerscale=mk,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
#    plt.tight_layout()
#    os.chdir(path_save)
#    plt.savefig( 'Re_phi' + label1 + '.png', format='png')   
#
#    graph(title,labelx,'Im($\phi$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#    plt.plot(listx,listy_im_ana,'.',ms = ms,color = 'purple',label = 'PP analytical')
#    plt.plot(listx,listy_im_pole_aprox,'.',ms = ms,color = 'darkred',label = 'PP numerical')
#    plt.plot(listx,listy_im_num,'.-',ms = ms,color = 'lightseagreen',label = 'full numerical')
#    plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
#    plt.tight_layout()
#    os.chdir(path_save)
#    plt.savefig( 'Im_phi' + label1 + '.png', format='png')   
#    
if plot_vs_E == 1: 
    listy_re_ana = []
    listy_im_ana = []

    listy_re_num = []
    listy_im_num = []
   
    listy_re_pole_aprox = []
    listy_im_pole_aprox = []
   
    
    
    for value in listx: 

        y_re_ana = function_real_ana(R0,value)
        y_im_ana = function_imag_ana(R0,value)        


        y_re_num = function_real_num(R0,value)
        y_im_num = function_imag_num(R0,value)  
        
        
        y_re_pole_aprox = function_real_pole_aprox(R0,value)
        y_im_pole_aprox = function_imag_pole_aprox(R0,value)  
        
        
        listy_re_ana.append(y_re_ana)
        listy_im_ana.append(y_im_ana)
        
        listy_re_num.append(y_re_num)
        listy_im_num.append(y_im_num)   
        
        listy_re_pole_aprox.append(y_re_pole_aprox)
        listy_im_pole_aprox.append(y_im_pole_aprox)   
        
    
#%%
        
graph(title,labelx,'Re($\phi$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(listx,listy_re_ana,'.',ms = ms,color = 'purple',label = 'PP analytical')
plt.plot(listx,listy_re_pole_aprox,'.',ms = ms,color = 'darkred',label = 'PP numerical')
plt.plot(listx,listy_re_num,'.-',ms = ms,color = 'lightseagreen',label = 'full numerical')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
os.chdir(path_save)
plt.savefig( 'Re_phi' + label1 + '.png', format='png')   

graph(title,labelx,'Im($\phi$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(listx,listy_im_ana,'.',ms = ms,color = 'purple',label = 'PP analytical')
plt.plot(listx,listy_im_pole_aprox,'.',ms = ms,color = 'darkred',label = 'PP numerical')
plt.plot(listx,listy_im_num,'.-',ms = ms,color = 'lightseagreen',label = 'full numerical')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
os.chdir(path_save)
plt.savefig( 'Im_phi' + label1 + '.png', format='png')   
    
#%%
