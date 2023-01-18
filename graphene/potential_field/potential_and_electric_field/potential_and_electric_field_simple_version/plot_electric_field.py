
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
path_constants =  path_basic.replace('/potential_field/potential_and_electric_field/potential_and_electric_field_simple_version','')
path_save = path_basic + '/' + 'electric_field'
if not os.path.exists(path_save):
    print('Creating folder to save graphs')
    os.mkdir(path_save)

err = 'fieldE_direct_numerical.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from electric_field_ana import electric_field_ana_f
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
cota = 75 #nanometros
epsi1,epsi2 = 1,1
hbmu,hbgama = 0.3,0.0001
zp = 0.05
phi0 = 0
z0 = zp

px = 1
py = 0
pz = 0

title2 = r'$\hbar\mu$ = %.2feV, $\hbar\gamma$ = %.4feV' %(hbmu,hbgama) 
title4 = r'$z_p$=%inm, px = %i, py = %i, pz = %i' %(zp*1e3,px,py,pz)
title5 = r'z = %inm, $\varphi$ = %i' %(z0*1e3,phi0)
labelp = r'_px%i_py%i_pz%i' %(px,py,pz)

N = 300
if plot_vs_R == 1: 
    E0 = 60 # meV 
    labelx = 'R [nm]'  
    title5 = title5 + ', ' + 'E = %.2f meV' %(E0)
    label1 = 'vs_R_E%i' %(E0) + labelp
    listx = np.linspace(500,6000,N)
    
else:
    R0 = 5000 #nanos
    
    # z0 = 0.06*1e3
    labelx = 'E [meV]'   
    title5 = title5 + ', ' + 'R = %inm' %(R0)
    label1 = 'vs_z%i' %(z0*1e3) + labelp
    listx = np.linspace(0.0001,2,N)
    listx = np.linspace(20,65,N)
  
title = title2 + '\n' + title4 + '\n'  + title5

def function_real_ana(R_nano,energy0):
    omegac0 = energy0*1e-3/aux 
    R = R_nano*1e-3 

    rta = electric_field_ana_f(omegac0,epsi1,epsi2,hbmu,hbgama,phi0,R,z0,zp,px,py,pz)
    # print(rta)`
    return rta.real

def function_imag_ana(R_nano,energy0):
    omegac0 = energy0*1e-3/aux 
    R = R_nano*1e-3 

    rta = electric_field_ana_f(omegac0,epsi1,epsi2,hbmu,hbgama,phi0,R,z0,zp,px,py,pz)
    
    return rta.imag

def function_abs_ana(R_nano,energy0):
    omegac0 = energy0*1e-3/aux 
    R = R_nano*1e-3 

    rta = electric_field_ana_f(omegac0,epsi1,epsi2,hbmu,hbgama,phi0,R,z0,zp,px,py,pz)
    
    return np.abs(rta)

#%%
    
tamfig = (4.5,3.5)
tamlegend = 10
tamletra = 11
tamtitle = 10
tamnum = 9
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
    listy_abs_ana = []
    
    for value in listx: 

        y_re_ana = function_real_ana(value,E0)
        y_im_ana = function_imag_ana(value,E0)        
        y_abs_ana = function_abs_ana(value,E0)
        
        listy_re_ana.append(y_re_ana)
        listy_im_ana.append(y_im_ana)
        listy_abs_ana.append(y_abs_ana)    
    
    graph(title,labelx,'Re($E_R$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    plt.plot(listx,listy_re_ana,'.',ms = ms,color = 'purple')
    plt.legend(loc = 'best',markerscale=1,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'Re_ER' + label1 + '.png', format='png')   

    graph(title,labelx,'Im($E_R$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    plt.plot(listx,listy_im_ana,'.',ms = ms,color = 'purple')
    plt.legend(loc = 'best',markerscale=1.5,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'Im_ER' + label1 + '.png', format='png')   


    graph(title,labelx,'|$E_R$|',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    plt.plot(listx,listy_abs_ana,'.',ms = ms,color = 'purple')
    plt.legend(loc = 'best',markerscale=1.5,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'abs_ER' + label1 + '.png', format='png')   
    
if plot_vs_E == 1: 
    listy_re_ana = []
    listy_im_ana = []
    
    for value in listx: 

        y_re_ana = function_real_ana(R0,value)
        y_im_ana = function_imag_ana(R0,value)        
        
        listy_re_ana.append(y_re_ana)
        listy_im_ana.append(y_im_ana)
        
    graph(title,labelx,'Re($E_R$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    plt.plot(listx,listy_re_ana,'.-',ms = ms,color = 'purple')
    plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'Re_ER' + label1 + '.png', format='png')   
   
    graph(title,labelx,'Im($E_R$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    plt.plot(listx,listy_im_ana,'.-',ms = ms,color = 'purple')
    plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'Im_ER' + label1 + '.png', format='png')   
    
#%%
