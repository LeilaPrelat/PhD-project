
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

plot_vs_z = 0
plot_vs_E = 1

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/potential_field/potential_and_electric_field/potential_and_electric_field_simple_version/potential_1_dipolo_analytical_vs_num','')
path_save = path_basic + '/' + 'potential_dif_terms'
if not os.path.exists(path_save):
    print('Creating folder to save graphs')
    os.mkdir(path_save)

err = 'fieldE_direct_numerical.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from potential_dif_terms import electric_potential_ana1, electric_potential_num1, electric_potential_ana2, electric_potential_num2
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
y0,x0 = 0.01,0.01

xD = 0
yD = 0
zD = 0

px = 1
py = 0
pz = 0

title2 = r'$\hbar\mu$ = %.2feV, $\hbar\gamma$ = %.4feV' %(hbmu,hbgama) 
title4 = r'$z_p$=%inm, px = %i, py = %i, pz = %i, xD = %inm, yD = %inm' %(zp*1e3,px,py,pz,xD*1e3,yD*1e3)
title5 = r'zD = %inm, x = %inm, y = %inm' %(zD*1e3,x0*1e3,y0*1e3)
labelp = r'_px%i_py%i_pz%i' %(px,py,pz)

N = 100
if plot_vs_z == 1: 
    E0 = 25 # meV 
    labelx = 'z [nm]'  
    title5 = title5 + ', ' + 'E = %.2f meV' %(E0)
    label1 = 'vs_z' + labelp
    listx = np.linspace(zp*1e3,-2*zp*1e3,N)
else:
    z0 = zp*1e3
    # z0 = 0.06*1e3
    labelx = 'E [meV]'   
    title5 = title5 + ', ' + 'z = %inm' %(z0)
    label1 = 'vs_E' + labelp
    listx = np.linspace(0.0001,2,N)
    listx = np.linspace(1,50,N)
  
title = title2 + '\n' + title4 + '\n'  + title5

def function_num1(z_nano,energy0):
    omegac0 = energy0*1e-3/aux 
    z = z_nano*1e-3 

    rta = electric_potential_num1(omegac0,epsi1,epsi2,hbmu,hbgama,x0,y0,z,zp,px,py,pz)
    # print(rta)
    return np.abs(rta)

def function_num2(z_nano,energy0):
    omegac0 = energy0*1e-3/aux 
    z = z_nano*1e-3 
    
    rta = electric_potential_num2(omegac0,epsi1,epsi2,hbmu,hbgama,x0,y0,z,zp,px,py,pz)
    return np.abs(rta)


def function_ana1(z_nano,energy0):
    omegac0 = energy0*1e-3/aux 
    z = z_nano*1e-3 

    rta = electric_potential_ana1(omegac0,epsi1,epsi2,hbmu,hbgama,x0,y0,z,zp,px,py,pz)
    
    return np.abs(rta)

def function_ana2(z_nano,energy0):
    omegac0 = energy0*1e-3/aux 
    z = z_nano*1e-3 
    
    rta = electric_potential_ana2(omegac0,epsi1,epsi2,hbmu,hbgama,x0,y0,z,zp,px,py,pz)
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

if plot_vs_z == 1: 
    listy_re_ana = []
    listy_im_ana = []
    listy_re_num = []
    listy_im_num = []
    
    for value in listx: 

        y_re_ana = function_real_ana(value,E0)
        y_im_ana = function_imag_ana(value,E0)        

        y_re_num = function_real_num(value,E0)
        y_im_num = function_imag_num(value,E0)
        
        listy_re_ana.append(y_re_ana)
        listy_im_ana.append(y_im_ana)
        listy_re_num.append(y_re_num)
        listy_im_num.append(y_im_num)    

    
    graph(title,labelx,'Re($\phi$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    plt.plot(listx,listy_re_ana,'.',ms = ms,color = 'purple',label = 'pole aprox')
    plt.plot(listx,listy_re_num,'.-',ms = ms,color = 'lightseagreen',label = 'sin aprox')
    plt.legend(loc = 'best',markerscale=1,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'Re_phi' + label1 + '.png', format='png')   

    graph(title,labelx,'Im($\phi$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    plt.plot(listx,listy_im_ana,'.',ms = ms,color = 'purple',label = 'pole aprox')
    plt.plot(listx,listy_im_num,'.-',ms = ms,color = 'lightseagreen',label = 'sin aprox')
    plt.legend(loc = 'best',markerscale=1.5,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'Im_phi' + label1 + '.png', format='png')   
    
if plot_vs_E == 1: 
    listy_ana1 = []
    listy_ana2 = []
    listy_num1 = []
    listy_num2 = []
    
    for value in listx: 

        y_ana1 = function_ana1(z0,value)
        y_ana2 = function_ana2(z0,value)        

        y_num1 = function_num1(z0,value)
        y_num2 = function_num2(z0,value)
        
        listy_ana1.append(y_ana1)
        listy_ana2.append(y_ana2)
        listy_num1.append(y_num1)
        listy_num2.append(y_num2)   

 
    graph(title,labelx,r'$\phi$ 1 term',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    plt.plot(listx,listy_ana1,'.',ms = ms,color = 'purple',label = 'pole aprox')
    plt.plot(listx,listy_num1,'.-',ms = ms,color = 'lightseagreen',label = 'sin aprox')
    plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'phi_term1' + label1 + '.png', format='png')   
   
    graph(title,labelx,r'$\phi$ 2 term',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    plt.plot(listx,listy_ana2,'.',ms = ms,color = 'purple',label = 'pole aprox')
    plt.plot(listx,listy_num2,'.-',ms = ms,color = 'lightseagreen',label = 'sin aprox')
    plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'phi_term2' + label1 + '.png', format='png')   
    

#%%
