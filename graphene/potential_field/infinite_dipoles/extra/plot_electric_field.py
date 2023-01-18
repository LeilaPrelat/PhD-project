
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

plot_vs_y = 0
plot_vs_E = 1

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/potential_field/infinite_dipoles','')
path_save = path_basic + '/' + 'many_dipoles_E'
if not os.path.exists(path_save):
    print('Creating folder to save graphs')
    os.mkdir(path_save)

err = 'many_potential_integrals.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from electric_field_many_dipoles import Ex_many_dipoles_num, Ex_many_dipoles_ana, Ex_many_dipoles_pole_aprox
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
b = -0.01

x0,z0 = 0.001, zp 

#Ndipoles = 50
#list_xD = np.linspace(-0.01,0.01,Ndipoles)
a = 0.05

px = 1
py = 0
pz = 0

int_v = 10

title2 = r'$\hbar\mu$ = %.2feV, $\hbar\gamma$ = %.4feV, v = c/%i' %(hbmu,hbgama,int_v) 
title3 = r'$z_p$=%inm, px=%i, py=%i, pz=%i' %(zp*1e3,px,py,pz)
title4 = r'b = %inm, z=%inm, x=%inm' %(b*1e3,z0*1e3,x0*1e3)

labelp = r'_px%i_py%i_pz%i_a%inm' %(px,py,pz,a*1e3)

N = 250
if plot_vs_y == 1: 
    E0 = 30 # meV 
    labelx = 'y [nm]'  
    title = title2 + '\n' +  title3 + '\n' + title4 + ', E = %.2f meV' %(E0)
    label1 = 'vs_y' 
    listx = np.linspace(10,500,N)
else:
    
    y0 = 0.05*1e3
#    y0 = 0
    # z0 = 0.06*1e3
    labelx = 'E [meV]'   
    title = title2 + '\n' + title3 + '\n' + title4  + ', y = %inm' %(y0)
    label1 = 'vs_E' 
    listx = np.linspace(80,20,N)

#%%

def function_real_num(y_nano,energy0):
    omegac0 = energy0*1e-3/aux 
    y = y_nano*1e-3 

    rta = Ex_many_dipoles_num(omegac0,epsi1,epsi2,hbmu,hbgama,x0,y,z0,a,zp,int_v,px,py,pz)
    # print(rta)
    return rta.real

def function_imag_num(y_nano,energy0):
    omegac0 = energy0*1e-3/aux 
    y = y_nano*1e-3 
    
    rta = Ex_many_dipoles_num(omegac0,epsi1,epsi2,hbmu,hbgama,x0,y,z0,a,zp,int_v,px,py,pz)
    return rta.imag

#%%
    
def function_real_ana(y_nano,energy0):
    omegac0 = energy0*1e-3/aux 
    y = y_nano*1e-3 

    rta = Ex_many_dipoles_ana(omegac0,epsi1,epsi2,hbmu,hbgama,x0,y,z0,a,zp,int_v,px,py,pz)
    # print(rta)
    return rta.real

def function_imag_ana(y_nano,energy0):
    omegac0 = energy0*1e-3/aux 
    y = y_nano*1e-3 
    
    rta = Ex_many_dipoles_ana(omegac0,epsi1,epsi2,hbmu,hbgama,x0,y,z0,a,zp,int_v,px,py,pz)
    return rta.imag

#%%
    
def function_real_pole_aprox(y_nano,energy0):
    omegac0 = energy0*1e-3/aux 
    y = y_nano*1e-3 

    rta = Ex_many_dipoles_pole_aprox(omegac0,epsi1,epsi2,hbmu,hbgama,x0,y,z0,a,zp,int_v,px,py,pz)
    # print(rta)
    return rta.real

def function_imag_pole_aprox(y_nano,energy0):
    omegac0 = energy0*1e-3/aux 
    y = y_nano*1e-3 
    
    rta = Ex_many_dipoles_pole_aprox(omegac0,epsi1,epsi2,hbmu,hbgama,x0,y,z0,a,zp,int_v,px,py,pz)
    return rta.imag

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

if plot_vs_y == 1: 
    
    listy_re_ana = []    
    listy_im_ana = []

    listy_re_num = []
    listy_im_num = []
 
    listy_re_pole_approx = []   
    listy_im_pole_approx = []   

    for value in listx: 

        y_re_ana = function_real_ana(value,E0)
        y_im_ana = function_imag_ana(value,E0)        
       
        y_re_num = function_real_num(value,E0)
        y_im_num = function_imag_num(value,E0)
        
        y_re_pole_approx = function_real_pole_aprox(value,E0)
        y_im_pole_approx = function_imag_pole_aprox(value,E0)


        listy_re_ana.append(y_re_ana)
        listy_im_ana.append(y_im_ana)

       
        listy_re_num.append(y_re_num)
        listy_im_num.append(y_im_num)

        listy_re_pole_approx.append(y_re_pole_approx)
        listy_im_pole_approx.append(y_im_pole_approx)

     
    graph(title,labelx,'Re($\phi$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    plt.plot(listx,listy_re_ana,'.',ms = ms,color = 'purple',label = 'analytical')
#    plt.plot(listx,listy_re_anaG1_v2,'.',ms = ms,color = 'pink',label = 'analytical v2')
    plt.plot(listx,listy_re_num,'.',ms = 4,color = 'lightseagreen',label = 'numerical')
    plt.plot(listx,listy_re_pole_approx,'.-',ms = 2,color = 'darkred',label = 'pole aprox')
    plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'Re_phi' + label1 + '.png', format='png')   
   
    graph(title,labelx,'Im($\phi$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    plt.plot(listx,listy_im_ana,'.',ms = ms,color = 'purple',label = 'analytical')
#    plt.plot(listx,listy_im_anaG1_v2,'.',ms = ms,color = 'pink',label = 'analytical v2')
    plt.plot(listx,listy_im_num,'.',ms = 4,color = 'lightseagreen',label = 'numerical')
    plt.plot(listx,listy_im_pole_approx,'.-',ms = 2,color = 'darkred',label = 'pole aprox')
    plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'Im_phi' + label1 + '.png', format='png')  

if plot_vs_E == 1: 
    
    listy_re_ana = []    
    listy_im_ana = []

    listy_re_num = []
    listy_im_num = []
 
    listy_re_pole_approx = []   
    listy_im_pole_approx = []   


    for value in listx: 

        y_re_ana = function_real_ana(y0,value)
        y_im_ana = function_imag_ana(y0,value)        
        

        y_re_num = function_real_num(y0,value)
        y_im_num = function_imag_num(y0,value)
        
        y_re_pole_approx = function_real_pole_aprox(y0,value)
        y_im_pole_approx = function_imag_pole_aprox(y0,value)
        

        listy_re_ana.append(y_re_ana)
        listy_im_ana.append(y_im_ana)

       
        listy_re_num.append(y_re_num)
        listy_im_num.append(y_im_num)

        listy_re_pole_approx.append(y_re_pole_approx)
        listy_im_pole_approx.append(y_im_pole_approx)

#        
    graph(title,labelx,'Re($\phi$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    plt.plot(listx,listy_re_ana,'.',ms = ms,color = 'purple',label = 'analytical')
#    plt.plot(listx,listy_re_anaG1_v2,'.',ms = ms,color = 'pink',label = 'analytical v2')
    plt.plot(listx,listy_re_num,'.',ms = 4,color = 'lightseagreen',label = 'numerical')
    plt.plot(listx,listy_re_pole_approx,'.-',ms = 2,color = 'darkred',label = 'pole aprox')
    plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'Re_phi' + label1 + '.png', format='png')   
   
    graph(title,labelx,'Im($\phi$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    plt.plot(listx,listy_im_ana,'.',ms = ms,color = 'purple',label = 'analytical')
#    plt.plot(listx,listy_im_anaG1_v2,'.',ms = ms,color = 'pink',label = 'analytical v2')
    plt.plot(listx,listy_im_num,'.',ms = 4,color = 'lightseagreen',label = 'numerical')
    plt.plot(listx,listy_im_pole_approx,'.-',ms = 2,color = 'darkred',label = 'pole aprox')
    plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'Im_phi' + label1 + '.png', format='png')  

#%%
