
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
path_save = path_basic + '/' + 'integrals'
if not os.path.exists(path_save):
    print('Creating folder to save graphs')
    os.mkdir(path_save)

err = 'many_potential_integrals.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from many_potential_integrals import G1_ana, G1_ana_v3, G1_num, G1_pole_aprox, G2_ana, G2_ana_v3, G2_num, G2_pole_aprox, G1_ana_v2, G3_ana, G3_ana_v3, G3_num, G3_num_v3, G3_pole_aprox
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

int_v = 10

title2 = r'$\hbar\mu$ = %.2feV, $\hbar\gamma$ = %.4feV, v = c/%i' %(hbmu,hbgama,int_v) 
title4 = r'$z_p$=%inm' %(zp*1e3)

N = 250
if plot_vs_y == 1: 
    E0 = 30 # meV 
    labelx = 'y [nm]'  
    title4 = title4 + ', ' + 'E = %.2f meV' %(E0)
    label1 = 'vs_y_E0%i' %(E0) 
    listx = np.linspace(10,2000,N)
else:
    
    y0 = 5*1e3
#    y0 = 0
    # z0 = 0.06*1e3
    labelx = 'E [meV]'   
    title4 = title4 + ', ' + 'y = %inm' %(y0)
    label1 = 'vs_E' 
    listx = np.linspace(80,20,N)
  
title = title2 + '\n' + title4 

#%%

def function_G1_real_num(y_nano,energy0):
    omegac0 = energy0*1e-3/aux 
    y = y_nano*1e-3 

    rta = G1_num(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,zp,y)
    # print(rta)
    return rta.real

def function_G1_imag_num(y_nano,energy0):
    omegac0 = energy0*1e-3/aux 
    y = y_nano*1e-3 
    
    rta = G1_num(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,zp,y)
    return rta.imag

def function_G2_real_num(y_nano,energy0):
    omegac0 = energy0*1e-3/aux 
    y = y_nano*1e-3 

    rta = G2_num(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,zp,y)
    # print(rta)
    return rta.real

def function_G2_imag_num(y_nano,energy0):
    omegac0 = energy0*1e-3/aux 
    y = y_nano*1e-3 
    
    rta = G2_num(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,zp,y)
    return rta.imag


def function_G3_real_num(y_nano,energy0):
    omegac0 = energy0*1e-3/aux 
    y = y_nano*1e-3 

    rta = G3_num_v3(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,zp,y)
    # print(rta)
    return rta.real

def function_G3_imag_num(y_nano,energy0):
    omegac0 = energy0*1e-3/aux 
    y = y_nano*1e-3 
    
    rta = G3_num_v3(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,zp,y)
    return rta.imag


#%%
    
def function_G1_real_ana(y_nano,energy0):
    omegac0 = energy0*1e-3/aux 
    y = y_nano*1e-3 

    rta = G1_ana(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,zp,y)
    # print(rta)
    return rta.real

def function_G1_imag_ana(y_nano,energy0):
    omegac0 = energy0*1e-3/aux 
    y = y_nano*1e-3 
    
    rta = G1_ana(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,zp,y)
    return rta.imag

    
def function_G1_v2_real_ana(y_nano,energy0):
    omegac0 = energy0*1e-3/aux 
    y = y_nano*1e-3 

    rta = G1_ana_v2(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,zp,y)
    # print(rta)
    return rta.real

def function_G1_v2_imag_ana(y_nano,energy0):
    omegac0 = energy0*1e-3/aux 
    y = y_nano*1e-3 
    
    rta = G1_ana_v2(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,zp,y)
    return rta.imag

def function_G1_v3_real_ana(y_nano,energy0):
    omegac0 = energy0*1e-3/aux 
    y = y_nano*1e-3 

    rta = G1_ana_v3(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,zp,y)
    # print(rta)
    return rta.real

def function_G1_v3_imag_ana(y_nano,energy0):
    omegac0 = energy0*1e-3/aux 
    y = y_nano*1e-3 
    
    rta = G1_ana_v3(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,zp,y)
    return rta.imag

def function_G2_real_ana(y_nano,energy0):
    omegac0 = energy0*1e-3/aux 
    y = y_nano*1e-3 

    rta = G2_ana(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,zp,y)
    # print(rta)
    return rta.real

def function_G2_imag_ana(y_nano,energy0):
    omegac0 = energy0*1e-3/aux 
    y = y_nano*1e-3 
    
    rta = G2_ana(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,zp,y)
    return rta.imag

def function_G2_v3_real_ana(y_nano,energy0):
    omegac0 = energy0*1e-3/aux 
    y = y_nano*1e-3 

    rta = G2_ana_v3(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,zp,y)
    # print(rta)
    return rta.real

def function_G2_v3_imag_ana(y_nano,energy0):
    omegac0 = energy0*1e-3/aux 
    y = y_nano*1e-3 
    
    rta = G2_ana_v3(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,zp,y)
    return rta.imag


def function_G3_real_ana(y_nano,energy0):
    omegac0 = energy0*1e-3/aux 
    y = y_nano*1e-3 

    rta = G3_ana(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,zp,y)
    # print(rta)
    return rta.real

def function_G3_imag_ana(y_nano,energy0):
    omegac0 = energy0*1e-3/aux 
    y = y_nano*1e-3 
    
    rta = G3_ana(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,zp,y)
    return rta.imag


def function_G3_v3_real_ana(y_nano,energy0):
    omegac0 = energy0*1e-3/aux 
    y = y_nano*1e-3 

    rta = G3_ana_v3(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,zp,y)
    # print(rta)
    return rta.real

def function_G3_v3_imag_ana(y_nano,energy0):
    omegac0 = energy0*1e-3/aux 
    y = y_nano*1e-3 
    
    rta = G3_ana_v3(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,zp,y)
    return rta.imag


#%%
    
def function_G1_real_pole_aprox(y_nano,energy0):
    omegac0 = energy0*1e-3/aux 
    y = y_nano*1e-3 

    rta = G1_pole_aprox(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,zp,y)
    # print(rta)
    return rta.real

def function_G1_imag_pole_aprox(y_nano,energy0):
    omegac0 = energy0*1e-3/aux 
    y = y_nano*1e-3 
    
    rta = G1_pole_aprox(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,zp,y)
    return rta.imag


def function_G2_real_pole_aprox(y_nano,energy0):
    omegac0 = energy0*1e-3/aux 
    y = y_nano*1e-3 

    rta = G2_pole_aprox(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,zp,y)
    # print(rta)
    return rta.real

def function_G2_imag_pole_aprox(y_nano,energy0):
    omegac0 = energy0*1e-3/aux 
    y = y_nano*1e-3 
    
    rta = G2_pole_aprox(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,zp,y)
    return rta.imag


def function_G3_real_pole_aprox(y_nano,energy0):
    omegac0 = energy0*1e-3/aux 
    y = y_nano*1e-3 

    rta = G3_pole_aprox(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,zp,y)
    # print(rta)
    return rta.real

def function_G3_imag_pole_aprox(y_nano,energy0):
    omegac0 = energy0*1e-3/aux 
    y = y_nano*1e-3 
    
    rta = G3_pole_aprox(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,zp,y)
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
    
    listy_re_anaG1 = []    
    listy_im_anaG1 = []


    listy_re_anaG1_v2 = []    
    listy_im_anaG1_v2 = []
    
    listy_re_numG1 = []
    listy_im_numG1 = []
 
    listy_re_pole_approx_G1 = []   
    listy_im_pole_approx_G1 = []   
    

    listy_re_anaG2 = []
    listy_im_anaG2 = []
    
    listy_re_numG2 = []
    listy_im_numG2 = []
 
    listy_re_pole_approx_G2 = []   
    listy_im_pole_approx_G2 = []   

    listy_re_anaG3 = []
    listy_im_anaG3 = []
    
    listy_re_numG3 = []
    listy_im_numG3 = []
 
    listy_re_pole_approx_G3 = []   
    listy_im_pole_approx_G3 = []   

    for value in listx: 

        y_re_anaG1 = function_G1_real_ana(value,E0)
        y_im_anaG1 = function_G1_imag_ana(value,E0)        


        y_re_anaG1_v2 = function_G1_v2_real_ana(value,E0)
        y_im_anaG1_v2 = function_G1_v2_imag_ana(value,E0)    
       
        y_re_numG1 = function_G1_real_num(value,E0)
        y_im_numG1 = function_G1_imag_num(value,E0)
        
        y_re_pole_approxG1 = function_G1_real_pole_aprox(value,E0)
        y_im_pole_approxG1 = function_G1_imag_pole_aprox(value,E0)

        y_re_anaG2 = function_G2_real_ana(value,E0)
        y_im_anaG2 = function_G2_imag_ana(value,E0)        
    
        y_re_numG2 = function_G2_real_num(value,E0)
        y_im_numG2 = function_G2_imag_num(value,E0)
        
        y_re_pole_approxG2 = function_G2_real_pole_aprox(value,E0)
        y_im_pole_approxG2 = function_G2_imag_pole_aprox(value,E0)


        y_re_anaG3 = function_G3_real_ana(value,E0)
        y_im_anaG3 = function_G3_imag_ana(value,E0)        
    
        y_re_numG3 = function_G3_real_num(value,E0)
        y_im_numG3 = function_G3_imag_num(value,E0)
        
        y_re_pole_approxG3 = function_G3_real_pole_aprox(value,E0)
        y_im_pole_approxG3 = function_G3_imag_pole_aprox(value,E0)

        listy_re_anaG1.append(y_re_anaG1)
        listy_im_anaG1.append(y_im_anaG1)


        listy_re_anaG1_v2.append(y_re_anaG1_v2)
        listy_im_anaG1_v2.append(y_im_anaG1_v2)
       
        listy_re_numG1.append(y_re_numG1)
        listy_im_numG1.append(y_im_numG1)

        listy_re_pole_approx_G1.append(y_re_pole_approxG1)
        listy_im_pole_approx_G1.append(y_im_pole_approxG1)
        
        
        listy_re_anaG2.append(y_re_anaG2)
        listy_im_anaG2.append(y_im_anaG2)
        
        listy_re_numG2.append(y_re_numG2)
        listy_im_numG2.append(y_im_numG2)

        listy_re_pole_approx_G2.append(y_re_pole_approxG2)
        listy_im_pole_approx_G2.append(y_im_pole_approxG2)
        

        listy_re_anaG3.append(y_re_anaG3)
        listy_im_anaG3.append(y_im_anaG3)
        
        listy_re_numG3.append(y_re_numG3)
        listy_im_numG3.append(y_im_numG3)

        listy_re_pole_approx_G3.append(y_re_pole_approxG3)
        listy_im_pole_approx_G3.append(y_im_pole_approxG3)

     
    graph(title,labelx,'Re($G_1$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    plt.plot(listx,listy_re_anaG1,'.',ms = ms,color = 'purple',label = 'analytical')
#    plt.plot(listx,listy_re_anaG1_v2,'.',ms = ms,color = 'pink',label = 'analytical v2')
    plt.plot(listx,listy_re_numG1,'.',ms = 4,color = 'lightseagreen',label = 'numerical')
    plt.plot(listx,listy_re_pole_approx_G1,'.-',ms = 2,color = 'darkred',label = 'pole aprox')
    plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'Re_G1' + label1 + '.png', format='png')   
   
    graph(title,labelx,'Im($G_1$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    plt.plot(listx,listy_im_anaG1,'.',ms = ms,color = 'purple',label = 'analytical')
#    plt.plot(listx,listy_im_anaG1_v2,'.',ms = ms,color = 'pink',label = 'analytical v2')
    plt.plot(listx,listy_im_numG1,'.',ms = 4,color = 'lightseagreen',label = 'numerical')
    plt.plot(listx,listy_im_pole_approx_G1,'.-',ms = 2,color = 'darkred',label = 'pole aprox')
    plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'Im_G1' + label1 + '.png', format='png')  



#        
    graph(title,labelx,'Re($G_2$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    plt.plot(listx,listy_re_anaG2,'.',ms = ms,color = 'purple',label = 'analytical')
    plt.plot(listx,listy_re_numG2,'.',ms = 4,color = 'lightseagreen',label = 'numerical')
    plt.plot(listx,listy_re_pole_approx_G2,'.-',ms = 2,color = 'darkred',label = 'pole aprox')
    plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'Re_G2' + label1 + '.png', format='png')   
   
    graph(title,labelx,'Im($G_2$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    plt.plot(listx,listy_im_anaG2,'.',ms = ms,color = 'purple',label = 'analytical')
    plt.plot(listx,listy_im_numG2,'.',ms = 4,color = 'lightseagreen',label = 'numerical')
    plt.plot(listx,listy_im_pole_approx_G2,'.-',ms = 2,color = 'darkred',label = 'pole aprox')
    plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'Im_G2' + label1 + '.png', format='png')  



    graph(title,labelx,'Re($G_3$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    plt.plot(listx,listy_re_anaG3,'.',ms = ms,color = 'purple',label = 'analytical')
    plt.plot(listx,listy_re_numG3,'.',ms = 4,color = 'lightseagreen',label = 'numerical')
    plt.plot(listx,listy_re_pole_approx_G3,'.-',ms = 2,color = 'darkred',label = 'pole aprox')
    plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'Re_G3' + label1 + '.png', format='png')   
   
    graph(title,labelx,'Im($G_3$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    plt.plot(listx,listy_im_anaG3,'.',ms = ms,color = 'purple',label = 'analytical')
    plt.plot(listx,listy_im_numG3,'.',ms = 4,color = 'lightseagreen',label = 'numerical')
    plt.plot(listx,listy_im_pole_approx_G3,'.-',ms = 2,color = 'darkred',label = 'pole aprox')
    plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'Im_G3' + label1 + '.png', format='png')  

if plot_vs_E == 1: 
    
    
    listy_re_anaG1 = []
    listy_im_anaG1 = []

    listy_re_anaG1_v2 = []
    listy_im_anaG1_v2 = []


    listy_re_anaG1_v3 = []    
    listy_im_anaG1_v3 = []
 
    listy_re_numG1 = []
    listy_im_numG1 = []
 
    listy_re_pole_approx_G1 = []   
    listy_im_pole_approx_G1 = []   
    

    listy_re_anaG2 = []
    listy_im_anaG2 = []


    listy_re_anaG2_v3 = []
    listy_im_anaG2_v3 = []

   
    listy_re_numG2 = []
    listy_im_numG2 = []
 
    listy_re_pole_approx_G2 = []   
    listy_im_pole_approx_G2 = []   

    listy_re_anaG3 = []
    listy_im_anaG3 = []

    listy_re_anaG3_v3 = []
    listy_im_anaG3_v3 = []

    
    listy_re_numG3 = []
    listy_im_numG3 = []
 
    listy_re_pole_approx_G3 = []   
    listy_im_pole_approx_G3 = []   

    for value in listx: 

        y_re_anaG1 = function_G1_real_ana(y0,value)
        y_im_anaG1 = function_G1_imag_ana(y0,value)        
        
        y_re_anaG1_v2 = function_G1_v2_real_ana(y0,value)
        y_im_anaG1_v2 = function_G1_v2_imag_ana(y0,value)        

        y_re_anaG1_v3 = function_G1_v3_real_ana(y0,value)
        y_im_anaG1_v3 = function_G1_v3_imag_ana(y0,value)        


        y_re_numG1 = function_G1_real_num(y0,value)
        y_im_numG1 = function_G1_imag_num(y0,value)
        
        y_re_pole_approxG1 = function_G1_real_pole_aprox(y0,value)
        y_im_pole_approxG1 = function_G1_imag_pole_aprox(y0,value)


#
        y_re_anaG2 = function_G2_real_ana(y0,value)
        y_im_anaG2 = function_G2_imag_ana(y0,value)        

        y_re_anaG2_v3 = function_G2_v3_real_ana(y0,value)
        y_im_anaG2_v3 = function_G2_v3_imag_ana(y0,value)    
#
        y_re_anaG3 = function_G3_real_ana(y0,value)
        y_im_anaG3 = function_G3_imag_ana(y0,value)  

   
        y_re_numG2 = function_G2_real_num(y0,value)
        y_im_numG2 = function_G2_imag_num(y0,value)
        
        y_re_pole_approxG2 = function_G2_real_pole_aprox(y0,value)
        y_im_pole_approxG2 = function_G2_imag_pole_aprox(y0,value)


        y_re_anaG3 = function_G3_real_ana(y0,value)
        y_im_anaG3 = function_G3_imag_ana(y0,value)        

        y_re_anaG3_v3 = function_G3_v3_real_ana(y0,value)
        y_im_anaG3_v3 = function_G3_v3_imag_ana(y0,value)    

    
        y_re_numG3 = function_G3_real_num(y0,value)
        y_im_numG3 = function_G3_imag_num(y0,value)
        
        y_re_pole_approxG3 = function_G3_real_pole_aprox(y0,value)
        y_im_pole_approxG3 = function_G3_imag_pole_aprox(y0,value)


        listy_re_anaG1.append(y_re_anaG1)
        listy_im_anaG1.append(y_im_anaG1)
        
        listy_re_anaG1_v2.append(y_re_anaG1_v2)
        listy_im_anaG1_v2.append(y_im_anaG1_v2)
        
        
        listy_re_anaG1_v3.append(y_re_anaG1_v3)
        listy_im_anaG1_v3.append(y_im_anaG1_v3)
                
        
    
        listy_re_numG1.append(y_re_numG1)
        listy_im_numG1.append(y_im_numG1)

        listy_re_pole_approx_G1.append(y_re_pole_approxG1)
        listy_im_pole_approx_G1.append(y_im_pole_approxG1)
        
        
        listy_re_anaG2.append(y_re_anaG2)
        listy_im_anaG2.append(y_im_anaG2)


        listy_re_anaG2_v3.append(y_re_anaG2_v3)
        listy_im_anaG2_v3.append(y_im_anaG2_v3)

        listy_re_numG2.append(y_re_numG2)
        listy_im_numG2.append(y_im_numG2)

        listy_re_pole_approx_G2.append(y_re_pole_approxG2)
        listy_im_pole_approx_G2.append(y_im_pole_approxG2)
#        

        listy_re_anaG3.append(y_re_anaG3)
        listy_im_anaG3.append(y_im_anaG3)

        listy_re_anaG3_v3.append(y_re_anaG3_v3)
        listy_im_anaG3_v3.append(y_im_anaG3_v3)
        
        
        listy_re_numG3.append(y_re_numG3)
        listy_im_numG3.append(y_im_numG3)

        listy_re_pole_approx_G3.append(y_re_pole_approxG3)
        listy_im_pole_approx_G3.append(y_im_pole_approxG3)
#        
       
    graph(title,labelx,'Re($G_1$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
 #   plt.plot(listx,listy_re_anaG1,'.',ms = ms,color = 'purple',label = 'analytical')
    plt.plot(listx,listy_re_anaG1_v3,'.',ms = ms,color = 'purple',label = 'PP analytical')
#    plt.plot(listx,listy_re_anaG1_v2,'.',ms = ms,color = 'pink',label = 'analytical v2')
    plt.plot(listx,listy_re_numG1,'.',ms = 4,color = 'lightseagreen',label = 'full numerical')
    plt.plot(listx,listy_re_pole_approx_G1,'.-',ms = 2,color = 'darkred',label = 'PP numerical')
    plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'Re_G1' + label1 + '.png', format='png')   
   
    graph(title,labelx,'Im($G_1$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
 #   plt.plot(listx,listy_im_anaG1,'.',ms = ms,color = 'purple',label = 'analytical')
    plt.plot(listx,listy_im_anaG1_v3,'.',ms = ms,color = 'purple',label = 'PP analytical')
#    plt.plot(listx,listy_im_anaG1_v2,'.',ms = ms,color = 'pink',label = 'analytical v2')
    plt.plot(listx,listy_im_numG1,'.',ms = 4,color = 'lightseagreen',label = 'full numerical')
    plt.plot(listx,listy_im_pole_approx_G1,'.-',ms = 2,color = 'darkred',label = 'PP numerical')
    plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'Im_G1' + label1 + '.png', format='png')  



        
    graph(title,labelx,'Re($G_2$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#    plt.plot(listx,listy_re_anaG2,'.',ms = ms,color = 'purple',label = 'analytical')
    plt.plot(listx,listy_re_anaG2_v3,'.',ms = ms,color = 'purple',label = 'PP analytical')
    plt.plot(listx,listy_re_numG2,'.',ms = 4,color = 'lightseagreen',label = 'full numerical')
    plt.plot(listx,listy_re_pole_approx_G2,'.-',ms = 2,color = 'darkred',label = 'PP numerical')
    plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'Re_G2' + label1 + '.png', format='png')   
   
    graph(title,labelx,'Im($G_2$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#    plt.plot(listx,listy_im_anaG2,'.',ms = ms,color = 'purple',label = 'analytical')
    plt.plot(listx,listy_im_anaG2_v3,'.',ms = ms,color = 'purple',label = 'PP analytical')
    plt.plot(listx,listy_im_numG2,'.',ms = 4,color = 'lightseagreen',label = 'full numerical')
    plt.plot(listx,listy_im_pole_approx_G2,'.-',ms = 2,color = 'darkred',label = 'PP numerical')
    plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'Im_G2' + label1 + '.png', format='png')  



    graph(title,labelx,'Re($G_3$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#    plt.plot(listx,listy_re_anaG3,'.',ms = ms,color = 'purple',label = 'analytical')
    plt.plot(listx,listy_re_anaG3_v3,'.-',ms = 4,color = 'purple',label = 'PP analytical')
    plt.plot(listx,listy_re_numG3,'.',ms = 4,color = 'lightseagreen',label = 'full numerical')
    plt.plot(listx,listy_re_pole_approx_G3,'.-',ms = 2,color = 'darkred',label =  'PP numerical')
    plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'Re_G3' + label1 + '.png', format='png')   
   
    graph(title,labelx,'Im($G_3$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#    plt.plot(listx,listy_im_anaG3,'.',ms = ms,color = 'purple',label = 'analytical')
    plt.plot(listx,listy_im_anaG3_v3,'.-',ms = 4,color = 'purple',label = 'PP analytical')
    plt.plot(listx,listy_im_numG3,'.',ms = 4,color = 'lightseagreen',label = 'full numerical')
    plt.plot(listx,listy_im_pole_approx_G3,'.-',ms = 2,color = 'darkred',label =  'PP numerical')
    plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'Im_G3' + label1 + '.png', format='png')  
  



#%%
