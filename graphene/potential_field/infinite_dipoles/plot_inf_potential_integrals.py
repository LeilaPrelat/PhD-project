
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

plot_vs_y = 1
plot_vs_E = 0

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/potential_field/infinite_dipoles','')
path_save = path_basic + '/' + 'integrals_inf_potential'
if not os.path.exists(path_save):
    print('Creating folder to save graphs')
    os.mkdir(path_save)

err = 'many_potential_integrals.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from inf_potential_integrals_final import F1_ana, F1_num, F1_pole_aprox, F2_ana, F2_num, F2_pole_aprox, F3_ana, F3_num, F3_pole_aprox 
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

epsi1,epsi2 = 1,1
hbmu,hbgama = 0.3,0.0001
zp = 0.05 ## micrones

int_v = 10
a = 50*1e-3 ## micrones
n = 0
z = 0 ## micrones

title2 = r'$\hbar\mu$ = %.2f eV, v = c/%i, $z_{\rm p}$ = %i nm' %(hbmu,int_v,zp*1e3) 
title4 = r'n = %i, a = %i nm, z = %i nm' %(n,a*1e3,z*1e3)

N = 250
if plot_vs_y == 1: 
    E0 = 30 # meV 
    labelx = 'y (nm)'  
    title4 = title4 + ', ' + 'E = %.2f meV' %(E0)
    label1 = 'vs_y_E0%i' %(E0) 
    listx = np.linspace(10,10000,N)
else:
    
    y0 = 5*1e3
#    y0 = 0
    # z0 = 0.06*1e3
    labelx = r'$\hbar\omega$ (meV)'   
    title4 = title4 + ', ' + 'y = %i nm' %(y0)
    label1 = 'vs_E' 
    listx = np.linspace(20,80,N)
  
title = title2 + '\n' + title4 

#%%

def function_G1_real_num(y_nano,energy0):
    omegac0 = energy0*1e-3/aux 
    y = y_nano*1e-3 

    rta = F1_num(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,zp,a,n,y,z)
    # print(rta)
    return rta.real

def function_G1_imag_num(y_nano,energy0):
    omegac0 = energy0*1e-3/aux 
    y = y_nano*1e-3 
    
    rta = F1_num(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,zp,a,n,y,z)
    return rta.imag

def function_G2_real_num(y_nano,energy0):
    omegac0 = energy0*1e-3/aux 
    y = y_nano*1e-3 

    rta = F2_num(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,zp,a,n,y,z)
    # print(rta)
    return rta.real

def function_G2_imag_num(y_nano,energy0):
    omegac0 = energy0*1e-3/aux 
    y = y_nano*1e-3 
    
    rta = F2_num(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,zp,a,n,y,z)
    return rta.imag


def function_G3_real_num(y_nano,energy0):
    omegac0 = energy0*1e-3/aux 
    y = y_nano*1e-3 

    rta = F3_num(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,zp,a,n,y,z)
    # print(rta)
    return rta.real

def function_G3_imag_num(y_nano,energy0):
    omegac0 = energy0*1e-3/aux 
    y = y_nano*1e-3 
    
    rta = F3_num(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,zp,a,n,y,z)
    return rta.imag


#%%
    
def function_G1_real_ana(y_nano,energy0):
    omegac0 = energy0*1e-3/aux 
    y = y_nano*1e-3 

    rta = F1_ana(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,zp,a,n,y,z)
    # print(rta)
    return rta.real

def function_G1_imag_ana(y_nano,energy0):
    omegac0 = energy0*1e-3/aux 
    y = y_nano*1e-3 
    
    rta = F1_ana(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,zp,a,n,y,z)
    return rta.imag


def function_G2_real_ana(y_nano,energy0):
    omegac0 = energy0*1e-3/aux 
    y = y_nano*1e-3 

    rta = F2_ana(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,zp,a,n,y,z)
    # print(rta)
    return rta.real

def function_G2_imag_ana(y_nano,energy0):
    omegac0 = energy0*1e-3/aux 
    y = y_nano*1e-3 
    
    rta = F2_ana(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,zp,a,n,y,z)
    return rta.imag

def function_G3_real_ana(y_nano,energy0):
    omegac0 = energy0*1e-3/aux 
    y = y_nano*1e-3 

    rta = F3_ana(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,zp,a,n,y,z)
    # print(rta)
    return rta.real

def function_G3_imag_ana(y_nano,energy0):
    omegac0 = energy0*1e-3/aux 
    y = y_nano*1e-3 
    
    rta = F3_ana(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,zp,a,n,y,z)
    return rta.imag

#%%
    
def function_G1_real_pole_aprox(y_nano,energy0):
    omegac0 = energy0*1e-3/aux 
    y = y_nano*1e-3 

    rta = F1_pole_aprox(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,zp,a,n,y,z)
    # print(rta)
    return rta.real

def function_G1_imag_pole_aprox(y_nano,energy0):
    omegac0 = energy0*1e-3/aux 
    y = y_nano*1e-3 
    
    rta = F1_pole_aprox(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,zp,a,n,y,z)
    return rta.imag


def function_G2_real_pole_aprox(y_nano,energy0):
    omegac0 = energy0*1e-3/aux 
    y = y_nano*1e-3 

    rta = F2_pole_aprox(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,zp,a,n,y,z)
    # print(rta)
    return rta.real

def function_G2_imag_pole_aprox(y_nano,energy0):
    omegac0 = energy0*1e-3/aux 
    y = y_nano*1e-3 
    
    rta = F2_pole_aprox(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,zp,a,n,y,z)
    return rta.imag


def function_G3_real_pole_aprox(y_nano,energy0):
    omegac0 = energy0*1e-3/aux 
    y = y_nano*1e-3 

    rta = F3_pole_aprox(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,zp,a,n,y,z)
    # print(rta)
    return rta.real

def function_G3_imag_pole_aprox(y_nano,energy0):
    omegac0 = energy0*1e-3/aux 
    y = y_nano*1e-3 
    
    rta = F3_pole_aprox(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,zp,a,n,y,z)
    return rta.imag

#%%
    
tamfig = [3.5, 3]
tamletra = 9
tamtitle  = 9
tamnum = 7
tamlegend = 7
labelpady = 2
labelpadx = 3
pad = 2.5
mk = 2
ms = 2
hp = 0.5
length_marker = 1.5
dpi = 500


def graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad):
    plt.figure(figsize=tamfig)
    plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
    plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
    plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in", pad = pad)
#    plt.tick_params(labelsize = tamnum, length = 4 , width=1, direction="in", pad = pad) ### para cleo europe
    plt.title(title,fontsize=int(tamtitle*0.9))

    return  
 

#%%

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

if plot_vs_y == 1: 
    
    for value in listx: 

        y_re_anaG1 = function_G1_real_ana(value,E0)
        y_im_anaG1 = function_G1_imag_ana(value,E0)        
       
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


if plot_vs_E == 1: 

    for value in listx: 

        y_re_anaG1 = function_G1_real_ana(y0,value)
        y_im_anaG1 = function_G1_imag_ana(y0,value)        
        

        y_re_numG1 = function_G1_real_num(y0,value)
        y_im_numG1 = function_G1_imag_num(y0,value)
        
        y_re_pole_approxG1 = function_G1_real_pole_aprox(y0,value)
        y_im_pole_approxG1 = function_G1_imag_pole_aprox(y0,value)


#
        y_re_anaG2 = function_G2_real_ana(y0,value)
        y_im_anaG2 = function_G2_imag_ana(y0,value)        

#
        y_re_anaG3 = function_G3_real_ana(y0,value)
        y_im_anaG3 = function_G3_imag_ana(y0,value)  

   
        y_re_numG2 = function_G2_real_num(y0,value)
        y_im_numG2 = function_G2_imag_num(y0,value)
        
        y_re_pole_approxG2 = function_G2_real_pole_aprox(y0,value)
        y_im_pole_approxG2 = function_G2_imag_pole_aprox(y0,value)


        y_re_anaG3 = function_G3_real_ana(y0,value)
        y_im_anaG3 = function_G3_imag_ana(y0,value)        
    
        y_re_numG3 = function_G3_real_num(y0,value)
        y_im_numG3 = function_G3_imag_num(y0,value)
        
        y_re_pole_approxG3 = function_G3_real_pole_aprox(y0,value)
        y_im_pole_approxG3 = function_G3_imag_pole_aprox(y0,value)


        listy_re_anaG1.append(y_re_anaG1)
        listy_im_anaG1.append(y_im_anaG1)
                        
        
    
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
#        

        listy_re_anaG3.append(y_re_anaG3)
        listy_im_anaG3.append(y_im_anaG3)
        
        
        listy_re_numG3.append(y_re_numG3)
        listy_im_numG3.append(y_im_numG3)

        listy_re_pole_approx_G3.append(y_re_pole_approxG3)
        listy_im_pole_approx_G3.append(y_im_pole_approxG3)
#        
    #%%  
    
    
graph(title,labelx,'Re{$F_1$}',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
 #   plt.plot(listx,listy_re_anaG1,'.',ms = ms,color = 'purple',label = 'analytical')
plt.plot(listx,listy_re_anaG1,'.',ms = ms,color = 'purple',label = 'PP analytical')
#    plt.plot(listx,listy_re_anaG1_v2,'.',ms = ms,color = 'pink',label = 'analytical v2')
plt.plot(listx,listy_re_numG1,'.',ms = ms,color = 'lightseagreen',label = 'full numerical')
plt.plot(listx,listy_re_pole_approx_G1,'.-',ms = ms,color = 'darkred',label = 'PP numerical')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=0.1,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
os.chdir(path_save)
plt.savefig( 'Re_G1' + label1 + '.png', format='png')   
   
graph(title,labelx,'Im{$F_1$}',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
 #   plt.plot(listx,listy_im_anaG1,'.',ms = ms,color = 'purple',label = 'analytical')
plt.plot(listx,listy_im_anaG1,'.',ms = ms,color = 'purple',label = 'PP analytical')
#    plt.plot(listx,listy_im_anaG1_v2,'.',ms = ms,color = 'pink',label = 'analytical v2')
plt.plot(listx,listy_im_numG1,'.',ms = ms,color = 'lightseagreen',label = 'full numerical')
plt.plot(listx,listy_im_pole_approx_G1,'.-',ms = ms,color = 'darkred',label = 'PP numerical')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=0.1,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
os.chdir(path_save)
plt.savefig( 'Im_G1' + label1 + '.png', format='png')  



    
graph(title,labelx,'Re{$F_2$} ($\mu$m$^{-1}$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#    plt.plot(listx,listy_re_anaG2,'.',ms = ms,color = 'purple',label = 'analytical')
plt.plot(listx,listy_re_anaG2,'.',ms = ms,color = 'purple',label = 'PP analytical')
plt.plot(listx,listy_re_numG2,'.',ms = ms,color = 'lightseagreen',label = 'full numerical')
plt.plot(listx,listy_re_pole_approx_G2,'.-',ms = ms,color = 'darkred',label = 'PP numerical')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=0.1,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
os.chdir(path_save)
plt.savefig( 'Re_G2' + label1 + '.png', format='png')   
   
graph(title,labelx,'Im{$F_2$} ($\mu$m$^{-1}$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#    plt.plot(listx,listy_im_anaG2,'.',ms = ms,color = 'purple',label = 'analytical')
plt.plot(listx,listy_im_anaG2,'.',ms = ms,color = 'purple',label = 'PP analytical')
plt.plot(listx,listy_im_numG2,'.',ms = ms,color = 'lightseagreen',label = 'full numerical')
plt.plot(listx,listy_im_pole_approx_G2,'.-',ms = ms,color = 'darkred',label = 'PP numerical')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=0.1,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
os.chdir(path_save)
plt.savefig( 'Im_G2' + label1 + '.png', format='png')  



graph(title,labelx,'Re{$F_3$} ($\mu$m$^{-1}$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#    plt.plot(listx,listy_re_anaG3,'.',ms = ms,color = 'purple',label = 'analytical')
plt.plot(listx,listy_re_anaG3,'.-',ms = ms,color = 'purple',label = 'PP analytical')
plt.plot(listx,listy_re_numG3,'.',ms = ms,color = 'lightseagreen',label = 'full numerical')
plt.plot(listx,listy_re_pole_approx_G3,'.-',ms = ms,color = 'darkred',label =  'PP numerical')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=0.1,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
os.chdir(path_save)
plt.savefig( 'Re_G3' + label1 + '.png', format='png')   
   
graph(title,labelx,'Im{$F_3$} ($\mu$m$^{-1}$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#    plt.plot(listx,listy_im_anaG3,'.',ms = ms,color = 'purple',label = 'analytical')
plt.plot(listx,listy_im_anaG3,'.-',ms = ms,color = 'purple',label = 'PP analytical')
plt.plot(listx,listy_im_numG3,'.',ms = ms,color = 'lightseagreen',label = 'full numerical')
plt.plot(listx,listy_im_pole_approx_G3,'.-',ms = ms,color = 'darkred',label =  'PP numerical')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=0.1,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
os.chdir(path_save)
plt.savefig( 'Im_G3' + label1 + '.png', format='png')  
  



#%%
