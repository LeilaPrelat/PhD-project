
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

plot_vs_x = 0
plot_vs_y = 1

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/potential_field/potential_and_electric_field/potential_and_electric_field_simple_version','')
path_save = path_basic + '/' + 'potential_vs_distance'
if not os.path.exists(path_save):
    print('Creating folder to save graphs')
    os.mkdir(path_save)

err = 'fieldE_direct_numerical.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from potential import electric_potential_ana, electric_potential_num
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


int_v = 10

zp = 0.05
b = -0.01

z0 = zp

xD = 0
yD = 0
zD = 0

px = 1
py = 1
pz = 0

omega0THz = 45
omega0 = omega0THz*1e12 
energy0_pol = omega0*hb
#R = 10 # 10 nm en unidades de micrometros
kappa_factor_omega0 = 0.1
kappa_r_factor= 0.5
 
title1 = r'v = c/%i, $\kappa$ = %.2f$\omega_0$, $\kappa_r$ = %.2f$\kappa$, $E_0$=%.4fmeV' %(int_v, kappa_factor_omega0, kappa_r_factor, energy0_pol*1e3)     
title2 = r'$\hbar\mu$ = %.2feV, $\hbar\gamma$ = %.4feV' %(hbmu,hbgama) 
title4 = r'$z_p$=%inm, px=%i, py=%i, pz=%i, xD=%inm, yD=%inm' %(zp*1e3,px,py,pz,xD*1e3,yD*1e3)
title5 = r'b = %inm, zD=%inm, z=%inm' %(b*1e3,zD*1e3,z0*1e3)
labelp = r'_px%i_py%i_pz%i_omega0THz%i' %(px,py,pz,omega0THz)

N = 400
if plot_vs_x == 1: 
    energy0 = 75*1e-3 # meV 
    y0 = 0.01
    omegac0 = energy0/(hb*c)
#    lambbda = 2*np.pi/omegac
    labelp = r'_px%i_py%i_pz%i_E%imeV' %(px,py,pz,energy0*1e3)
    labelx = 'x [nm]'  
    title5 = title5 + ', E = %.2f meV, y=%inm' %(energy0*1e3,y0*1e3)
    label1 = 'vs_x' + labelp
else:
    energy0 = 30*1e-3 # meV 
    x0 = 0
    
    omegac0 = energy0/(hb*c)
#    lambbda = 2*np.pi/omegac
    labelp = r'_px%i_py%i_pz%i_E%imeV' %(px,py,pz,energy0*1e3)
    labelx = 'y [nm]'  
    title5 = title5 + ', E = %.2f meV, x=%inm' %(energy0*1e3,x0*1e3)
    label1 = 'vs_y' + labelp

listx = np.linspace(10,500,N)  
title = title1 + '\n'  + title2  + '\n'  + title4 + '\n' + title5


def function_real_num(x0_nano,y0_nano):
    y0 = y0_nano*1e-3 
    x0 = x0_nano*1e-3 

    rta = electric_potential_num(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,px,py,pz,b,zp,x0,y0,z0,xD,yD,zD,omega0,kappa_factor_omega0,kappa_r_factor)
    # print(rta)
    return rta.real

def function_imag_num(x0_nano,y0_nano):
    y0 = y0_nano*1e-3 
    x0 = x0_nano*1e-3 
    rta = electric_potential_num(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,px,py,pz,b,zp,x0,y0,z0,xD,yD,zD,omega0,kappa_factor_omega0,kappa_r_factor)
    return rta.imag

def function_abs_num(x0_nano,y0_nano):
    y0 = y0_nano*1e-3 
    x0 = x0_nano*1e-3 
    
    rta = electric_potential_num(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,px,py,pz,b,zp,x0,y0,z0,xD,yD,zD,omega0,kappa_factor_omega0,kappa_r_factor)
    return np.abs(rta)

def function_real_ana(x0_nano,y0_nano):
    y0 = y0_nano*1e-3 
    x0 = x0_nano*1e-3 
    rta = electric_potential_ana(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,px,py,pz,b,zp,x0,y0,z0,xD,yD,zD,omega0,kappa_factor_omega0,kappa_r_factor)
    
    return rta.real

def function_imag_ana(x0_nano,y0_nano):
    y0 = y0_nano*1e-3 
    x0 = x0_nano*1e-3 
    rta =electric_potential_ana(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,px,py,pz,b,zp,x0,y0,z0,xD,yD,zD,omega0,kappa_factor_omega0,kappa_r_factor)
    return rta.imag


def function_abs_ana(x0_nano,y0_nano):
    y0 = y0_nano*1e-3 
    x0 = x0_nano*1e-3 
    rta =electric_potential_ana(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,px,py,pz,b,zp,x0,y0,z0,xD,yD,zD,omega0,kappa_factor_omega0,kappa_r_factor)
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

if plot_vs_x == 1: 
    listy_re_ana = []
    listy_im_ana = []
    
    listy_re_num = []
    listy_im_num = []
    
    for value in listx: 

        y_re_ana = function_real_ana(value,y0)
        y_im_ana = function_imag_ana(value,y0)        

        y_re_num = function_real_num(value,y0)
        y_im_num = function_imag_num(value,y0)
       
        listy_re_ana.append(y_re_ana)
        listy_im_ana.append(y_im_ana)
        
        listy_re_num.append(y_re_num)
        listy_im_num.append(y_im_num)   

 
    graph(title,labelx,'Re($\phi$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    plt.plot(listx,listy_re_ana,'.-',ms = ms,color = 'purple',label = 'analytical')
    plt.plot(listx,listy_re_num,'.-',ms = ms,color = 'lightseagreen',label = 'numerical')
    plt.legend(loc = 'best',markerscale=1.5,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
#    plt.yscale('log')
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'Re_phi' + label1 + '.png', format='png')   
   
    graph(title,labelx,'Im($\phi$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    plt.plot(listx,listy_im_ana,'.',ms = ms,color = 'purple',label = 'analytical')
#    plt.plot(listx,listy_im_num,'.-',ms = 3,color = 'lightseagreen',label = 'sin aprox')
    plt.legend(loc = 'best',markerscale=1.5,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
#    plt.yscale('log')
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'Im_phi_ana' + label1 + '.png', format='png')   

    graph(title,labelx,'Im($\phi$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    plt.plot(listx,listy_im_ana,'.',ms = ms,color = 'purple',label = 'analytical')
    plt.plot(listx,listy_im_num,'.-',ms = ms,color = 'lightseagreen',label = 'numerical')
    plt.legend(loc = 'best',markerscale=1.5,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
#    plt.yscale('log')
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'Im_phi_num' + label1 + '.png', format='png')       

    
if plot_vs_y == 1: 
    
    listy_re_ana = []
    listy_im_ana = []
    listy_abs_ana = []
    
    
    listy_re_num = []
    listy_im_num = []
    listy_abs_num = []    
    
    
    for value in listx: 

        y_re_ana = function_real_ana(x0,value)
        y_im_ana = function_imag_ana(x0,value)        

        y_re_num = function_real_num(x0,value)
        y_im_num = function_imag_num(x0,value)
        
        
        listy_re_ana.append(y_re_ana)
        listy_im_ana.append(y_im_ana)
        listy_abs_ana.append(np.abs(y_re_ana + 1j*y_im_ana))
        
        listy_re_num.append(y_re_num)
        listy_im_num.append(y_im_num)   
        listy_abs_num.append(np.abs(y_re_ana + 1j*y_im_ana))

 
    graph(title,labelx,'Re($\phi$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    plt.plot(listx,listy_re_ana,'.-',ms = ms,color = 'purple',label = 'analytical')
    plt.plot(listx,listy_re_num,'.-',ms = ms,color = 'lightseagreen',label = 'numerical')
    plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
#    plt.yscale('log')
    plt.tight_layout()
    os.chdir(path_save)
#    plt.xlim([-10,600])
    plt.savefig( 'Re_phi_sin_zoom' + label1 + '.png', format='png')   

    graph(title,labelx,'Im($\phi$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    plt.plot(listx,listy_im_ana,'.',ms = ms,color = 'purple',label = 'analytical')
    plt.plot(listx,listy_im_num,'.-',ms = ms,color = 'lightseagreen',label = 'numerical')
    plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
#    plt.yscale('log')
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'Im_phi_sin_zoom' + label1 + '.png', format='png')    



    graph(title,labelx,'$|\phi|$',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    plt.plot(listx,listy_abs_ana,'.',ms = ms,color = 'purple',label = 'analytical')
    plt.plot(listx,listy_abs_num,'.-',ms = ms,color = 'lightseagreen',label = 'numerical')
    plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
#    plt.yscale('log')
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'Abs_phi_sin_zoom' + label1 + '.png', format='png')    
  
#%%
