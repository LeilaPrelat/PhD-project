
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

plot_vs_b = 0
plot_vs_E = 1

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/potential_field/potential/potential_and_electric_field_with_dipole_moment_formula','')
path_save = path_basic + '/' + 'dipole_moment_sin_alpha_eff'
if not os.path.exists(path_save):
    print('Creating folder to save graphs')
    os.mkdir(path_save)

err = 'fieldE_direct_numerical.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from dipole_moment_sin_alpha_eff import dipole_moment_ana, dipole_moment_num
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

int_v = 10
#v = c/int_v
#omega = 0.7*1e12
cota = 75 #nanometros
epsi1,epsi2 = 1,1
hbmu,hbgama = 0.3,0.0001
zp = 0.05
b = -0.01

y0,x0 = 0,0 

xD = 0.005
yD = 0.005
zD = 0

px = 1
py = 1
pz = 1

# omega0THz = 21
# omega0 = omega0THz*1e12 
# energy0_pol = omega0*hb
# #R = 10 # 10 nm en unidades de micrometros
# kappa_factor_omega0 = 0.1
# kappa_r_factor= 0.5
 
title1 = r'v = c/%i' %(int_v)     
title2 = r'$\hbar\mu$ = %.2feV, $\hbar\gamma$ = %.4feV' %(hbmu,hbgama) 
title4 = r'$z_p$=%inm, px=%i, py=%i, pz=%i, xD=%inm, yD=%inm' %(zp*1e3,px,py,pz,xD*1e3,yD*1e3)
title5 = r'b = %inm, zD=%inm, x=%inm, y=%inm' %(b*1e3,zD*1e3,x0*1e3,y0*1e3)
labelp = r'_px%i_py%i_pz%i' %(px,py,pz)

N = 100
if plot_vs_b == 1: 
    E0 = 17 # meV 
    labelx = 'b [nm]'  
    title5 = title5 + ', ' + 'E = %.2f meV' %(E0)
    label1 = 'vs_b' + labelp
    listx = np.linspace(-0.05*1e3,-2*zp*1e3,N)
else:
    z0 = zp*1e3
    # z0 = 0.06*1e3
    labelx = 'E [meV]'   
    title5 = title5 + ', ' + 'z = %inm' %(z0)
    label1 = 'vs_E' + labelp
    listx = np.linspace(1,80,N)
    
title = title1 + ', '  + title2 + '\n'  + title4 + '\n' + title5

def function_num_re(b_nano,energy0):
    omegac0 = energy0*1e-3/aux 
    b = b_nano*1e-3 

    px_f,py_f,pz_f = dipole_moment_num(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,px,py,pz,b,zp)
    # print(rta)
    rta = px_f + py_f + pz_f 
    
    return rta.real

def function_num_im(b_nano,energy0):
    omegac0 = energy0*1e-3/aux 
    b = b_nano*1e-3 

    px_f,py_f,pz_f = dipole_moment_num(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,px,py,pz,b,zp)
    # print(rta)
    rta = px_f + py_f + pz_f 
    
    return rta.imag

def function_num_abs(b_nano,energy0):
    omegac0 = energy0*1e-3/aux 
    b = b_nano*1e-3 

    px_f,py_f,pz_f = dipole_moment_num(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,px,py,pz,b,zp)
    # print(rta)
    rta = px_f + py_f + pz_f 
    
    return np.abs(rta)

def function_ana_re(b_nano,energy0):
    omegac0 = energy0*1e-3/aux 
    b = b_nano*1e-3 

    px_f,py_f,pz_f  = dipole_moment_ana(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,px,py,pz,b,zp)
    
    rta = px_f + py_f + pz_f 
    
    return rta.real

def function_ana_im(b_nano,energy0):
    omegac0 = energy0*1e-3/aux 
    b = b_nano*1e-3 

    px_f,py_f,pz_f  = dipole_moment_ana(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,px,py,pz,b,zp)
    
    rta = px_f + py_f + pz_f 
    
    return rta.imag

def function_ana_abs(b_nano,energy0):
    omegac0 = energy0*1e-3/aux 
    b = b_nano*1e-3 

    px_f,py_f,pz_f  = dipole_moment_ana(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,px,py,pz,b,zp)
    
    rta = px_f + py_f + pz_f 
    
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
ms = 5
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

if plot_vs_b == 1: 
    
    listy_im_ana = []
    listy_im_num = []
    
    listy_re_ana = []
    listy_re_num = []
    
    for value in listx: 

        y_re_ana = function_ana_re(value,E0)       
        y_re_num = function_num_re(value,E0)
        
        listy_re_ana.append(y_re_ana)
        listy_re_num.append(y_re_num)

        y_im_ana = function_ana_im(value,E0)       
        y_im_num = function_num_im(value,E0)
        
        listy_im_ana.append(y_im_ana)
        listy_im_num.append(y_im_num)

    graph(title,labelx,r'$|p|/\alpha_{eff}$',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    plt.plot(listx,listy_re_ana,'.',ms = ms,color = 'purple',label = 'pole aprox')
    plt.plot(listx,listy_re_num,'.-',ms = 3,color = 'lightseagreen',label = 'sin aprox')
    plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
    plt.tight_layout()
    # plt.yscale('log')
    os.chdir(path_save)
    plt.savefig( 'p_tot' + label1 + '.png', format='png')   
    
if plot_vs_E == 1: 
 
    listy_im_ana = []
    listy_im_num = []
    
    listy_re_ana = []
    listy_re_num = []

    listy_abs_ana = []
    listy_abs_num = []
    
    for value in listx: 

        y_re_ana = function_ana_re(z0,value)       
        y_re_num = function_num_re(z0,value)
        
        listy_re_ana.append(y_re_ana)
        listy_re_num.append(y_re_num)

        y_im_ana = function_ana_im(z0,value)       
        y_im_num = function_num_im(z0,value)
        
        listy_im_ana.append(y_im_ana)
        listy_im_num.append(y_im_num)

        y_abs_ana = function_ana_abs(z0,value)       
        y_abs_num = function_num_abs(z0,value)
        
        listy_abs_ana.append(y_abs_ana)
        listy_abs_num.append(y_abs_num)


    graph(title,labelx,r'Re$(p/\alpha_{eff}$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    plt.plot(listx,listy_re_ana,'.',ms = ms,color = 'purple',label = 'pole aprox')
    plt.plot(listx,listy_re_num,'.-',ms = 3,color = 'lightseagreen',label = 'sin aprox')
    plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
    plt.tight_layout()
#    plt.yscale('log')
    os.chdir(path_save)
    plt.savefig( 'Re_p' + label1 + '.png', format='png')   



    graph(title,labelx,r'Im$(p/\alpha_{eff}$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    plt.plot(listx,listy_im_ana,'.',ms = ms,color = 'purple',label = 'pole aprox')
    plt.plot(listx,listy_im_num,'.-',ms = 3,color = 'lightseagreen',label = 'sin aprox')
    plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
    plt.tight_layout()
#    plt.yscale('log')
    os.chdir(path_save)
    plt.savefig( 'Im_p' + label1 + '.png', format='png')   


    graph(title,labelx,r'$|p/\alpha_{eff}$|',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    plt.plot(listx,listy_abs_ana,'.',ms = ms,color = 'purple',label = 'pole aprox')
    plt.plot(listx,listy_abs_num,'.-',ms = 3,color = 'lightseagreen',label = 'sin aprox')
    plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
    plt.tight_layout()
#    plt.yscale('log')
    os.chdir(path_save)
    plt.savefig( 'Abs_p' + label1 + '.png', format='png')   

#%%
