
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

plot_vs_z = 1
plot_vs_E = 0

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/potential_field/potential/potential_and_electric_field_with_dipole_moment_formula','')
path_save = path_basic + '/' + 'electric_field'
if not os.path.exists(path_save):
    print('Creating folder to save graphs')
    os.mkdir(path_save)

err = 'fieldE_direct_numerical.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from electric_field import electric_field_tot_ana, electric_field_tot_num
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

int_v = 300
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
py = 0
pz = 0

omega0THz = 0.05
omega0 = omega0THz*1e12 
energy0_pol = omega0*hb
#R = 10 # 10 nm en unidades de micrometros
kappa_factor_omega0 = 0.1
kappa_r_factor= 0.5
 
title1 = r'v = c/%i $\mu$m/s, $\kappa$ = %.2f$\omega_0$, $\kappa_r$ = %.2f$\kappa$, $E_0$=%.2feV' %(int_v, kappa_factor_omega0, kappa_r_factor, energy0_pol)   
title2 = r'$\hbar\mu$ = %.2feV, $\hbar\gamma$ = %.4feV' %(hbmu,hbgama) 
title3 = r'$z_p$=%inm, px=%i, py=%i, pz=%i' %(zp*1e3,px,py,pz)
title4 = r'$z_p$=%inm, px=%i, py=%i, pz=%i, xD=%inm, yD=%inm' %(zp*1e3,px,py,pz,xD*1e3,yD*1e3)
title5 = r'b = %inm, zD=%inm, x=%inm, y=%inm' %(b*1e3,zD*1e3,x0*1e3,y0*1e3)
labelp = r'_px%i_py%i_pz%i' %(px,py,pz)

N = 100
if plot_vs_z == 1: 
    E0 = 0.46 #eV 
    labelx = 'z [nm]'  
    title5 = title5 + ', ' + 'E = %.2f eV' %(E0)
    label1 = 'vs_z' + labelp
    listx = np.linspace(zp*1e3,-2*zp*1e3,N)
else:
    z0 = zp*1e3
    # z0 = 0.06*1e3
    labelx = 'E [eV]'   
    title5 = title5 + ', ' + 'z = %inm' %(z0)
    label1 = 'vs_E' + labelp
    listx = np.linspace(0.46,2,N)
    
title = title1 + '\n'  + title2 + ', ' + title3 + '\n'  + title4 + '\n' + title5

def function_real_num(z_nano,energy0):
    omegac0 = energy0/aux 
    z = z_nano*1e-3 

    rta = electric_field_tot_num(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,px,py,pz,b,zp,x0,y0,z,xD,yD,zD,omega0,kappa_factor_omega0,kappa_r_factor)
    # print(rta)
    return rta

def function_real_ana(z_nano,energy0):
    omegac0 = energy0/aux 
    z = z_nano*1e-3 

    rta = electric_field_tot_ana(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,px,py,pz,b,zp,x0,y0,z,xD,yD,zD,omega0,kappa_factor_omega0,kappa_r_factor)
    
    return rta

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

if plot_vs_z == 1: 
    listy_re_ana = []
    listy_re_num = []
    
    for value in listx: 

        y_re_ana = function_real_ana(value,E0)     
        y_re_num = function_real_num(value,E0)
        
        listy_re_ana.append(y_re_ana)
        listy_re_num.append(y_re_num)  

    
    graph(title,labelx,'$|E|^2$',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    plt.plot(listx,listy_re_ana,'.',ms = ms,color = 'purple',label = 'pole aprox')
    plt.plot(listx,listy_re_num,'.-',ms = 3,color = 'lightseagreen',label = 'sin aprox')
    plt.legend(loc = 'best',markerscale=1,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'Etot' + label1 + '.png', format='png')   

    
if plot_vs_E == 1: 
    
    listy_re_ana = []
    listy_re_num = []    
    
    for value in listx: 

        y_re_ana = function_real_ana(z0,value)     
        y_re_num = function_real_num(z0,value)
        
        listy_re_ana.append(y_re_ana)
        listy_re_num.append(y_re_num)
 
    graph(title,labelx,'$|E|^2$',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    plt.plot(listx,listy_re_ana,'.',ms = ms,color = 'purple',label = 'pole aprox')
    plt.plot(listx,listy_re_num,'.-',ms = 3,color = 'lightseagreen',label = 'sin aprox')
    plt.legend(loc = 'best',markerscale=1,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'Etot' + label1 + '.png', format='png')   


#%%
