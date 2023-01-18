
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

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/potential_field/potential/potential_and_electric_field_with_dipole_moment_formula','')
path_save = path_basic + '/' + 'dipole_moment_integrands'
if not os.path.exists(path_save):
    print('Creating folder to save graphs')
    os.mkdir(path_save)

err = 'fieldE_direct_numerical.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from dipole_moment_integrands import INT1_pole_aprox, INT2_pole_aprox, INT3_pole_aprox, INT1_num, INT2_num, INT3_num
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
epsi1,epsi3 = 1,1
hbmu,hbgama = 0.3,0.0001
zp = 0.05
b = -0.01


d_nano = 10

#omega0THz = 65
#omega0 = omega0THz*1e12 
energy0_pol = 43
omega0 = energy0_pol*1e-3/hb 
#R = 10 # 10 nm en unidades de micrometros
kappa_factor_omega0 = 0.1
kappa_r_factor= 0.5

energy0 = 43
 
title1 = r'v = c/%i, $\kappa$ = %.2f$\omega_0$, $\kappa_r$ = %.2f$\kappa$, $\hbar\omega_0$=%i meV' %(int_v, kappa_factor_omega0, kappa_r_factor, energy0_pol)     
#title2 = r'$\hbar\mu$ = %.2feV, $\hbar\gamma$ = %.4feV' %(hbmu,hbgama) 
#title3 = r'$z_p$=%inm, px=%i, py=%i, pz=%i' %(zp*1e3,px,py,pz)
title4 = r'$z_p$=%i nm, b = %i nm, \hbar\omega = %i meV' %(zp*1e3,b*1e3,energy0)
labelp = r'_E0_%i' %(energy0_pol)

N = 100

    # z0 = 0.06*1e3
labelx = r'$\alpha_y$'   
label1 = 'vs_E' + labelp
listx = np.linspace(0.001,80,N)

title = title1  + '\n'  + title4 

def Int1_num(energy0,u):
    omegac0 = energy0*1e-3/aux 

    px_f = INT1_num(omegac0,epsi1,epsi3,d_nano,int_v,b,zp,omega0,kappa_factor_omega0,kappa_r_factor,u)
    # print(rta)

    return px_f


def Int2_num(energy0,u):
    omegac0 = energy0*1e-3/aux 

    px_f = INT2_num(omegac0,epsi1,epsi3,d_nano,int_v,b,zp,omega0,kappa_factor_omega0,kappa_r_factor,u)
    # print(rta)

    return px_f


def Int3_num(energy0,u):
    omegac0 = energy0*1e-3/aux 

    px_f = INT3_num(omegac0,epsi1,epsi3,d_nano,int_v,b,zp,omega0,kappa_factor_omega0,kappa_r_factor,u)
    # print(rta)

    return px_f

#%%
    
def Int1_pole_aprox(energy0,u):
    omegac0 = energy0*1e-3/aux 

    px_f = INT1_pole_aprox(omegac0,epsi1,epsi3,d_nano,int_v,b,zp,omega0,kappa_factor_omega0,kappa_r_factor,u)
    # print(rta)

    return px_f


def Int2_pole_aprox(energy0,u):
    omegac0 = energy0*1e-3/aux 

    px_f = INT2_pole_aprox(omegac0,epsi1,epsi3,d_nano,int_v,b,zp,omega0,kappa_factor_omega0,kappa_r_factor,u)
    # print(rta)

    return px_f


def Int3_pole_aprox(energy0,u):
    omegac0 = energy0*1e-3/aux 

    px_f = INT3_pole_aprox(omegac0,epsi1,epsi3,d_nano,int_v,b,zp,omega0,kappa_factor_omega0,kappa_r_factor,u)
    # print(rta)

    return px_f


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


listy_re_num1 = []
listy_re_num2 = []
listy_re_num3 = []

listy_im_num1 = []
listy_im_num2 = []
listy_im_num3 = []


listy_re_pole_aprox1 = []
listy_re_pole_aprox2 = []
listy_re_pole_aprox3 = []

listy_im_pole_aprox1 = []
listy_im_pole_aprox2 = []
listy_im_pole_aprox3 = []


for value in listx: 


    y_num1 = Int1_num(energy0,value)
    y_pole1 = Int1_pole_aprox(energy0,value)
    


    y_num2 = Int2_num(energy0,value)
    y_pole2 = Int2_pole_aprox(energy0,value)



    y_num3 = Int3_num(energy0,value)
    y_pole3 = Int3_pole_aprox(energy0,value)







    listy_re_num1.append(np.real(y_num1))
    listy_re_num2.append(np.real(y_num2))
    listy_re_num3.append(np.real(y_num3))

    listy_im_num1.append(np.imag(y_num1))
    listy_im_num2.append(np.imag(y_num2))
    listy_im_num3.append(np.imag(y_num3))



    listy_re_pole_aprox1.append(np.real(y_pole1))
    listy_re_pole_aprox2.append(np.real(y_pole2))
    listy_re_pole_aprox3.append(np.real(y_pole3))

    listy_im_pole_aprox1.append(np.imag(y_pole1))
    listy_im_pole_aprox2.append(np.imag(y_pole2))
    listy_im_pole_aprox3.append(np.imag(y_pole3))


#%%
graph(title,labelx,r'Real(I1)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#    plt.plot(listx,listy_re_ana,'.',ms = ms,color = 'purple',label = 'analytical')

plt.plot(listx,listy_re_num1,'.',ms = ms,color = 'lightseagreen',label = 'full numerical')
plt.plot(listx,listy_re_pole_aprox1,'.-',ms = 3,color = 'darkred',label = 'PP numerical')
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
#plt.yscale('log')
os.chdir(path_save)
plt.savefig( 'I1_real' + label1 + '.png', format='png')   


graph(title,labelx,r'Imag(I1)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#    plt.plot(listx,listy_re_ana,'.',ms = ms,color = 'purple',label = 'analytical')

plt.plot(listx,listy_im_num1,'.',ms = ms,color = 'lightseagreen',label = 'full numerical')
plt.plot(listx,listy_im_pole_aprox1,'.-',ms = 3,color = 'darkred',label = 'PP numerical')
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
#plt.yscale('log')
os.chdir(path_save)
plt.savefig( 'I1_im' + label1 + '.png', format='png')   


#%%


graph(title,labelx,r'Real(I2)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#    plt.plot(listx,listy_re_ana,'.',ms = ms,color = 'purple',label = 'analytical')

plt.plot(listx,listy_re_num2,'.',ms = ms,color = 'lightseagreen',label = 'full numerical')
plt.plot(listx,listy_re_pole_aprox2,'.-',ms = 3,color = 'darkred',label = 'PP numerical')
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
#plt.yscale('log')
os.chdir(path_save)
plt.savefig( 'I2_real' + label1 + '.png', format='png')   


graph(title,labelx,r'Imag(I2)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#    plt.plot(listx,listy_re_ana,'.',ms = ms,color = 'purple',label = 'analytical')

plt.plot(listx,listy_im_num2,'.',ms = ms,color = 'lightseagreen',label = 'full numerical')
plt.plot(listx,listy_im_pole_aprox2,'.-',ms = 3,color = 'darkred',label = 'PP numerical')
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
#plt.yscale('log')
os.chdir(path_save)
plt.savefig( 'I2_im' + label1 + '.png', format='png')   

#%%


graph(title,labelx,r'Real(I3)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#    plt.plot(listx,listy_re_ana,'.',ms = ms,color = 'purple',label = 'analytical')

plt.plot(listx,listy_re_num3,'.',ms = ms,color = 'lightseagreen',label = 'full numerical')
plt.plot(listx,listy_re_pole_aprox3,'.-',ms = 3,color = 'darkred',label = 'PP numerical')
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
#plt.yscale('log')
os.chdir(path_save)
plt.savefig( 'I3_real' + label1 + '.png', format='png')   


graph(title,labelx,r'Imag(I3)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#    plt.plot(listx,listy_re_ana,'.',ms = ms,color = 'purple',label = 'analytical')

plt.plot(listx,listy_im_num3,'.',ms = ms,color = 'lightseagreen',label = 'full numerical')
plt.plot(listx,listy_im_pole_aprox3,'.-',ms = 3,color = 'darkred',label = 'PP numerical')
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
#plt.yscale('log')
os.chdir(path_save)
plt.savefig( 'I3_im' + label1 + '.png', format='png')   
