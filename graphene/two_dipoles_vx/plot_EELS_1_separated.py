
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
path_constants =  path_basic.replace('/two_dipoles_vx','')
path_save = path_basic + '/' + 'EELS_1_separated'
if not os.path.exists(path_save):
    print('Creating folder to save graphs')
    os.mkdir(path_save)

err = 'EELS_1_separated.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from EELS_1_separated import EELS_1_ana, EELS_2_ana, EELS_parallel_f_dir
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



#omega0THz = 65
#omega0 = omega0THz*1e12 
energy0_pol1 = 40
omega01 = energy0_pol1*1e-3/hb 

energy0_pol2 = 45
omega02 = energy0_pol2*1e-3/hb 



kappa_factor_omega01 = 0.1
kappa_r_factor1 = 0.5

kappa_factor_omega02 = 0.1
kappa_r_factor2 = 0.5

L = 1000*1e-3

x1 = -0.01
x2 = 0.01

z1 = 0
z2 = 0
 
title1 = r'$\hbar\omega_{01}$=%i meV, $\hbar\omega_{02}$=%i meV, v = c/%i, L = %i nm' %(energy0_pol1,energy0_pol2,int_v,L*1e3)          
#title2 = r'$\hbar\mu$ = %.2feV, $\hbar\gamma$ = %.4feV' %(hbmu,hbgama) 
#title3 = r'$z_p$=%inm, px=%i, py=%i, pz=%i' %(zp*1e3,px,py,pz)
title4 = r'$x_1$ = %i nm, $x_2$ = %i nm, b = %i nm, $z_p$ = %i nm' %(x1*1e3,x2*1e3,b*1e3,zp*1e3)


N = 100

    # z0 = 0.06*1e3
labelx = r'$\hbar\omega$ [meV]'   
label1 = '_vs_E_x1_%i_x2_%i' %(x1*1e3,x2*1e3) 
listx = np.linspace(43,47,N)

title = title1   + '\n' + title4 

def function_parallel(energy0):
    omegac0 = energy0*1e-3/aux 

    px_f  = EELS_parallel_f_dir(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,b,zp,L)
    # print(rta)
    return np.real(px_f)

def function_dip_1(energy0):
    omegac0 = energy0*1e-3/aux 

    px_f  = EELS_1_ana(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,b,zp,x1,z1,x2,z2,omega01,kappa_factor_omega01,kappa_r_factor1,omega02,kappa_factor_omega02,kappa_r_factor2)
       
    return np.real(px_f)

def function_dip_2(energy0):
    omegac0 = energy0*1e-3/aux 

    px_f  = EELS_2_ana(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,b,zp,x1,z1,x2,z2,omega01,kappa_factor_omega01,kappa_r_factor1,omega02,kappa_factor_omega02,kappa_r_factor2)
       
    return np.real(px_f)


#%%
    
tamfig = (5.5,3.5)
tamlegend = 12
tamletra = 12
tamtitle = 11
tamnum = 10
labelpady = 2
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

listy_re_1 = []
listy_re_2 = []
listy_re_parallel = []

for value in listx: 

#    y_re_ana = function_ana(value)       
    y_1 = function_dip_1(value)  
    y_2 = function_dip_2(value)  
    y_parallel = function_parallel(value)
 
#    listy_re_ana.append(y_re_ana)
    
    listy_re_1.append(np.real(y_1))
    listy_re_2.append(np.real(y_2))
    listy_re_parallel.append(np.real(y_parallel))



#Emax = listx[np.argmax(listy_re_num)]
#print(Emax)

#%%

#mini = np.min([np.min(listy_re_ana),np.min(listy_re_num)])
#print(mini)
hspace = 0.5
wspace = 0.1  
fig,axs = plt.subplots(2,1, sharex=True, facecolor='w', figsize = tamfig)
plt.subplots_adjust(hspace =hspace,wspace = wspace)
labely = r'$\Gamma^{EELS}$ [eV$^{-1}$]'
#graph(title,labelx,r'$\Gamma^{EELS}$ [eV$^{-1}$]',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#    plt.plot(listx,listy_re_ana,'.',ms = ms,color = 'purple',label = 'analytical')
axs[0].plot(listx,np.array(listy_re_parallel)/hb,'.-',ms = ms,color = 'lightseagreen')
axs[1].plot(listx,np.array(listy_re_2)/hb,'.-',ms = ms,color = 'darkred',label = r'$\Gamma_{dip,2}$')
axs[1].plot(listx,np.array(listy_re_1)/hb,'.-',ms = ms,color = 'purple',label = r'$\Gamma_{dip,1}$')


for i in [0,1]:
    axs[i].minorticks_on()
    axs[i].tick_params(labelsize = tamnum,pad = pad)

axs[0].set_ylabel(r'$\Gamma^{EELS}_\parallel$ [eV$^{-1}$]',fontsize = tamletra,labelpad=labelpady) 
axs[1].set_ylabel(r'$\Gamma^{EELS}_{dip}$ [eV$^{-1}$]',fontsize = tamletra,labelpad=labelpady) 
axs[1].set_xlabel(labelx,fontsize = tamletra,labelpad=labelpadx)
axs[0].set_title(title,fontsize = tamnum)
fig.legend(loc = [0.8,0.27],markerscale=mk,fontsize=tamlegend,frameon=0.01,handletextpad=0.2, handlelength=length_marker)

#plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
fig.tight_layout()
#plt.yscale('log')
os.chdir(path_save)
plt.savefig( 'EELS' + label1 + '.png', format='png')   

#%%
