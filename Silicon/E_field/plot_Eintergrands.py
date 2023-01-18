
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
path_ctes =  path_basic.replace('/' + 'E_field','')
path_save = path_basic + '/' + 'E_integrands_pols'
if not os.path.exists(path_save):
    print('Creating folder to save graphs')
    os.mkdir(path_save)

err = 'E_integrand_pols.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from E_integrand_pols import Ex_pols_num, Ex_pols_pole_aprox, Ey_pols_num, Ey_pols_pole_aprox
except ModuleNotFoundError:
    print(err)

try:
    sys.path.insert(1, path_ctes)
    from constants import constantes
except ModuleNotFoundError:
    print('constants.py no se encuentra en ' + path_ctes)

pi,hb,c,alfac,mu1,mu2,mu3 = constantes()
aux = c*hb
aux2 = 1e12/c

#%%

print('Definir parametros del problema')

#v = c/int_v
#omega = 0.7*1e12
cota = 75 #nanometros

epsilon2 = 12
d = 20*1e-3
zp = 0.05
z0 = zp
phi0 = 0

px,py,pz = 1,0,0

title0 = r'px = %i, py = %i, pz = %i, $\varphi$=%i' %(px,py,pz,phi0)
title1 = r'z = %inm, $z_p$=%inm, $\epsilon_2$ = %i, d = %inm' %(z0*1e3,zp*1e3,epsilon2,d*1e3)
labelp = r'_d%inm_px%i_py%i_pz%i' %(d*1e3,px,py,pz)
title = title1 + '\n' + title0

N = 200
omega_omegaWG0 = 2.8 # meV 
R0 = 50 #nanos

# z0 = 0.06*1e3
labelx = r'$\alpha_\parallel$'   
title = title + ', ' + 'R = %inm' %(R0)  + ', ' + r'$\omega/\omega_{WG}$ = %.2f' %(omega_omegaWG0)
label1 = 'vs_alpha_parallel' + labelp
listx = np.linspace(1.001,np.sqrt(epsilon2),N)
listx = np.linspace(1.001,np.sqrt(epsilon2),N)
listx = np.linspace(0.001,np.sqrt(epsilon2),N)

omegaWG = (np.pi/d)/(np.sqrt(epsilon2-1)) # el c esta metido en aux

#%%

def k_parallel_air(omega):

    return omega/c

def k_parallel_medium(omega,epsilon2):

    return omega*np.sqrt(epsilon2)/c



def Ex_real_num(alpha_parallel):
#    omega0 = energy0*1e-3/aux 
#    omega_omegaWG = omega0/omegaWG

    rta = Ex_pols_num(omega_omegaWG0,d,R0,phi0,z0,zp,px,py,pz,alpha_parallel)
    
    return rta.real

def Ex_imag_num(alpha_parallel):
#    omega0 = energy0*1e-3/aux    
#    omega_omegaWG = omega0/omegaWG
    
    rta = Ex_pols_num(omega_omegaWG0,d,R0,phi0,z0,zp,px,py,pz,alpha_parallel)
    return rta.imag

def Ey_real_num(alpha_parallel):
#    omega0 = energy0*1e-3/aux    
#    omega_omegaWG = omega0/omegaWG

    rta = Ey_pols_num(omega_omegaWG0,d,R0,phi0,z0,zp,px,py,pz,alpha_parallel)
    
    return rta.real

def Ey_imag_num(alpha_parallel):
#    omega0 = energy0*1e-3/aux    
#    omega_omegaWG = omega0/omegaWG
    
    rta = Ey_pols_num(omega_omegaWG0,d,R0,phi0,z0,zp,px,py,pz,alpha_parallel)
    return rta.imag

#%%


def Ex_real_pole_aprox(alpha_parallel):
#    omega0 = energy0*1e-3/aux 
#    omega_omegaWG = omega0/omegaWG

    rta = Ex_pols_pole_aprox(omega_omegaWG0,d,R0,phi0,z0,zp,px,py,pz,alpha_parallel)
    
    return rta.real

def Ex_imag_pole_aprox(alpha_parallel):
#    omega0 = energy0*1e-3/aux 
#    omega_omegaWG = omega0/omegaWG
    
    rta = Ex_pols_pole_aprox(omega_omegaWG0,d,R0,phi0,z0,zp,px,py,pz,alpha_parallel)
    return rta.imag


def Ey_real_pole_aprox(alpha_parallel):
#    omega0 = energy0*1e-3/aux 
#    omega_omegaWG = omega0/omegaW

    rta = Ey_pols_pole_aprox(omega_omegaWG0,d,R0,phi0,z0,zp,px,py,pz,alpha_parallel)
    
    return rta.real

def Ey_imag_pole_aprox(alpha_parallel):
#    omega0 = energy0*1e-3/aux 
#    omega_omegaWG = omega0/omegaWG
    
    rta = Ey_pols_pole_aprox(omega_omegaWG0,d,R0,phi0,z0,zp,px,py,pz,alpha_parallel)
    return rta.imag

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



listy_re_numG1 = []
listy_im_numG1 = []
   
listy_re_pole_aproxG1 = []
listy_im_pole_aproxG1 = []
   

listy_re_numG2 = []
listy_im_numG2 = []
   
listy_re_pole_aproxG2 = []
listy_im_pole_aproxG2 = []
   

for value in listx: 

    y_re_numG1 = Ex_real_num(value)
    y_im_numG1 = Ex_imag_num(value)  
    
    
    y_re_pole_aproxG1 = Ex_real_pole_aprox(value)
    y_im_pole_aproxG1 = Ex_imag_pole_aprox(value)  
    
    
    listy_re_numG1.append(y_re_numG1)
    listy_im_numG1.append(y_im_numG1)   
    
    listy_re_pole_aproxG1.append(y_re_pole_aproxG1)
    listy_im_pole_aproxG1.append(y_im_pole_aproxG1)      


    y_re_numG2 = Ey_real_num(value)
    y_im_numG2 = Ey_imag_num(value)  
    
    
    y_re_pole_aproxG2 = Ey_real_pole_aprox(value)
    y_im_pole_aproxG2 = Ey_imag_pole_aprox(value)  
    
    
    listy_re_numG2.append(y_re_numG2)
    listy_im_numG2.append(y_im_numG2)   
    
    listy_re_pole_aproxG2.append(y_re_pole_aproxG2)
    listy_im_pole_aproxG2.append(y_im_pole_aproxG2)   
    

graph(title,labelx,'Integrand Re($E_x$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(listx,listy_re_pole_aproxG1,'.',ms = ms,color = 'darkred',label = 'PP numerical')
plt.plot(listx,listy_re_numG1,'.-',ms = ms,color = 'lightseagreen',label = 'full numerical')
plt.legend(loc = 'lower left',markerscale=mk,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
os.chdir(path_save)
plt.savefig( 'Re_Ex' + label1 + '.png', format='png')   

graph(title,labelx,'Integrand Im($E_x$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(listx,listy_im_pole_aproxG1,'.',ms = ms,color = 'darkred',label = 'PP numerical')
plt.plot(listx,listy_im_numG1,'.-',ms = ms,color = 'lightseagreen',label = 'full numerical')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
os.chdir(path_save)
plt.savefig( 'Im_Ex' + label1 + '.png', format='png')   


graph(title,labelx,'Integrand Re($E_y$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(listx,listy_re_pole_aproxG2,'.',ms = ms,color = 'darkred',label = 'PP numerical')
plt.plot(listx,listy_re_numG2,'.-',ms = ms,color = 'lightseagreen',label = 'full numerical')
plt.legend(loc = 'lower left',markerscale=mk,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
os.chdir(path_save)
plt.savefig( 'Re_Ey' + label1 + '.png', format='png')   

graph(title,labelx,'Integrand Im($E_y$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(listx,listy_im_pole_aproxG2,'.',ms = ms,color = 'darkred',label = 'PP numerical')
plt.plot(listx,listy_im_numG2,'.-',ms = ms,color = 'lightseagreen',label = 'full numerical')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
os.chdir(path_save)
plt.savefig( 'Im_Ey' + label1 + '.png', format='png')   



    
#%%
