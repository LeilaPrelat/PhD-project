
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


comparar_otro_hBN = 0

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/potential_field/potential_and_electric_field_with_dipole_moment_formula','')
path_save = path_basic + '/' + 'dipole_moment'
if not os.path.exists(path_save):
    print('Creating folder to save graphs')
    os.mkdir(path_save)

err = 'fieldE_direct_numerical.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from dipole_moment_dif_sign import dipole_moment_ana_resonance_v1, dipole_moment_ana_resonance_v2, dipole_moment_pole_aprox_resonance_v1, dipole_moment_pole_aprox_resonance_v2, dipole_moment_num_resonance
except ModuleNotFoundError:
    print(err)

try:
    sys.path.insert(1, path_constants)
    from Silica_epsilon import epsilon_Silica
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

zp = 0.05
b = -0.01

d_thickness_disk_nano = 0.1
D_disk_nano = 50 

d_nano_film = 1

#omega0THz = 65
#omega0 = omega0THz*1e12 

 
    
#title2 = r'$\hbar\mu$ = %.2feV, $\hbar\gamma$ = %.4feV' %(hbmu,hbgama) 
#title3 = r'$z_p$=%inm, px=%i, py=%i, pz=%i' %(zp*1e3,px,py,pz)
title1 = r'$v$ = c/%i, $z_{\rm p}$ = %i nm, $b$ = %i nm, $d_{\rm film}$ = %.2f nm' %(int_v, zp*1e3,b*1e3,d_nano_film)
title2 = r'$d_{\rm disk}$ = %.2f nm, $D$ = %.2f nm' %(d_thickness_disk_nano,D_disk_nano)
labelp = r'_res_d%.2fnm' %(d_nano_film)

N = 75

    # z0 = 0.06*1e3
labelx = r'$\hbar\omega$ [eV]'   
labelp = 'vs_E' + labelp

x1 = 0.09260651629072682 
x2 = 0.10112781954887218
x3 = 0.17030075187969923
x4 = 0.19937343358395992

x0 = 0.08684904004367106 
xf =  0.8064559363534132

listx = np.linspace(0.09,0.195,N)
listx = np.linspace(x1+1e-3,0.1975,N)
listx = np.linspace(x1 + 1e-3, x4 - 1e-3,N)
#listx = np.linspace(x0+1e-3,xf - 1e-3,N)
#listx = np.linspace(x1,x2,N)

## entre x1 y x2, y entre x3 y x4 ##

title =  title1 + '\n' + title2  



#%%

def function_num(energy0):
    omegac0 = energy0/aux 

    px_f,py_f,pz_f = dipole_moment_num_resonance(omegac0,epsilon_Silica,d_nano_film,d_thickness_disk_nano,D_disk_nano,int_v,b,zp)
    # print(rta)
    rta = np.abs(px_f)**2 + np.abs(py_f)**2 + np.abs(pz_f)**2 

    return np.sqrt(rta)


#%%

def function_ana_v1(energy0):
    omegac0 = energy0/aux 

    px_f,py_f,pz_f  = dipole_moment_ana_resonance_v1(omegac0,epsilon_Silica,d_nano_film,d_thickness_disk_nano,D_disk_nano,int_v,b,zp)
    
    rta = np.abs(px_f)**2 + np.abs(py_f)**2 + np.abs(pz_f)**2 
    
    return np.sqrt(rta)




def function_ana_v2(energy0):
    omegac0 = energy0/aux 

    px_f,py_f,pz_f  = dipole_moment_ana_resonance_v2(omegac0,epsilon_Silica,d_nano_film,d_thickness_disk_nano,D_disk_nano,int_v,b,zp)
    
    rta = np.abs(px_f)**2 + np.abs(py_f)**2 + np.abs(pz_f)**2 
    
    return np.sqrt(rta)

#%%

def function_pole_approx_v1(energy0):
    omegac0 = energy0/aux 

    px_f,py_f,pz_f  = dipole_moment_pole_aprox_resonance_v1(omegac0,epsilon_Silica,d_nano_film,d_thickness_disk_nano,D_disk_nano,int_v,b,zp)
    
    rta = np.abs(px_f)**2 + np.abs(py_f)**2 + np.abs(pz_f)**2 
    
    return np.sqrt(rta)




def function_pole_approx_v2(energy0):
    omegac0 = energy0/aux 

    px_f,py_f,pz_f  = dipole_moment_pole_aprox_resonance_v2(omegac0,epsilon_Silica,d_nano_film,d_thickness_disk_nano,D_disk_nano,int_v,b,zp)
    
    rta = np.abs(px_f)**2 + np.abs(py_f)**2 + np.abs(pz_f)**2 
    
    return np.sqrt(rta)


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

listy_ana_v1 = []
listy_ana_v2 = []

listy_pole_v1 = []
listy_pole_v2 = []

listy_num = []

for value in listx: 
    

    y_ana_v1 = function_ana_v1(value)     
    y_ana_v2 = function_ana_v2(value)  
    
    y_pole_v1 = function_pole_approx_v1(value)
    y_pole_v2 = function_pole_approx_v2(value)
    
    y_num = function_num(value)
    
#    listy_re_ana.append(y_re_ana)
    
    listy_ana_v1.append(y_ana_v1)
    listy_ana_v2.append(y_ana_v2)
    listy_num.append(y_num)
    
    listy_pole_v1.append(y_pole_v1)
    listy_pole_v2.append(y_pole_v2)

#Emax = listx[np.argmax(listy_re_num)]
#print(Emax)

#%%

label1 = r'$r_{\rm p} = k_\parallel/(k_\parallel - k_{\rm p})$'    
label2 = r'$r_{\rm p} = k_{\rm p}/(k_\parallel - k_{\rm p})$'    



if comparar_otro_hBN == 1:
    
    os.chdir(path_save)
    label1 = 'vs_E_res_d%.2fnm' %(d_nano_film) 
       
    tabla = np.loadtxt( 'dip_mom_hBN_num_' + label1 + '.txt', delimiter='\t', skiprows=1)
    tabla = np.transpose(tabla)
    [listx_num_S,listy_re_num_S] = tabla 
    
    
    tabla = np.loadtxt( 'dip_mom_hBN_ana_' + label1 + '.txt', delimiter='\t', skiprows=1)
    tabla = np.transpose(tabla)
    [listx_ana_S,listy_re_ana2_S] = tabla 
    
    
    tabla = np.loadtxt( 'dip_mom_hBN_pole_aprox_' + label1 + '.txt', delimiter='\t', skiprows=1)
    tabla = np.transpose(tabla)
    [listx_pole_S,listy_re_pole_S] = tabla 



graph(title,labelx,r'$|p|$/(e/2$\pi v$) ($\mu$m$^{2}$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#plt.plot(listx[0:-2],listy_ana_v1[0:-2],'.-',ms = ms,color = 'purple',label = 'PP ana')
#plt.plot(listx,listy_ana_v2,'--',ms = ms,color = 'purple',label = 'PP ana 2 ' +  label2)
plt.plot(listx,listy_num,'.',ms = ms,color = 'darkred',label = 'full numerical')
plt.plot(listx,listy_pole_v1,'.-',ms = 3,color = 'darkred',label = 'PP num' )
plt.plot(listx,listy_pole_v2,'--',ms = 3,color = 'purple',label = 'PP num 2'  )
if comparar_otro_hBN == 1:
    plt.plot(listx_ana_S,listy_re_ana2_S,'--',ms = ms,color = 'blue',label = 'PP ana sphere')
    plt.plot(listx_num_S,listy_re_num_S,'--',ms = ms,color = 'lightseagreen',label = 'full num sphere')

ejey_aux = np.linspace(np.min([np.min(listy_pole_v1),np.min(listy_num)]), np.max([np.max(listy_pole_v1[0:-2]),np.max(listy_num[0:-2])]) , 10)
for x in [x1,x2,x3]:
    plt.plot(x*np.ones(10), ejey_aux,'--',color = 'grey' )

if d_nano_film == 0.1:
    plt.ylim(1e-1,1e5)
plt.plot(listx[-1]*np.ones(10), ejey_aux,'--',color = 'grey' )    
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0.05,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
plt.yscale('log')
os.chdir(path_save)
plt.savefig( 'p_tot' + labelp + '.png', format='png')   



#%%
