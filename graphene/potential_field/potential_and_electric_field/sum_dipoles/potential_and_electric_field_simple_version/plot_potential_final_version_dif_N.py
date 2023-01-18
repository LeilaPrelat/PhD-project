
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
path_constants =  path_basic.replace('/potential_field/potential_and_electric_field/sum_dipoles/potential_and_electric_field_simple_version','')
path_save = path_basic + '/' + 'sum_potential_dif_Ndip'
if not os.path.exists(path_save):
    print('Creating folder to save graphs')
    os.mkdir(path_save)

err = 'fieldE_direct_numerical.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from potential_final_version_dif_N import potential_final_ana
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
z0 = zp
x0 = 0
a = 50*1e-3 # 50 nanometros0

px = 1
py = 0
pz = 0
Ndip = 15

list_Ndip = [200,1000,1500]
int_v = 10

title5 = r'a = %i nm, v = c/%i' %(a*1e3,int_v)
labelp = r'_px%i_py%i_pz%i_a%i_N%i' %(px,py,pz,a,Ndip)

N = 250
if plot_vs_y == 1: 
    E0 = 43 # meV 
    labelx = 'y [nm]'  

    #### longitud de propagacion del plasmon ####
    
    listx = np.linspace(500,2000,N)

    title5 = title5 + ', ' + 'E = %i meV' %(E0)
    label1 = 'vs_y' + labelp
    
else:
    y0 = 3000 #nanos
    # z0 = 0.06*1e3
    labelx = 'E [meV]'   
    title5 = title5 + ', ' + 'y = %inm' %(y0)
    label1 = 'vs_E' + labelp
    listx = np.linspace(0.0001,2,N)
    listx = np.linspace(20,55,N)


#title2 = r'$\hbar\mu$ = %.2feV, $\hbar\gamma$ = %.4feV' %(hbmu,hbgama) 
title4 = r'z = $z_p$=%i nm, px = %i, py = %i, pz = %i, x = %i nm,' %(zp*1e3,px,py,pz,x0*1e3)
#title = title2 + '\n' + title4 + '\n'  + title5
title = title4 + '\n'  + title5

def function_real_ana(y_nano,energy0,Ndip):
    omegac0 = energy0*1e-3/aux 
    y = y_nano*1e-3 

    rta = potential_final_ana(omegac0,epsi1,epsi2,hbmu,hbgama,x0,y,z0,a,zp,int_v,px,py,pz,Ndip)
    
    return rta.real

def function_imag_ana(y_nano,energy0,Ndip):
    omegac0 = energy0*1e-3/aux 
    y = y_nano*1e-3 
    
    rta = potential_final_ana(omegac0,epsi1,epsi2,hbmu,hbgama,x0,y,z0,a,zp,int_v,px,py,pz,Ndip)
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

if plot_vs_y == 1: 

    listy_re_ana_tot = []
    listy_im_ana_tot = []

    for Ndip in list_Ndip:
        
        print(Ndip)
        
        listy_re_ana = []
        listy_im_ana = []
    
        for value in listx:
        
        
            y_re_ana = function_real_ana(value,E0,Ndip)
            y_im_ana = function_imag_ana(value,E0,Ndip)    
            
            listy_re_ana.append(y_re_ana)
            listy_im_ana.append(y_im_ana)
        
        
        listy_re_ana_tot.append(listy_re_ana)
        listy_im_ana_tot.append(listy_im_ana)
        

if plot_vs_E == 1: 
    
    
    listy_re_ana_tot = []
    listy_im_ana_tot = []

    for Ndip in list_Ndip:
        
        print(Ndip)
        
        listy_re_ana = []
        listy_im_ana = []
    
        for value in listx:
        
        
            y_re_ana = function_real_ana(y0,value,Ndip)
            y_im_ana = function_imag_ana(y0,value,Ndip)    
            
            listy_re_ana.append(y_re_ana)
            listy_im_ana.append(y_im_ana)
        
        
        listy_re_ana_tot.append(listy_re_ana)
        listy_im_ana_tot.append(listy_im_ana)
        

#%%     
    
graph(title,labelx,'Re($\phi$) analytical',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
for k in range(len(list_Ndip)): 
    plt.plot(listx,listy_re_ana_tot[k],'.-',ms = ms,label = 'Ndip = %i' %(list_Ndip[k]))
plt.legend(loc = 'best',markerscale=1.5,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
os.chdir(path_save)
plt.savefig( 'Re_phi' + label1 + '.png', format='png')   

graph(title,labelx,'Im($\phi$) analytical',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
for k in range(len(list_Ndip)): 
    plt.plot(listx,listy_im_ana_tot[k],'.-',ms = ms,label = 'Ndip = %i' %(list_Ndip[k]))
plt.legend(loc = 'best',markerscale=1.5,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
os.chdir(path_save)
plt.savefig( 'Im_phi' + label1 + '.png', format='png')   
    
#%%
