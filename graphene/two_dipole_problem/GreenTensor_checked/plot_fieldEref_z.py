
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: leila

graficar el campo externo directo con la convencion de z hacia abajo

"""

import numpy as np
import sys
import os 
import matplotlib.pyplot as plt
import seaborn as sns

from scipy import special

#%%

save_graphs = 1 #guardar los graficos 2D del campo


plot_vs_E = 1
plot_vs_x = 0
sns.set()    

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/GreenTensor_checked','')
err = 'fieldE_ref_z.py no se encuentra en ' + path_basic

try:
    sys.path.insert(1, path_basic)
    from fieldE_ref_z import Efield_ref_num,Efield_ref_ana,Efield_ref_fresnel
except ModuleNotFoundError:
    print(err)


path_save = path_basic + '/' + 'fieldE_ref_z'

if save_graphs==1:    
    if not os.path.exists(path_save):
        print('Creating folder to save graphs')
        os.mkdir(path_save)

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
epsi1 = 1
epsi2 = 1

zp = 50*1e-3 #micrones # posicion del plano (z apunta hacia abajo)
z = zp
xe = -10*1e-3

hbmu, hbgama = 0.3, 0.0001 #potencial quimico in ev, freq de collision in ev

N = 100

title1 = r'z = %i nm, $x_e$ = $x_1$ = %i nm' %(z*1e3,xe*1e3)
title2 = r'$z_p$ = %i nm, v = c/%i' %(zp*1e3,int_v)
########################
if plot_vs_E == 1:
#    omegaTHz0 = 0.7
#    omegaTHz0 = 1.511
    

    # omegaTHz = 1.663 - 0.152j
    
    R0 = 1000
    listx = np.linspace(25,55,N)

    labelx = '$\hbar\omega$ [meV]'
    label1 = '_vs_E'
    title2 = title2 +  r', x = %i nm' %(R0)


elif plot_vs_x == 1:
#    omegaTHz0 = 0.7
#    omegaTHz0 = 1.511
    

    # omegaTHz = 1.663 - 0.152j
    
    E0 = 40
    listx = np.linspace(150,2000,N)

    labelx = 'x [nm]'
    label1 = '_vs_x'
    title2 = title2 +  r', $\hbar\omega$ = %i meV' %(E0)

title = title1 + '\n' + title2


def Efield_ANA_QE_function(energy0,x_nano):
    
    omegac0 = energy0*1e-3/aux 
    x = x_nano*1e-3
    
    return Efield_ref_ana(omegac0,epsi1,epsi2,hbmu,hbgama,x,z,xe,zp,int_v)

def Efield_NUM_QE_function(energy0,x_nano):
    
    omegac0 = energy0*1e-3/aux 
    x = x_nano*1e-3
    
    return Efield_ref_num(omegac0,epsi1,epsi2,hbmu,hbgama,x,z,xe,zp,int_v)

def Efield_NUM_fresnel_function(energy0,x_nano):
    
    omegac0 = energy0*1e-3/aux 
    x = x_nano*1e-3
    
    return Efield_ref_fresnel(omegac0,epsi1,epsi2,hbmu,hbgama,x,z,xe,zp,int_v)



#%%
    
tamfig = (4.5,3.5)
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

def graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad):
    plt.figure(figsize=tamfig)
    plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
    plt.ylabel(labely,fontsize=tamletra,labelpad =labelpadx)
    plt.tick_params(labelsize = tamnum, pad = pad)
    plt.title(title,fontsize=int(tamtitle*0.9))

    return   

# color_rgb1 = (1.0, 0.8, 0.6)
# color_rgb2 = (0.99, 0.74, 0.71)
# color_rgb3 = (0.529, 0.992, 0.945)
# color_rgb4 = (0.41, 0.21, 0.61)

color_rgb1 = (1.0, 0.01, 0.24)
color_rgb2 = (0.8, 0.58, 0.46)
color_rgb3 = (0.55, 0.71, 0)
color_rgb4 = (0.6, 0.73, 0.89)

#%%


total1_re = []
total1_im = []


total2_re = []
total2_im = []

total3_re = []
total3_im = []


if plot_vs_x == 1:
    for value in listx: 
        
        tot2 = Efield_ANA_QE_function(E0,value)
        total2_re.append(tot2.real)
        total2_im.append(tot2.imag)
    
            
        tot1 = Efield_NUM_QE_function(E0,value)   
        total1_re.append(tot1.real)
        total1_im.append(tot1.imag)


        tot3 = Efield_NUM_fresnel_function(E0,value)   
        total3_re.append(tot3.real)
        total3_im.append(tot3.imag)


        print(value)
    

if plot_vs_E == 1:
    for value in listx:
        
        tot2 = Efield_ANA_QE_function(value,R0)
        total2_re.append(tot2.real)
        total2_im.append(tot2.imag)
                

        tot1 = Efield_NUM_QE_function(value,R0)   
        total1_re.append(tot1.real)
        total1_im.append(tot1.imag)


        tot3 = Efield_NUM_fresnel_function(value,R0)   
        total3_re.append(tot3.real)
        total3_im.append(tot3.imag)


        print(value)
    
#%%


graph(title,labelx,r'Re($E_{ref}$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(listx,total2_re,'.-',ms = ms,color = 'purple',label = 'analytical PP')
# plt.yscale('log')
plt.plot(listx,total1_re,'.',ms = ms,color = 'darkred',label = 'numerical PP')
plt.plot(listx,total3_re,'.-',ms = ms,color = 'lightseagreen', alpha = 0.5,label = 'full numerical')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
if save_graphs==1:
    os.chdir(path_save)
    plt.tight_layout()
    plt.savefig( 'Re_Eref' + label1 + '.png', format='png')

graph(title,labelx,r'Im($E_{ref}$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(listx,total2_im,'.-',ms = ms,color = 'purple',label = 'analytical PP')
# plt.yscale('log')
plt.plot(listx,total1_im,'.',ms = ms,color = 'darkred',label = 'numerical PP')
plt.plot(listx,total3_im,'.-',ms = ms,color = 'lightseagreen', alpha = 0.5,label = 'full numerical')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
if save_graphs==1:
    plt.tight_layout()
    plt.savefig( 'Im_Eref' + label1 + '.png', format='png')
    
    
#%%