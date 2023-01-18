
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
err = 'fieldE_direct_z.py no se encuentra en ' + path_basic

try:
    sys.path.insert(1, path_basic)
    from fieldE_direct_z import Efield_NUM_QE,Efield_ANA
except ModuleNotFoundError:
    print(err)


path_save = path_basic + '/' + 'fieldE_direct_z'

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
#v = c/den
#omega = 0.7*1e12

N = 100
epsi1 = 1



z = 0.025

title = r'z = %i nm, v = c/%i' %(z*1e3,int_v)

########################
if plot_vs_E == 1:
#    omegaTHz0 = 0.7
#    omegaTHz0 = 1.511
    

    # omegaTHz = 1.663 - 0.152j
    
    R0 = 1000
    listx = np.linspace(15,65,N)

    labelx = '$\hbar\omega$ [meV]'
    label1 = '_vs_E'
    title = title +  r', R = %i nm' %(R0)


elif plot_vs_x == 1:
#    omegaTHz0 = 0.7
#    omegaTHz0 = 1.511
    

    # omegaTHz = 1.663 - 0.152j
    
    E0 = 43
    listx = np.linspace(50,1000,N)

    labelx = 'R [nm]'
    label1 = '_vs_R'
    title = title +  r', $\hbar\omega$ = %i meV' %(E0)
#title = title + '\n' + name_this_py 


def Efield_ANA_QE_function(energy0,R_nano):
    
    omegac0 = energy0*1e-3/aux 
    R = R_nano*1e-3
    
    return Efield_ANA(omegac0,epsi1,R,z,int_v)

def Efield_NUM_QE_function(energy0,R_nano):
    
    omegac0 = energy0*1e-3/aux 
    R = R_nano*1e-3
    
    return Efield_NUM_QE(omegac0,epsi1,R,z,int_v)


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

if plot_vs_x == 1:
    for value in listx: 
        
        tot2 = Efield_ANA_QE_function(E0,value)
        total2_re.append(tot2.real)
        total2_im.append(tot2.imag)
    

            
        tot1 = Efield_NUM_QE_function(E0,value)
   
        total1_re.append(tot1.real)
        total1_im.append(tot1.imag)
    

if plot_vs_E == 1:
    for value in listx:
        
        tot2 = Efield_ANA_QE_function(value,R0)
        total2_re.append(tot2.real)
        total2_im.append(tot2.imag)
    

            
        tot1 = Efield_NUM_QE_function(value,R0)
   
        total1_re.append(tot1.real)
        total1_im.append(tot1.imag)

    
    
#%%


graph(title,labelx,r'Re($E_{dir}$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(listx,total2_re,'-',color = 'darkred',label = 'analytical')
# plt.yscale('log')
plt.plot(listx,total1_re,'.',ms = ms,color = 'green',label = 'numerical')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
if save_graphs==1:
    os.chdir(path_save)
    plt.tight_layout()
    plt.savefig( 'Re_Edir' + label1 + '.png', format='png')

graph(title,labelx,r'Im($E_{dir}$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(listx,total2_im,'-',color = 'darkred',label = 'analytical')
# plt.yscale('log')
plt.plot(listx,total1_im,'.',ms = ms,color = 'green',label = 'numerical')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
if save_graphs==1:
    plt.tight_layout()
    plt.savefig( 'Im_Edir' + label1 + '.png', format='png')
    
    
#%%
