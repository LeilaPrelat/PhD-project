
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: leila

graficar el campo externo reflejado en el plano con la convencion de z hacia abajo
"""

import numpy as np
import sys
import os 
import matplotlib.pyplot as plt
import seaborn as sns

#%%

save_graphs = 1 #guardar los graficos 2D del campo
plot_sin_QE = 1
sns.set()    

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/External_Efield/x_0__y_0__z_0','')
path_save = path_basic + '/' + 'fieldE_ref'

if save_graphs==1:    
    if not os.path.exists(path_save):
        print('Creating folder to save graphs')
        os.mkdir(path_save)

err = 'fieldE_ref_numerical.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from fieldE_ref_numerical import Efield_NUM_QE,Efield_NUM 
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

int_v = 400
v = c/int_v
b = -2 #electron position in z < 0 (arriba del plano)
#omega = 0.7*1e12
omegaTHz = 0.7
omegac = omegaTHz*aux2 
list_OmegaTHz = np.linspace(0.01,2.01,251)
list_OmegaTHz = np.linspace(0.001,0.2,200)
epsi1,epsi2 = 1,1
hbmu,hbgama,zp = 0.7,0.0001,2

#%%
    
labelx = '$\omega$ [THz]'
label1 = '_b%i' %(b)

tamfig = (4.5,3.5)
tamlegend = 10
tamletra = 11
tamtitle = 10
tamnum = 9
labelpady = -1.5
labelpadx = 0.8
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

color_rgb1 = (1.0, 0.01, 0.24)
color_rgb2 = (0.8, 0.58, 0.46)
color_rgb3 = (0.55, 0.71, 0)
color_rgb4 = (0.6, 0.73, 0.89)

#%%

title1 = r'$\epsilon_1$ = %i, $\epsilon_2$ = %i, v = c/%i $\mu$m/s, b = %i$\mu$m' %(epsi1,epsi2,int_v,b)
title2 = r'$\hbar\mu$ = %.2feV, $\hbar\gamma$ = %.4feV, $z_p=%.2f\mu$m' %(hbmu,hbgama,zp)
title = title1 + '\n'  + title2 

list_v2_re = []
list_v2_im = []

list_v3_re = []
list_v3_im = []

j=0

for omegaTHz in list_OmegaTHz: 
    
    omegaTHz = np.round(omegaTHz,8)
    omegac = omegaTHz*aux2
   
    v2 = Efield_NUM_QE(omegac,epsi1,epsi2,hbmu,hbgama,int_v,b,zp)

    list_v2_re.append(v2.real)
    list_v2_im.append(v2.imag)

    if plot_sin_QE == 1:
        v3 = Efield_NUM(omegac,epsi1,epsi2,hbmu,hbgama,int_v,b,zp)

        list_v3_re.append(v3.real)
        list_v3_im.append(v3.imag)

    print(j)
    j = j + 1
    
#%%

#list_OmegaTHz = list_OmegaTHz[0:j+1]

#%%

os.chdir(path_save)
np.savetxt('list_OmegaTHz' + label1 + '.txt', list_OmegaTHz, fmt='%1.11e', delimiter='\t', header = title)
header1 = 'omega [THz]     Re(E_{ref,x}) QE num' + ', ' +  title
tabla1 = np.array([list_OmegaTHz,list_v2_re])
if plot_sin_QE == 1:
    header1 = 'omega [THz]     Re(E_{ref,x}) QE num     Re(E_{ref,x}) sin QE num' + ', ' +  title
    tabla1 = np.array([list_OmegaTHz,list_v2_re,list_v3_re])   
tabla1 = np.transpose(tabla1)
 
#####################################

header2 = 'omega [THz]     Im(E_{ref,x}) QE num' + ', ' +  title
tabla2 = np.array([list_OmegaTHz,list_v2_im])
if plot_sin_QE == 1:
    header2 = 'omega [THz]     Im(E_{ref,x}) QE num     Im(E_{ref,x}) sin QE num' + ', ' +  title
    tabla2 = np.array([list_OmegaTHz,list_v2_im,list_v3_im])
tabla2 = np.transpose(tabla2)

#%%
   
print('Graficar el External field tensor')
    
graph(title,labelx,r'Re$(E_{ref,x})$',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(list_OmegaTHz,list_v2_re,'.',ms = ms,color = 'lightseagreen',label = 'numerical QE')
if plot_sin_QE == 1:
    plt.plot(list_OmegaTHz,list_v3_re,'--',color = 'darkviolet',label = 'numerical')
# plt.yscale('log')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
if save_graphs==1:
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'ReE_ref' + label1 + '.png', format='png')
    np.savetxt('ReE_ref' + label1+ '.txt', tabla1, fmt='%1.11e', delimiter='\t', header = header1)

graph(title,labelx,r'Im$(E_{ref,x})$',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(list_OmegaTHz,list_v2_im,'.',ms = ms,color = 'lightseagreen',label = 'numerical QE')
if plot_sin_QE == 1:
    plt.plot(list_OmegaTHz,list_v3_im,'--',color = 'darkviolet',label = 'numerical')
# plt.yscale('log')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
if save_graphs==1:
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'ImE_ref' + label1+ '.png', format='png')
    np.savetxt('ImE_ref'+ label1 + '.txt', tabla2, fmt='%1.11e', delimiter='\t', header = header2)

#%%