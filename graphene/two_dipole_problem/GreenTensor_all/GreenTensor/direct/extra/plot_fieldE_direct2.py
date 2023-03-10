
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: leila

graficar green tensor Gnn con la misma convencion 
que el paper 370 (z hacia abajo, zplane>0, particula en r = 0, medio 1 arriba del plano, 
                  medio 2 abajo del plano)
"""

import numpy as np
import sys
import os 
import matplotlib.pyplot as plt
import seaborn as sns


#%%

save_graphs = 1 #guardar los graficos 2D del campo
plot_sin_QE = 0
plot_num_QE = 1
sns.set()    

name_this_py = os.path.basename(__file__)

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('GreenTensor/' + 'direct','')
path_save = path_basic + '/' + 'fieldE_direct'

if save_graphs==1:    
    if not os.path.exists(path_save):
        print('Creating folder to save graphs')
        os.mkdir(path_save)

err = 'fieldE_direct_numerical.py no se encuentra en ' + path_basic
err2 = 'fieldE_direct_analyticalQE.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from fieldE_direct_numerical2 import Efield_NUM_QE_2terms,Efield_NUM_2terms 
except ModuleNotFoundError:
    print(err)

try:
    sys.path.insert(1, path_basic)
    from fieldE_direct_analyticalQE import Efield_ANA_QE_2terms
except ModuleNotFoundError:
    print(err2)

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
#v = c/den
b = -20 #electron position in z < 0 (arriba del plano)
#omega = 0.7*1e12
omegaTHz = 0.7
omegac = omegaTHz*aux2 
list_OmegaTHz = np.linspace(0.01,2.01,251)
list_OmegaTHz = np.linspace(0.01,0.08,36)
#list_b = np.linspace(-20,-200,100)
epsi1 = 1

#%%
    
labelx = '$\omega$ [THz]'

tamfig = (4.5,3.5)
tamlegend = 10
tamletra = 11
tamtitle = 10
tamnum = 9
labelpady = -1.5
labelpadx = 0
pad = 0.5
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

color_rgb1 = (1.0, 0.8, 0.6)
color_rgb2 = (0.99, 0.74, 0.71)
color_rgb3 = (0.529, 0.992, 0.945)
color_rgb4 = (0.41, 0.21, 0.61)

#%%

title1 = r'$\epsilon_1$ = %i, v = c/%i $\mu$m/s, b = %i$\mu$m' %(epsi1,int_v,b)
<<<<<<< HEAD
title = title1 + ', ' + name_this_py?? 
=======
title = title1 + ', ' + name_this_py
>>>>>>> 9d609f242713b6f4e8a1bb0984bdd69647f4be64

list_v1_re = []
list_v1_im = []

list_v2_re = []
list_v2_im = []

list_v3_re = []
list_v3_im = []

j=0

for omegaTHz in list_OmegaTHz: 
    
    omegaTHz = np.round(omegaTHz,8)
    
    omegac = omegaTHz*aux2
    v1_05,v1_06 = Efield_ANA_QE_2terms(omegac,epsi1,int_v,b)
    v1 = v1_05 + v1_06
    
    list_v1_re.append(v1.real)
    list_v1_im.append(v1.imag)

    if plot_num_QE == 1:    
        v2 = Efield_NUM_QE_2terms(omegac,epsi1,int_v,b)

        list_v2_re.append(v2.real)
        list_v2_im.append(v2.imag)
    
    if plot_sin_QE == 1:
        v3 = Efield_NUM_2terms(omegac,epsi1,int_v,b)

        list_v3_re.append(v3.real)
        list_v3_im.append(v3.imag)

    print(j)
    j = j + 1
    
#%%

#list_OmegaTHz = list_OmegaTHz[0:j+1]

#%%

if plot_num_QE == 1 and plot_sin_QE == 0:
    header1 = 'omega [THz]     Re(E_{dir}) analytic     Re(E_{dir}) QE num' + ', ' +  title
    tabla1 = np.array([list_OmegaTHz,list_v1_re,list_v2_re])
elif plot_num_QE == 1 and plot_sin_QE == 1:
    header1 = 'omega [THz]     Re(E_{dir}) analytic     Re(E_{dir}) QE num     Re(E_{dir}) sin QE num' + ', ' +  title
    tabla1 = np.array([list_OmegaTHz,list_v1_re,list_v2_re,list_v3_re])   
tabla1 = np.transpose(tabla1)
 
#####################################

if plot_num_QE == 1 and plot_sin_QE == 0:
    header2 = 'omega [THz]     Im(E_{dir}) analytic     Im(E_{dir}) QE num' + ', ' +  title
    tabla2 = np.array([list_OmegaTHz,list_v1_im,list_v2_im])
elif plot_num_QE == 1 and plot_sin_QE == 1:
    header2 = 'omega [THz]     Im(E_{dir}) analytic     Im(E_{dir}) QE num     Im(E_{dir}) sin QE num' + ', ' +  title
    tabla2 = np.array([list_OmegaTHz,list_v1_im,list_v2_im,list_v3_im])
tabla2 = np.transpose(tabla2)

#%%
   
print('Graficar el External field tensor')
    
graph(title,labelx,r'Re$(E_{dir})$',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(list_OmegaTHz,list_v1_re,'-',color = 'black',label = 'analytical QE')
if plot_num_QE == 1:
    plt.plot(list_OmegaTHz,list_v2_re,'.',ms = ms,color = color_rgb1,label = 'numerical QE')
if plot_sin_QE == 1:
    plt.plot(list_OmegaTHz,list_v3_re,'--',color = 'green',label = 'numerical')
# plt.yscale('log')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
if save_graphs==1:
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'Re_Edir' + '.png', format='png')
    np.savetxt('Re_Edir' + '.txt', tabla1, fmt='%1.11e', delimiter='\t', header = header1)

graph(title,labelx,r'Im$(E_{dir})$',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(list_OmegaTHz,list_v1_im,'-',color = 'black',label = 'analytical QE')
if plot_num_QE == 1:
    plt.plot(list_OmegaTHz,list_v2_im,'.',ms = ms,color = color_rgb1,label = 'numerical QE')
if plot_sin_QE == 1:
    plt.plot(list_OmegaTHz,list_v3_im,'--',color = 'green',label = 'numerical')
# plt.yscale('log')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
if save_graphs==1:
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'Im_Edir' + '.png', format='png')
    np.savetxt('Im_Edir' + '.txt', tabla2, fmt='%1.11e', delimiter='\t', header = header2)

#%%