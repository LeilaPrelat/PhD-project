
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: leila

graficar el campo externo reflejado en el plano con la convencion de z hacia abajo
en z = 0
"""

import numpy as np
import sys
import os 
import matplotlib.pyplot as plt
import seaborn as sns

#%%

save_graphs = 1 #guardar los graficos 2D del campo
plot_vs_omegaTHz = 1 # graficar vs omega THz
sns.set()    
Gradshtein_correcto = 1 # gradshtein es correcto, lo chequee con sympy : ver revisar_Gradshtein2.py

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/External_Efield/fieldE_ref/extra','')
if Gradshtein_correcto == 1:
    path_save = path_basic + '/' + 'fieldE_ref/formula_Gradshtein'
else:
    path_save = path_basic + '/' + 'fieldE_ref/formula_diferente'
#print('Importar modulos necesarios para este codigo')

if save_graphs==1:    
    if not os.path.exists(path_save):
        print('Creating folder to save graphs')
        os.mkdir(path_save)

err = 'fieldE_ref_numerical.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from revisar_Gradshtein import revisarG_num, revisarG_ana
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

#omega = 0.7*1e12
cota = 2
epsi1,epsi2 = 1,1
hbmu,hbgama = 0.3,0.0001

list_OmegaTHz = np.linspace(0.001,0.201,21)
#list_OmegaTHz = np.linspace(0.001,0.7,26)
list_OmegaTHz = np.linspace(1,7,51)

list_E = np.array(list_OmegaTHz)*1e12*hb*1e3
label1 = '_vsOmega'
title1 = r'$\epsilon_1$ = %i, v = c/%i $\mu$m/s' %(epsi1,int_v) 
labelx = r'E [meV]'
title2A = r'$\epsilon_2$ = %i, $\hbar\mu$ = %.2feV, $\hbar\gamma$ = %.4feV' %(epsi2,hbmu,hbgama) 
title = title1 + '\n'  + title2A   

#%%

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

listy_num_re = []
listy_num_im = []
    
listy_ana_re = []
listy_ana_im = []
    
j = 0 
for omegaTHz in list_OmegaTHz:
    omegaTHz = np.round(omegaTHz,9)
    omegac = omegaTHz*aux2  

                
    rta = revisarG_ana(omegac,hbmu,hbgama,int_v,Gradshtein_correcto)
    listy_ana_re.append(rta.real)
    listy_ana_im.append(rta.imag)
    
    rta = revisarG_num(omegac,hbmu,hbgama,int_v)
    listy_num_re.append(rta.real)
    listy_num_im.append(rta.imag)

    j = j + 1
    print(j)
    
#%%

graph(title,labelx,r'Re$(E_{ref,x})$',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(list_E,listy_ana_re,'.',ms = ms,color = 'blue',label = 'analytical')
plt.plot(list_E,listy_num_re,'-',ms = ms,color = 'magenta',label = 'numeric')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
#plt.xlim([1,4.7])
#plt.ylim([-0.01,0.01])
plt.yscale('log')
if save_graphs==1:
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'Re_Gradshtein_' + label1 + '.png', format='png')
    

graph(title,labelx,r'Im$(E_{ref,x})$',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(list_E,listy_ana_im,'.',ms = ms,color = 'blue',label = 'analytical')
plt.plot(list_E,listy_num_im,'-',ms = ms,color = 'magenta',label = 'numeric')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
plt.yscale('log')
if save_graphs==1:
    plt.tight_layout()
    plt.savefig( 'Im_Gradshtein_' + label1 + '.png', format='png')      
                
#%%
    
