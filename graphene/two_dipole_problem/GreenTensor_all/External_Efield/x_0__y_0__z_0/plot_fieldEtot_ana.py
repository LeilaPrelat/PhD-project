
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: leila

graficar campo total approximacion analitica
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

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/External_Efield/x_0__y_0__z_0','')
path_save = path_basic + '/' + 'fieldE_tot'

if save_graphs==1:    
    if not os.path.exists(path_save):
        print('Creating folder to save graphs')
        os.mkdir(path_save)

err = 'fieldE_ref_numerical.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from fieldE_ref import Efield_ana1b
except ModuleNotFoundError:
    print(err)
err2 = 'fieldE_direct_analyticalQE.py no se encuentra en ' + path_basic
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

int_v = 400
v = c/int_v
b = -0.01 #electron position in z < 0 (arriba del plano)
#omega = 0.7*1e12
cota = 2
epsi1,epsi2 = 1,1
hbmu,hbgama = 0.3,0.0001
zp = 0.05
z0 = 0 # z tiene que ser 0 para que la parte directa tenga solucion analitica

list_OmegaTHz = np.linspace(0.01,2.01,51)
    
list_OmegaTHz = np.linspace(10,810,201)
list_E = np.array(list_OmegaTHz)*1e12*hb
label1 = '_b%i_vsOmega' %(b)
title1 = r'$\epsilon_1$ = %i, v = c/%i $\mu$m/s, b = %inm' %(epsi1,int_v,b*1e3) 
labelx = r'E [eV]'
title2A = r'$\epsilon_2$ = %i, $\hbar\mu$ = %.2feV, $\hbar\gamma$ = %.4feV, $z_p$=%inm, z = %inm' %(epsi2,hbmu,hbgama,zp*1e3,z0*1e3) 
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

list_v_re = []
list_v_im = []
j=0

for omegaTHz in list_OmegaTHz: 
    
    omegaTHz = np.round(omegaTHz,8)
    omegac = omegaTHz*aux2
   
    v1 = Efield_ana1b(omegac,epsi1,epsi2,hbmu,hbgama,int_v,b,zp,z0)
    v2 = Efield_ANA_QE_2terms(omegac,epsi1,int_v,b)

    v = v1 + v2
    list_v_re.append(v.real)
    list_v_im.append(v.imag)

    print(j)
    j = j + 1


#%%
   
print('Graficar el External field tensor')
    
graph(title,labelx,r'Re$(E_{x})$',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(list_E,list_v_re,'.',ms = ms,color = 'lightseagreen')# plt.yscale('log')
#plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
if save_graphs==1:
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'ReEtot' + label1 + '.png', format='png')

graph(title,labelx,r'Im$(E_{x})$',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(list_E,list_v_im,'.',ms = ms,color = 'lightseagreen')
# plt.yscale('log')
#plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
if save_graphs==1:
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'ImEtot' + label1+ '.png', format='png')

#%%