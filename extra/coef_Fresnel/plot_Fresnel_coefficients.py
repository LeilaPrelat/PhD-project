#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: leila

graficar los coeficientes de Fresnel 

"""

import numpy as np
import sys
import os 
import matplotlib.pyplot as plt
import seaborn as sns

#%%

save_graphs = 1 #guardar los graficos 2D del campo
sns.set() 
   
#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/' + 'coef_Fresnel','')

err = 'coef_Fresnel_graphene.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from coef_Fresnel_graphene import coefTM_TE
except ModuleNotFoundError:
    print(err)

try:
    sys.path.insert(1, path_constants)
    from constants import constantes
except ModuleNotFoundError:
    print('constants.py no se encuentra en ' + path_constants)

pi,hb,c,alfac,mu1,mu2 = constantes()
aux = c*hb

#%%

print('Definir parametros del problema')

hbgamma = 0.0001  # collision frequency of graphene in eV
hbmu = 0.3        # chemical potential of graphene in eV
omegaTHz = 0.7
omegac = omegaTHz*1e12/c 
epsi1 = 1
epsi2 = 1
title = r'$\epsilon_1$ = %i, $\epsilon_2$ = %i, $\mu$ = %.2feV, $\gamma$ = %.4feV' %(epsi1,epsi2,hbmu,hbgamma)
 
list1_alpha_parallel = np.linspace(0,0.98,99)
list2_alpha_parallel = np.linspace(1.02,25.02,201)

def omega_esp(alpha_parallel):
    n1 = np.sqrt(epsi1*mu1)
    epsi_mean = epsi1 + epsi2
    term1 = alfac*hbmu*n1/epsi_mean
    
    return term1/hb - 1j*hbgamma/hb

#%%
    
labelx = r'$\alpha_\parallel$'

tamfig = (4.5,3.5)
tamlegend = 10
tamletra = 11
tamtitle = 10
tamnum = 9
labelpady = -1.5
labelpadx = 0
pad = 0.5
hp = 0.1
mk = 2

def graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad):
    plt.figure(figsize=tamfig)
    plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
    plt.ylabel(labely,fontsize=tamletra,labelpad =labelpadx)
    plt.tick_params(labelsize = tamnum, pad = pad)
    plt.title(title,fontsize=int(tamtitle*0.9))

    return   

#%%
tot_alpha_parallel = []

list_rs_re = []
list_rp_re = []
list_esp_rs_re = []
list_esp_rp_re = []

list_rs_im = []
list_rp_im = []
list_esp_rs_im = []
list_esp_rp_im = []

omegaTHz = 0.7
omegac0 = omegaTHz*1e12/c 

for alpha_parallel in list1_alpha_parallel: 
    coef_rs, coef_rp = coefTM_TE(omegac0,epsi1,epsi2,hbmu,hbgamma,alpha_parallel)
    list_rs_re.append(coef_rs.real)
    list_rp_re.append(coef_rp.real)

    list_rs_im.append(coef_rs.imag)
    list_rp_im.append(coef_rp.imag)

    omegac = omega_esp(alpha_parallel)/c
    coef_rs, coef_rp = coefTM_TE(omegac,epsi1,epsi2,hbmu,hbgamma,alpha_parallel)
    list_esp_rs_re.append(coef_rs.real)
    list_esp_rp_re.append(coef_rp.real)
    
    list_esp_rs_im.append(coef_rs.imag)
    list_esp_rp_im.append(coef_rp.imag)
    
    tot_alpha_parallel.append(alpha_parallel)

for alpha_parallel in list2_alpha_parallel: 
    coef_rs, coef_rp = coefTM_TE(omegac0,epsi1,epsi2,hbmu,hbgamma,alpha_parallel)

    list_rs_re.append(coef_rs.real)
    list_rp_re.append(coef_rp.real)

    list_rs_im.append(coef_rs.imag)
    list_rp_im.append(coef_rp.imag)

    omegac = omega_esp(alpha_parallel)/c    
    coef_rs, coef_rp = coefTM_TE(omegac,epsi1,epsi2,hbmu,hbgamma,alpha_parallel)
    list_esp_rs_re.append(coef_rs.real)
    list_esp_rp_re.append(coef_rp.real)    

    list_esp_rs_im.append(coef_rs.imag)
    list_esp_rp_im.append(coef_rp.imag)    

    tot_alpha_parallel.append(alpha_parallel)


#%%

# ver el maximo del coeficiente de fresnel rp 

maxi = np.argmax(list_esp_rp_re)
alpha_maxi = tot_alpha_parallel[maxi]
omegaTHz_maxi_re = (omega_esp(alpha_maxi).real)*1e-12
omegaTHz_maxi_im = (omega_esp(alpha_maxi).imag)*1e-12

#%%
print('Graficar los coeficientes de Fresnel')
    
graph(title,labelx,r'Re $r_s$',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(tot_alpha_parallel,list_rs_re,'.',ms = 4,color = 'green',label = 'pol s')
plt.plot(tot_alpha_parallel,list_esp_rs_re,'.',ms = 4,color = 'black',label = r'pol s $\omega_p$')
# plt.yscale('log')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp)
if save_graphs==1:
    plt.tight_layout()
    plt.savefig( 'Recoef_rs' + '.png', format='png')

graph(title,labelx,r'Im $r_s$',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(tot_alpha_parallel,list_rs_im,'.',ms = 4,color = 'blue',label = 'pol s')
plt.plot(tot_alpha_parallel,list_esp_rs_im,'.',ms = 4,color = 'black',label = r'pol s $\omega_p$')
# plt.yscale('log')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp)
if save_graphs==1:
    plt.tight_layout()
    plt.savefig( 'Imcoef_rs' + '.png', format='png')

#####################################################################################################

graph(title,labelx,r'Re $r_p$',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(tot_alpha_parallel,list_rp_re,'.',ms = 4,color = 'green',label = 'pol p')
plt.plot(tot_alpha_parallel,list_esp_rp_re,'.',ms = 4,color = 'black',label =  r'pol p $\omega_p$')
plt.plot([],[],'.',ms = 0,color = 'white', label = r'max Re$\omega$ = %.3f THz'%(omegaTHz_maxi_re))
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp)
if save_graphs==1:
    plt.tight_layout()
    plt.savefig( 'Recoef_rp' + '.png', format='png')


graph(title,labelx,r'Im $r_p$',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(tot_alpha_parallel,list_rp_im,'.',ms = 4,color = 'blue',label = 'pol p')
plt.plot(tot_alpha_parallel,list_esp_rp_im,'.',ms = 4,color = 'black',label =  r'pol p $\omega_p$')
plt.plot([],[],'.',ms = 0,color = 'white', label = r'max Im$\omega$ = %.3f THz'%(omegaTHz_maxi_im))
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp)
if save_graphs==1:
    plt.tight_layout()
    plt.savefig( 'Imcoef_rp' + '.png', format='png')

#%%
