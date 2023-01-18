
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
import matplotlib.colors as colors

#%%

save_graphs = 1 #guardar los graficos 2D del campo

plot_3D = 1

plot_vs_Energy = 0

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('potential_field/potential_and_electric_field/potential_and_electric_field_simple_version/potential_1_dipolo_analytical_vs_num','')
path_save = path_basic + '/' + 'plots_rp_fresnel_vs_pole_aprox'

if save_graphs==1:    
    if not os.path.exists(path_save):
        print('Creating folder to save graphs')
        os.mkdir(path_save)

err = 'fieldE_direct_numerical.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from rp_fresnel_vs_pole_aprox import compare_fresnel_vs_pole_approx
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

#omega = 0.7*1e12
cota = 75 #nanometros
epsi1,epsi2 = 1,1
hbmu,hbgama = 0.3,0.0001

title = r'$\hbar\mu$ = %.2feV, $\hbar\gamma$ = %.4feV' %(hbmu,hbgama) 

listE = np.linspace(0.01,90,100)
listQ = np.linspace(0.01,90,100)

# z0 = 0.06*1e3
labelx = 'E [meV]'   
labely = r'$\alpha_\parallel$'
label1 = 'vs_E' 

if plot_3D == 1 :    #omega = 0.7*1e12 

    def function_abs(energy0_meV,u):
        
        omegac0 = energy0_meV*1e-3/aux 
        
        rta = compare_fresnel_vs_pole_approx(omegac0,epsi1,epsi2,hbmu,hbgama,u)
        return rta

else: 
    k_parallel = 0.1
    title = title + r', $k_\parallel$ = %.2f 1/$\mu$m' %(k_parallel)
    
    
    def function_abs(energy0_meV):
        
        omegac0 = energy0_meV*1e-3/aux 
        u = k_parallel/omegac0
        
        rta = compare_fresnel_vs_pole_approx(omegac0,epsi1,epsi2,hbmu,hbgama,u)
        return rta
    

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

#%%

if plot_3D == 1:

    list_x = listE
    list_y = listQ
    X, Y = np.meshgrid(list_x, list_y)
    limits = [min(list_x) , max(list_x), min(list_y) , max(list_y)]
    
    f1_real = np.vectorize(function_abs)
    Z1_real = f1_real(X, Y)
   
    print('Graficar el External field' + label1)
    
    graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)    
    im = plt.imshow(Z1_real, extent = limits, cmap=plt.cm.hot, interpolation='bilinear')
    cbar = plt.colorbar(im)
    cbar.ax.tick_params(labelsize = tamnum)
    cbar.set_label(r'|1-$r^{PP}_p/r_p$|',fontsize=tamlegend,labelpad = 1)
    

    if save_graphs == True:
        os.chdir(path_save)
        plt.tight_layout()
        plt.savefig('rp_fresnel_pp_3D' + label1 + '.png', format='png')     
        
else:
    listy = function_abs(listE)

    graph(title,labelx,r'|1-$r^{PP}_p/r_p$|',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    plt.plot(listE,listy,'.',ms = ms,color = 'lightseagreen')
    # plt.yscale('log')
    plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
    if save_graphs==1:
        plt.tight_layout()
        os.chdir(path_save)
        plt.savefig( 'rp_fresnel_pp' + label1 + '.png', format='png')    

#%%
