
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

save_graphs = 1 #guardar los graficos 2D del campo

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/potential_field/many_potential','')
path_save = path_basic + '/' + 'disp_relation'

if save_graphs==1:    
    if not os.path.exists(path_save):
        print('Creating folder to save graphs')
        os.mkdir(path_save)

try:
    sys.path.insert(1, path_constants)
    from graphene_sigma import sigma_DL
except ModuleNotFoundError:
    print('graphene_sigma.py no se encuentra en ' + path_basic)
    
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
cota = 1200 #nanometros
epsi1,epsi2 = 1,1
hbmu,hbgama = 0.3,0.0001
 


omega0THz = 25
omega0 = omega0THz*1e12
lambdda0 = c/(2*np.pi*omega0)
E0 =omega0THz*1e12*hb  #eV
omegac0 = E0/aux   

title = r'E = %.4feV' %(E0) + r', $\hbar\mu$ = %.2feV, $\hbar\gamma$ = %.4feV' %(hbmu,hbgama) 
labely = r'Im($r_p$)'
labelx = r'$k_\parallel$ [1/$\mu$m]'
k1 = omegac0*np.sqrt(epsi1*mu1)

list_k_parallel = np.linspace(0,k1,199)
     
def rp_imag(k_parallel):
    
    cond = alfac*sigma_DL(E0,hbmu,hbgama)*c
    den = 1 - 1j*omega0/(2*np.pi*k_parallel*cond)
    
    rta = 1/den
    return rta.imag
    
def recta(int_v,k_parallel):
    v = c/int_v
    kv = omega0/v
    y = k_parallel/kv
    return y 

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

list_int_v = [200,250,350]


listy = []
for k_parallel in list_k_parallel:
    rta = rp_imag(k_parallel)
    listy.append(rta) 
    
graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(list_k_parallel,listy,'.',ms = ms,color = 'purple')

for int_v in list_int_v: 
    rectay = recta(int_v,list_k_parallel)
    plt.plot(list_k_parallel,rectay,'-',ms = ms,label = 'v=c/%i $\mu$/s' %(int_v))

plt.legend(loc = 'best',markerscale=1,fontsize=tamlegend,frameon=False,handletextpad=0.05)
if save_graphs==1:
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'Im_rpE%.4f.png' %(E0), format='png')    
    
#%%
