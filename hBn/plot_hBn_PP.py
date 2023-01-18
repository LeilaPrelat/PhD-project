#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
conductividad del grafeno  
ref:
@BOOK{kubo2,
   author       = {R. A. Depine}, 
   year         = {2017},
   title        = {Graphene Optics: Electromagnetic solution of canonical problems}, 
   publisher    = {IOP Concise Physics.\ San Rafael, CA, USA: Morgan and Claypool Publishers}
}

"""

import numpy as np
import matplotlib.pyplot as plt
import os 
import sys

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')


#%%

try:
    sys.path.insert(1, path_basic)
    from hBn_PP import epsilon_x, epsilon_z
except ModuleNotFoundError:
    print('hBn_PP.py no se encuentra en ' + path_basic)


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
N = 400
listx = np.linspace(0.05, 0.25,N)

listy_re_z = []
listy_re_x = []

listy_im_z = []
listy_im_x = []

list_x_polaritons_coupling = []
list_ind = []
k = 0
for value in listx: 

#    y_re_ana = function_ana(value)       
    y_z = epsilon_z(value)  
    
    y_x = epsilon_x(value)

    
#    listy_re_ana.append(y_re_ana)
    
    listy_re_z.append(np.real(y_z))
    listy_re_x.append(np.real(y_x))

    listy_im_z.append(np.imag(y_z))
    listy_im_x.append(np.imag(y_x))


    if np.real(y_z)*np.real(y_x) < 0:
        
        list_x_polaritons_coupling.append(value)
        list_ind.append(k)
    k = k + 1 
    
#%%
list_ind2 = [list_ind[0]]
for ind in range(len(list_ind)-1):
    
    if list_ind[ind + 1 ] - list_ind[ind ] > 1:
        list_ind2.append(list_ind[ind ]) 
        list_ind2.append(list_ind[ind + 1 ]) 
list_ind2.append(list_ind[-1])    
#Emax = listx[np.argmax(listy_re_num)]
#print(Emax)

labelx = r'$\hbar\omega$ [eV]'
labely1 = r'Re{$\epsilon$}'
labely2 = r'Im{$\epsilon$}'
title = 'hBn'

#%%
graph(title,labelx,labely1,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#    plt.plot(listx,listy_re_ana,'.',ms = ms,color = 'purple',label = 'analytical')
plt.plot(listx,listy_re_x,'.-',ms = ms,color = 'lightseagreen',label = r'Re{$\epsilon_x$}')
plt.plot(listx,listy_re_z,'.-',ms = 3,color = 'darkred',label = r'Re{$\epsilon_z$}')
ejey_aux = np.linspace(np.min([np.min(listy_re_z),np.min(listy_re_x)]), np.max([np.max(listy_re_z),np.max(listy_re_x)]) , 10)
for ind in list_ind2:
    plt.plot(listx[ind]*np.ones(10), ejey_aux,'--',color = 'grey' )
    print(listx[ind])
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
os.chdir(path_basic)
plt.savefig( 'Re_epsilon_hBn' + '.png', format='png')   


#%%
graph(title,labelx,labely2,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#    plt.plot(listx,listy_re_ana,'.',ms = ms,color = 'purple',label = 'analytical')
plt.plot(listx,listy_im_x,'.-',ms = ms,color = 'lightseagreen',label = r'Im{$\epsilon_x$}')
plt.plot(listx,listy_im_z,'.-',ms = 3,color = 'darkred',label = r'Im{$\epsilon_z$}')
ejey_aux = np.linspace(np.min([np.min(listy_im_x),np.min(listy_im_z)]), np.max([np.max(listy_im_x),np.max(listy_im_z)]) , 10)
for ind in list_ind2:
    plt.plot(listx[ind]*np.ones(10), ejey_aux,'--',color = 'grey' )
    print(listx[ind])
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
os.chdir(path_basic)
plt.savefig( 'Im_epsilon_hBn' + '.png', format='png')   



#%%