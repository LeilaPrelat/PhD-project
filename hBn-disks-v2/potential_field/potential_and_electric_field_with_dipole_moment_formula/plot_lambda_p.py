
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
from scipy.signal import find_peaks
import seaborn as sns
sns.set()

#%%

plot_vs_E = 0
plot_vs_c = 0
plot_vs_zp = 1

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/potential_field//potential_and_electric_field_with_dipole_moment_formula','')
path_save = path_basic + '/' + 'decay_rate_film'
if not os.path.exists(path_save):
    print('Creating folder to save graphs')
    os.mkdir(path_save)

err = 'decay_rate_film3.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from decay_rate_film_resonance import EELS_film_ana_f, EELS_film_num_f, EELS_film_pole_aprox_f
except ModuleNotFoundError:
    print(err)



try:
    sys.path.insert(1, path_constants)
    from hBn_PP import epsilon_x
except ModuleNotFoundError:
    print('epsilon_z.py no se encuentra en ' + path_constants)


try:
    sys.path.insert(1, path_basic)
    from Silica_epsilon import epsilon_Silica
except ModuleNotFoundError:
    print('Silica_epsilon.py no se encuentra en ' + path_basic)
    

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


x1 = 0.09260651629072682 
x2 = 0.10112781954887218
x3 = 0.17030075187969923
x4 = 0.19937343358395992

x0 = 0.08684904004367106 
xf =  0.8064559363534132

#x0 = x1 + 1e-3
#xf = x4 - 1e-3

b = -0.01

d_nano = 0.4

def lambda_p(energy0):
    
    epsi_x = epsilon_x(energy0)
    epsi_HBN_par = epsi_x
    epsi_silica = epsilon_Silica(energy0)        
    
    omegac = energy0/aux
#    d_micros = d_nano*1e-3
    d_micro = d_nano*1e-3
    alfa_p = epsi_silica*2/(omegac*d_micro*(epsi_HBN_par-1))
    kp = alfa_p*omegac

    return (2*np.pi/kp)*1e3 ## en micro


N = 1e3
listx = np.linspace(x0, xf, N) 
#listx = np.linspace(x0, 0.8, N) 

listy_re = []
listy_im = []
list_re_lambda_pos = []
list_ind = []

k = 0
for x in listx : 
    y = lambda_p(x)
    listy_re.append(np.real(y))
    
    if np.real(y) > 0 : 
        list_re_lambda_pos.append(x)
        list_ind.append(k)
        
    listy_im.append(np.imag(y))
    k = k + 1


list_ind2 = [list_ind[0]]
for ind in range(len(list_ind)-1):
    
    if list_ind[ind + 1 ] - list_ind[ind ] > 1:
        list_ind2.append(list_ind[ind ]) 
        list_ind2.append(list_ind[ind + 1 ]) 
list_ind2.append(list_ind[-1])  

title = r'b = %i nm, d = %.2f nm' %(b*1e3,d_nano)
labelp = r'_d%inm' %(d_nano) 

N = 30
labelx = r'$\hbar\omega$ [eV]'

#%%
    
tamfig = (4.5,3.5)
tamlegend = 12
tamletra = 11
tamtitle = 10
tamnum = 9
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

ejey_aux = np.linspace(np.min(listy_re), np.max(listy_re) , 10)

graph(title,labelx, r'Re($\lambda_{\rm p}$) [nm]',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(listx,listy_re,'.-',ms = ms,color = 'purple', label = 'PP analytical')
#plt.plot(np.ones(10)*maxi, listy_aux,'-k')
for ind in list_ind2:
    plt.plot(listx[ind]*np.ones(10), ejey_aux,'--',color = 'grey' )
    print(listx[ind])
plt.legend(loc = 'best',markerscale=1.5,fontsize=tamlegend,frameon=0.1,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
#if plot_vs_c == 1:
#    plt.yscale('log')
os.chdir(path_save)
plt.savefig( 'lambda_p_re' + labelp + '.png', format='png')   

graph(title,labelx, r'Im($\lambda_{\rm p}$) [nm]',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(listx,listy_im,'.-',ms = ms,color = 'purple', label = 'PP analytical')
#plt.plot(np.ones(10)*maxi, listy_aux,'-k')
plt.legend(loc = 'best',markerscale=1.5,fontsize=tamlegend,frameon=0.1,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
#if plot_vs_c == 1:
#    plt.yscale('log')
os.chdir(path_save)
plt.savefig( 'lambda_p_im' + labelp + '.png', format='png')   


#%%


