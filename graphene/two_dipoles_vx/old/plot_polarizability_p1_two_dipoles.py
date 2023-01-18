#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: leila

graficar polarizabilidad
"""

import numpy as np
import sys
import os 
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.signal import find_peaks

sns.set()

#%%

save_graphs = 1 #guardar los graficos 2D del campo
    
#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/two_dipoles' ,'')
path_save = path_basic + '/' + 'p1_two_dipoles'

if save_graphs==1:    
    if not os.path.exists(path_save):
        print('Creating folder to save graphs')
        os.mkdir(path_save)

err = 'function_polarizability.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from function_polarizability import p1_2dip
except ModuleNotFoundError:
    print(err)
try:
    sys.path.insert(1, path_constants)
    from constants import constantes
except ModuleNotFoundError:
    print('constants.py no se encuentra en ' + path_constants)

pi,hb,c,alfac,mu1,mu2 = constantes()


#%%

print('Definir parametros del problema')

#omega = 0.7*1e12

epsi1 = 1
epsi2 = 1

b = -0.01
zp = 0.05

int_v = 10
hbmu,hbgama = 0.3,0.0001
zp = 0.05 #micrometros 
xD1,yD1,zD1 = -0.05,-0.05,0
xD2,yD2,zD2 = 0.05,0.05,0

omega0THz = 90
omega0 = omega0THz*1e12 
E0 = omega0*hb

kappa_factor_omega0 = 0.1
kappa_r_factor= 0.5

Nenergy = 300
list_meV = np.linspace(0.1,80,Nenergy) 

def f_alpha_function(omegac):
    return p1_2dip(omegac,epsi1,epsi2,hbmu,hbgama,int_v,zp,b,xD1,yD1,zD1,xD2,yD2,zD2,omega0,kappa_factor_omega0,kappa_r_factor)

#%%

title1 = r'$\kappa$ = %.2f$\omega_0$, $\kappa_r$ = %.2f$\kappa$, $E_0$ =%imeV' %(kappa_factor_omega0,kappa_r_factor,E0*1e3)    
title2 = r'v = c/%i $\mu$m/s, $\hbar\mu$ = %.2feV, $\hbar\gamma$ = %.4feV' %(int_v,hbmu,hbgama) 
title3 = r'$z_p$=%inm, b=%inm, xD1 = %inm, xD2 = %inm' %(zp*1e3, b*1e3, xD1*1e3,xD2*1e3)
title = title1 + '\n' + title2  + '\n' + title3

labelx = r'E [meV]'

tamfig = (4.5,3.5)
tamlegend = 10
tamletra = 11
tamtitle = 10
tamnum = 9
labelpady = 2
labelpadx = 2
pad = 0
length_marker = 1

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

list_alpha_re = []
list_alpha_im = []
list_alpha_tot = []

aux = c*hb
for E in list_meV: 

#    OmegaTHz = np.round(OmegaTHz,10)
    omegac = E*1e-3/aux

    rta1 = f_alpha_function(omegac)
    rta1_re = rta1.real
    rta1_im = rta1.imag
    
    list_alpha_re.append(rta1_re)
    list_alpha_im.append(rta1_im)
    list_alpha_tot.append(np.abs(rta1))    

#peaks_im, _ = find_peaks(list_alpha_im, height=0)
#maxis_x_im = []
#maxis_y_im = []
#
#for maxs in peaks_im: 
#    Emax = list_E[maxs]
#    omegac_max = Emax/aux    
#    maxis_x_im.append(Emax)
#    maxis_y_im.append(f_alpha_function(omegac_max).imag)
#
#maxi_tot_im = np.argmax(maxis_y_im)
#maxis_tot_x_im = maxis_x_im[maxi_tot_im]
#
#peaks_re, _ = find_peaks(-np.array(list_alpha_re), height=0)
##maxis_x_re = []
##maxis_y_re = []
#
#maxis_tot_x_re = list_E[int(peaks_re)]

#for maxs in peaks_re: 
#    Emax = list_E[maxs]
#    omegac_max = Emax/aux    
#    maxis_x_re.append(Emax)
#    maxis_y_re.append(f_alpha_function(omegac_max).real)
#
#maxi_tot_re = np.argmax(-np.array(maxis_y_re))
#maxis_tot_x_re = maxis_x_im[maxi_tot_re]
#plt.plot(peaks, list_y[peaks], "x")

#%%
print('Graficar el green tensor')
    
graph(title,labelx,r'Re($p_{1,x}$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(list_meV,list_alpha_re,'k-')
#plt.plot(np.ones(Nenergy)*maxis_tot_x_re, list_alpha_re,label = 'E = %.2feV'%(maxis_tot_x_re))
#plt.legend(loc = 'best',markerscale=1.5,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
if save_graphs==1:
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'Re_p1' + '.png', format='png')

graph(title,labelx,r'Im($p_{1,x}$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(list_meV,list_alpha_im,'k-')
#plt.plot(np.ones(Nenergy)*maxis_tot_x_im, list_alpha_im,label = 'E = %.2feV'%(maxis_tot_x_im))
#plt.legend(loc = 'best',markerscale=1.5,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
# plt.yscale('log')
if save_graphs==1:
    plt.tight_layout()
    plt.savefig( 'Im_p1' + '.png', format='png')


graph(title,labelx,r'|$p_{1,x}$|',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(list_meV,list_alpha_tot,'k-')
#plt.plot(np.ones(Nenergy)*maxis_tot_x_im, list_alpha_im,label = 'E = %.2feV'%(maxis_tot_x_im))
#plt.legend(loc = 'best',markerscale=1.5,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
# plt.yscale('log')
if save_graphs==1:
    plt.tight_layout()
    plt.savefig( '|p1|' + '.png', format='png')
    
#%%
