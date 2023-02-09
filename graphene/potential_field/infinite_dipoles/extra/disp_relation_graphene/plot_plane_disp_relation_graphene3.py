#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 19 20:56:58 2020

@author: leila

relacion de dispersion
solucion analitica
para un plano de Ag

"""
import os 
import sys    
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_graphene =  path_basic.replace('/' + 'plane','')
#print('Importar modulos necesarios para este codigo')

try:
    sys.path.insert(1, path_graphene)
    from graphene_sigma import sigma_DL
except ModuleNotFoundError:
    print('graphene_sigma.py no se encuentra en ' + path_basic)
    
try:
    sys.path.insert(1, path_graphene)
    from constants import constantes
except ModuleNotFoundError:
    print('constants_plane.py no se encuentra en ' + path_basic)

pi,hb,c,alfac,mu1,mu2 = constantes()
aux = c*hb

#%%


def theta(int_v,a_nm,Nmax,omega_n):
    

    omegac = omega_n/c    
    
 #   omegac = omegac - 0.5 ## ver si funciona con una freq menor 

    E = omegac*aux    
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    alfa_p = 1j*(epsi1 + epsi2)/(cond)
    kp = alfa_p*omegac    
    
    lambdda_p = 2*pi/kp
    
    a = a_nm*1e-3
    
    theta0 = np.arccos((omegac*int_v + 2*np.pi*Nmax/a)/np.real(kp))
    
    return theta0

def disp_relation_omega(kx,a,Nmax,epsi1,epsi2,hbmu):

    """
    Parameters
    ----------
    omega : frecuencia en Hz
    mu_c : chemical potential of graphene in eV
    gamma_in : frecuencia de colision en eV
    epsilon1 : permeabilidad electrica del medio de arriba
    epsilon2 : permeabilidad electrica del medio de abajo
    Returns
        relacion de dispersion para un plano de Ag
        (solucion analitica)
    -------
    """

    cte = (epsi1 + epsi2)
    
    kn = kx + 2*np.pi*Nmax/a
    
    
    return np.sqrt(kn*4*hbmu*alfac*(c/hb)/cte)


def omega_n_corte(int_v,a_nm,Nmax):  ## omega puede valer esto o menos para que se generen plasmones
    
    a = a_nm*1e-3
    
    aux = alfac*int_v*hbmu/(hb)
    rta = aux + 0.5*np.sqrt((2*aux)**2 + 16*alfac*c*hbmu*Nmax*np.pi/(a*hb))
#    
#    rta2 = aux - 0.5*np.sqrt((2*aux)**2 + 16*alfac*c*hbmu*Nmax*np.pi/(a*hb))
    
#    print(rta*hb*1e3)
    
    return rta


def omega_min(int_v,n,a):
    
    A_num = np.pi*(epsi1 + epsi2)*hb
    A_den = hbmu*alfac*c
    A = A_num/A_den
        
    
    omega1 = (2*A*v)**(-1) + 0.5*np.sqrt((2*A*v)**(-2)  + 8*np.pi*n/(A*a) )
    omega2 = (2*A*v)**(-1) - 0.5*np.sqrt((2*A*v)**(-2)  + 8*np.pi*n/(A*a) )
    if n > 0:
        
        omega_min_val = omega1 
    
    else:
        
        omega_min_val = np.max([omega1,omega2])
    
    return omega_min_val
    
def omega_max(int_v,n,a):
    A_num = np.pi*(epsi1 + epsi2)*hb
    A_den = hbmu*alfac*c
    A = A_num/A_den
        
    
    omega1 = (2*A*v)**(-1) + 0.5*np.sqrt((2*A*v)**(-2)  + 8*np.pi*n/(A*a) )
    omega2 = (2*A*v)**(-1) - 0.5*np.sqrt((2*A*v)**(-2)  + 8*np.pi*n/(A*a) )
    if n < 0:
        
        omega_max_val = np.min([omega1,omega2])
    
        return omega_max_val
    

#%%

sns.set()
path_save = path_basic + '/' + 'disp_relation_graphene'

tamfig = (4.5,3.5)
tamlegend = 10
tamletra = 10
tamtitle = 10
tamnum = 9
loc2 = [0,1]
pad = -2
lw = 1.5
hp = 0.3
mk = 2

#%%

epsi1, epsi2 = 1, 1
hbmu = 0.3
hbgama = 0.0001

a = 0.05
a_nm = a*1e3
Nmax = -1

intt_v = 0.1

list_colours = ['black', 'gold', 'orange', 'darksalmon', 'brown', 'darkred']
title1 = '$\epsilon_1 = %.1f$, $\epsilon_2 = %.1f$, $\gamma_{in}$ = %.2f meV' %(epsi1,epsi2,hbgama*1e3)
title2 = '$\mu_c$ = %.2f eV, a = %i nm, n = %i' %(hbmu,a*1e3,Nmax)

def n_min(int_v,a):
    A_num = np.pi*(epsi1 + epsi2)*hb
    A_den = hbmu*alfac*c
    A = A_num/A_den
    v_2 = (c/int_v)**2
    n_min = -a/(8*np.pi*A*v_2)
    return n_min

if Nmax < n_min(intt_v,a):
    raise TypeError('Wrong value for N')


Ntot = 200
list_kx = np.linspace(-10,10, Ntot) ##creo que va a estar en 1/micrones

list_y_re = []
list_y_im = []
for kx in list_kx:
    
    valuey = disp_relation_omega(kx,a,Nmax,epsi1,epsi2,hbmu)
    valuey = hb*valuey*1e3
    list_y_re.append(valuey.real)
    list_y_im.append(valuey.imag)
    

plt.figure(figsize=tamfig)    
plt.title(title1 + '\n' + title2 ,fontsize=tamtitle)
plt.plot(list_kx,list_y_re,'-',color = 'darkgreen',ms = lw, label = 'DR')
list_hb_omega_meV = np.linspace(np.min(list_y_re),np.max(list_y_re),Ntot)
j = 0
list_cortes_omega = []
list_cortes_k = []


v = c*intt_v

list_y2_re = []

for hb_omega_meV in list_hb_omega_meV:
    omega = hb_omega_meV*1e-3/hb
    
    valuey_e = omega/v
    list_y2_re.append(valuey_e.real)
    

hb_omega_corte = omega_n_corte(1/intt_v,a_nm,Nmax)*hb*1e3

k_corte = omega_n_corte(1/intt_v,a_nm,Nmax)/v

list_cortes_omega.append(hb_omega_corte)
list_cortes_k.append(k_corte)

print(theta(intt_v,a_nm,Nmax,list_cortes_omega[0]*1e-3/(c*hb)))
##  
if intt_v != 1:
    plt.plot(list_y2_re,list_hb_omega_meV,'-', color = list_colours[0],ms = lw, label = 'v= %.1fc' %(intt_v))
    plt.plot(list_cortes_k,list_cortes_omega,'.', color = 'red',ms = '3')
else:
    plt.plot(list_y2_re,list_hb_omega_meV,'-', color = list_colours[0], ms = lw, label = 'light')

listy_1 = np.ones(Ntot)*omega_min(intt_v,Nmax,a)*hb*1e3
if Nmax < 0:
    listy_2 = np.ones(Ntot)*omega_max(intt_v,Nmax,a)*hb*1e3
plt.plot(list_kx,listy_1,'--', color = 'black', ms = lw, label = '$\hbar\omega_{min}$ = %.2f meV' %(omega_min(intt_v,Nmax,a)*hb*1e3))
if Nmax < 0:
    plt.plot(list_kx,listy_2,'--', color = 'black', ms = lw, label = '$\hbar\omega_{max}$ = %.2f meV'%(omega_max(intt_v,Nmax,a)*hb*1e3))
plt.xlabel('$k_\parallel$ [1/$\mu$m]',fontsize=tamletra, labelpad = 0)
plt.ylabel('$\hbar\omega$ [meV]',fontsize=tamletra, labelpad = 0)
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
plt.tick_params(labelsize = tamnum,pad = pad)
plt.grid(1)
plt.tight_layout()
plt.xlim([0,10])
# os.chdir(path_save)
plt.savefig('disp_relation_graphene_real_vs_mEv_n%i.png'%(Nmax))

#%%


