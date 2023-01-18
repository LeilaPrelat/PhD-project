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
import math

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.interpolate import interp1d

sns.set()

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_ctes =  path_basic.replace('/' + 'finding_poles','')
#print('Importar modulos necesarios para este codigo')
path_save = path_basic + '/poles'
path_load = path_ctes
try:
    sys.path.insert(1, path_ctes)
    from constants import constantes
except ModuleNotFoundError:
    print('constants.py no se encuentra en ' + path_ctes)

pi,hb,c,alfac,mu1,mu2,mu3 = constantes()
aux = c*hb

save_plots = 1
save_txt = 1

 
epsilon2 = 12
d = 20*1e-3
omega_omegaWG = 2.5
title = '$\epsilon_2$ = %i, d = %inm' %(epsilon2,d*1e3)

if omega_omegaWG > 4 or omega_omegaWG < 5*1e-3:
    raise TypeError('Wrong value for omega/omegaWG')

omegaWG = (np.pi*c/d)/(np.sqrt(epsilon2-1)) 

#%%

tamfig = (4.5,3.5)
tamlegend = 10
tamletra = 12
tamtitle = 10
tamnum = 9
loc2 = [0,1]
pad = -2
lw = 1.5
hp = 0.3
mk = 2
labelpady = 0
labelpadx = 0
ms = 3


def k_parallel_air(omega):

    return omega/c

def k_parallel_medium(omega,epsilon2):

    return omega*np.sqrt(epsilon2)/c

#%%
 
def function_poles(type_sol,value):
    
    os.chdir(path_load)
    tabla_TM1A = np.loadtxt('sol_%s_d%inm.txt' %(type_sol,d*1e3), delimiter='\t', skiprows=1)
    tabla_TM1A = np.transpose(tabla_TM1A)
    [listx,listy] = tabla_TM1A
    
    if np.min(listx) <= value <= np.max(listx):
    
        return float(interp1d(listx, listy)(value))


def coef_fresnel_TM_pole_aprox(omega,d,k_parallel):
    k0 = omega/c
    eta = epsilon2
#    k_air = k_parallel_air(omega)
#    k_medium = k_parallel_medium(omega,epsilon2)   
    omega_omegaWG = omega/omegaWG
    number_of_poles = math.ceil(omega_omegaWG)    
    q = k_parallel/k0

    list_all_modes = ['TM1A','TM1B','TM2A','TM2B']
    list_poles = []
    for type_of_mode in list_all_modes:  
        rta = function_poles(type_of_mode,omega_omegaWG)
        if rta != None:
            list_poles.append(rta)

    if number_of_poles != len(list_poles):
        raise TypeError('TM: N(omega) != [omega/omegaWG]')
        
    chi_prima = d*np.sqrt(epsilon2*k0**2 - k_parallel**2)
    chi = d*np.sqrt(k_parallel**2 - k0**2)

    term0 = 0
    for j in range(len(list_poles)):
        polo_j = list_poles[j]
#            print(polo_j)
        
        den =  (polo_j*d)**2
        
        Rj_aux = 2*chi_prima*chi/den
        
        term1 = (chi**2 + chi_prima**2)/(chi*chi_prima)
        term2 =  ((eta*chi)**2 + chi_prima**2)/(eta*chi_prima)
        
        Rj_term_final = 1/(2*term1 + term2)
        
        Rj = Rj_aux*Rj_term_final
                
        qj_polo = polo_j/k0
        
        term0 = term0 + Rj*qj_polo/(q - qj_polo)  
        
    return term0


def coef_fresnel_TE_pole_aprox(omega,d,k_parallel):
    k0 = omega/c
    eta = 1
#    k_air = k_parallel_air(omega)
#    k_medium = k_parallel_medium(omega,epsilon2)   
    omega_omegaWG = omega/omegaWG
    number_of_poles = math.ceil(omega_omegaWG)    
    q = k_parallel/k0

    list_all_modes = ['TE1A','TE1B','TE2A','TE2B']
    list_poles = []
    for type_of_mode in list_all_modes:  
        rta = function_poles(type_of_mode,omega_omegaWG)
        if rta != None:
            list_poles.append(rta)

    if number_of_poles != len(list_poles):
        raise TypeError('TM: N(omega) != [omega/omegaWG]')
        
    chi_prima = d*np.sqrt(epsilon2*k0**2 - k_parallel**2)
    chi = d*np.sqrt(k_parallel**2 - k0**2)
    
    term0 = 0
    for j in range(len(list_poles)):
        polo_j = list_poles[j]
#            print(polo_j)
        
        den =  (polo_j*d)**2
        
        Rj_aux = 2*chi_prima*chi/den
        
        term1 = (chi**2 + chi_prima**2)/(chi*chi_prima)
        term2 =  ((eta*chi)**2 + chi_prima**2)/(eta*chi_prima)
        
        Rj_term_final = 1/(2*term1 + term2)
        
        Rj = Rj_aux*Rj_term_final
                
        qj_polo = polo_j/k0
        
        term0 = term0 + Rj*qj_polo/(q - qj_polo)  
        
    return term0

#%%


def coef_fresnel_TM_Fresnel(omega,d,u):
    k0 = omega/c
#    eta = epsilon2
#    k_air = k_parallel_air(omega)
#    k_medium = k_parallel_medium(omega,epsilon2)   

    kz1 = np.sqrt(k0**2 - u**2) if (k0 > u) else 1j*np.sqrt(u**2 - k0**2)
    kz2 = np.sqrt(epsilon2*k0**2 - u**2) if (k0*np.sqrt(epsilon2) > u) else 1j*np.sqrt(u**2 - epsilon2*k0**2)
    
    r12_p = (kz1*epsilon2 - kz2)/(kz1*epsilon2 + kz2)
    r21_p = (kz2 - kz1*epsilon2)/(kz1*epsilon2 + kz2)

    r23_p = (kz2 - kz1*epsilon2)/(kz2 + kz1*epsilon2)

    t12_p = 2*kz1*np.sqrt(epsilon2)/(kz1*epsilon2 + kz2)
    t21_p = 2*kz2*np.sqrt(epsilon2)/(kz1*epsilon2 + kz2)

    expo_coef = np.exp(2*1j*kz2*d)

    rp =  r12_p + t12_p*t21_p*r23_p*expo_coef/(1 - r21_p*r23_p*expo_coef)
        
    return rp


def coef_fresnel_TE_Fresnel(omega,d,u):
    k0 = omega/c
#    eta = 1
#    k_air = k_parallel_air(omega)
#    k_medium = k_parallel_medium(omega,epsilon2)   
    
    kz1 = np.sqrt(k0**2 - u**2) if (k0 > u) else 1j*np.sqrt(u**2 - k0**2)
    kz2 = np.sqrt(epsilon2*k0**2 - u**2) if (k0*np.sqrt(epsilon2) > u) else 1j*np.sqrt(u**2 - epsilon2*k0**2)
    
    r12_p = (kz1 - kz2)/(kz1 + kz2)
    r21_p = (kz2 - kz1)/(kz1 + kz2)

    r23_p = (kz2 - kz1)/(kz2 + kz1)

    t12_p = 2*kz1/(kz1 + kz2)
    t21_p = 2*kz2/(kz1 + kz2)

    expo_coef = np.exp(2*1j*kz2*d)

    rp =  r12_p + t12_p*t21_p*r23_p*expo_coef/(1 - r21_p*r23_p*expo_coef)
        
    return rp


#%%

omega = omega_omegaWG*omegaWG
k1 = omega/c
cte_new = c/omegaWG
#if omega_omegaWG <= 2.5:
#    cota_sup = 410
#    N = 150
#else:
#    cota_sup = 720
#    N = 200

k_air = k_parallel_air(omega)
k_medium = k_parallel_medium(omega,epsilon2)   
N = 200
list_k_parallel = np.linspace(k_air*0.95,k_medium*1.1,N)
list_k_parallel_c_omegaWG = np.array(list_k_parallel)*cte_new

title = title + r', $\omega/\omega_{WG}$ = %.1f' %(omega_omegaWG)
labelpng = '_vs_kparallel_d%inm_omega_omegaWG_%.1f' %(d*1e3,omega_omegaWG)
labelpng2 = '_vs_alphaparallel_d%inm_omega_omegaWG_%.1f' %(d*1e3,omega_omegaWG)

list_y_rpTM_imag = []
list_y_rpTM_real = []
list_y_rsTE_imag = []
list_y_rsTE_real = []

list_y_rpTM_imag_Fresnel = []
list_y_rpTM_real_Fresnel = []
list_y_rsTE_imag_Fresnel = []
list_y_rsTE_real_Fresnel = []

for k_parallel in list_k_parallel: 
    rp_valueTM = coef_fresnel_TM_pole_aprox(omega,d,k_parallel) 
    list_y_rpTM_real.append(np.real(rp_valueTM))
    list_y_rpTM_imag.append(np.imag(rp_valueTM))

    rs_valueTE = coef_fresnel_TE_pole_aprox(omega,d,k_parallel) 
    list_y_rsTE_real.append(np.real(rs_valueTE))
    list_y_rsTE_imag.append(np.imag(rs_valueTE))


    rp_valueTM_Fresnel = coef_fresnel_TM_Fresnel(omega,d,k_parallel) 
    list_y_rpTM_real_Fresnel.append(np.real(rp_valueTM_Fresnel))
    list_y_rpTM_imag_Fresnel.append(np.imag(rp_valueTM_Fresnel))

    rs_valueTE_Fresnel = coef_fresnel_TE_Fresnel(omega,d,k_parallel) 
    list_y_rsTE_real_Fresnel.append(np.real(rs_valueTE_Fresnel))
    list_y_rsTE_imag_Fresnel.append(np.imag(rs_valueTE_Fresnel))


#%%

maxi_rp1_real, mini_rp1_real = np.max(list_y_rpTM_real), np.min(list_y_rpTM_real)
maxi_rp2_real, mini_rp2_real = np.max(list_y_rpTM_real_Fresnel), np.min(list_y_rpTM_real_Fresnel)
maxi_rp_real = np.max([maxi_rp1_real,maxi_rp2_real])
mini_rp_real = np.min([mini_rp1_real,mini_rp2_real])
ejey_rp_real = np.linspace(mini_rp_real,maxi_rp_real,10)

maxi_rp1_imag, mini_rp1_imag = np.max(list_y_rpTM_imag), np.min(list_y_rpTM_imag)
maxi_rp2_imag, mini_rp2_imag = np.max(list_y_rpTM_imag_Fresnel), np.min(list_y_rpTM_imag_Fresnel)
maxi_rp_imag = np.max([maxi_rp1_imag,maxi_rp2_imag])
mini_rp_imag = np.min([mini_rp1_imag,mini_rp2_imag])
ejey_rp_imag = np.linspace(mini_rp_imag,maxi_rp_imag,10)

maxi_rs1_real, mini_rs1_real = np.max(list_y_rsTE_real), np.min(list_y_rsTE_real)
maxi_rs2_real, mini_rs2_real = np.max(list_y_rsTE_real_Fresnel), np.min(list_y_rsTE_real_Fresnel)
maxi_rs_real = np.max([maxi_rs1_real,maxi_rs2_real])
mini_rs_real = np.min([mini_rs1_real,mini_rs2_real])
ejey_rs_real = np.linspace(mini_rs_real,maxi_rs_real,10)


maxi_rs1_imag, mini_rs1_imag = np.max(list_y_rsTE_imag), np.min(list_y_rsTE_imag)
maxi_rs2_imag, mini_rs2_imag = np.max(list_y_rsTE_imag_Fresnel), np.min(list_y_rsTE_imag_Fresnel)
maxi_rs_imag = np.max([maxi_rs1_imag,maxi_rs2_imag])
mini_rs_imag = np.min([mini_rs1_imag,mini_rs2_imag])
ejey_rs_imag = np.linspace(mini_rs_imag,maxi_rs_imag,10)


#%%

plt.figure(figsize=tamfig)    
plt.title(title,fontsize=tamtitle)
plt.xlabel(r'$k_\parallel$c/$\omega_{WG}$',fontsize=tamletra, labelpad = 0)
plt.ylabel(r'Im($r_p$)',fontsize=tamletra, labelpad = 0)
plt.plot(list_k_parallel_c_omegaWG,list_y_rpTM_imag,'-',color = 'blue',ms = lw, label = 'TM Pole')
plt.plot(list_k_parallel_c_omegaWG,list_y_rpTM_imag_Fresnel,'-.',color = 'purple',ms = lw, label = 'TM Fresnel')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
plt.tick_params(labelsize = tamnum,pad = pad)
plt.plot(np.ones(10)*k_air*cte_new, ejey_rp_imag,color = 'grey')
plt.plot(np.ones(10)*k_medium*cte_new, ejey_rp_imag,color = 'black')
plt.grid(1)
plt.tight_layout()
os.chdir(path_save)
plt.savefig('Im_rp' + labelpng  + '.png')



plt.figure(figsize=tamfig)    
plt.title(title,fontsize=tamtitle)
plt.xlabel(r'$k_\parallel$c/$\omega_{WG}$',fontsize=tamletra, labelpad = 0)
plt.ylabel(r'Re($r_p$)',fontsize=tamletra, labelpad = 0)
plt.plot(list_k_parallel_c_omegaWG,list_y_rpTM_real,'-.',color = 'blue',ms = lw, label = 'TM Pole')
plt.plot(list_k_parallel_c_omegaWG,list_y_rpTM_real_Fresnel,'-.',color = 'purple',ms = lw, label = 'TM Fresnel')
plt.plot(np.ones(10)*k_air*cte_new, ejey_rp_real,color = 'grey')
plt.plot(np.ones(10)*k_medium*cte_new, ejey_rp_real,color = 'black')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
plt.tick_params(labelsize = tamnum,pad = pad)
#plt.yscale('log')
plt.grid(1)
plt.tight_layout()
plt.savefig('Re_rp' + labelpng  + '.png')


#
plt.figure(figsize=tamfig)    
plt.title(title,fontsize=tamtitle)
plt.xlabel(r'$k_\parallel$c/$\omega_{WG}$',fontsize=tamletra, labelpad = 0)
plt.ylabel(r'Im($r_s$)',fontsize=tamletra, labelpad = 0)
plt.plot(list_k_parallel_c_omegaWG,list_y_rsTE_imag,'-',color = 'red',ms = lw, label = 'TE Pole')
plt.plot(list_k_parallel_c_omegaWG,list_y_rsTE_imag_Fresnel,'-.',color = 'darkorange',ms = lw, label = 'TE Fresnel')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
plt.tick_params(labelsize = tamnum,pad = pad)
plt.plot(np.ones(10)*k_air*cte_new, ejey_rs_imag,color = 'grey')
plt.plot(np.ones(10)*k_medium*cte_new, ejey_rs_imag,color = 'black')
plt.grid(1)
plt.tight_layout()
# os.chdir(path_save)
plt.savefig('Im_rs' + labelpng  + '.png')


plt.figure(figsize=tamfig)    
plt.title(title,fontsize=tamtitle)
plt.xlabel(r'$k_\parallel$c/$\omega_{WG}$',fontsize=tamletra, labelpad = 0)
plt.ylabel(r'Re($r_s$)',fontsize=tamletra, labelpad = 0)
plt.plot(list_k_parallel_c_omegaWG,list_y_rsTE_real,'-.',color = 'red',ms = lw, label = 'TE Pole')
plt.plot(list_k_parallel_c_omegaWG,list_y_rsTE_real_Fresnel,'-.',color = 'darkorange',ms = lw, label = 'TE Fresnel')
plt.plot(np.ones(10)*k_air*cte_new, ejey_rs_real,color = 'grey')
plt.plot(np.ones(10)*k_medium*cte_new, ejey_rs_real,color = 'black')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
plt.tick_params(labelsize = tamnum,pad = pad)
plt.grid(1)
plt.tight_layout()
# os.chdir(path_save)
plt.savefig('Re_rs' + labelpng  + '.png')

#%%

list_alpha_parallel = np.array(list_k_parallel)/k1

plt.figure(figsize=tamfig)    
plt.title(title,fontsize=tamtitle)
plt.xlabel(r'$\alpha_\parallel$',fontsize=tamletra, labelpad = 0)
plt.ylabel(r'Im($r_p$)',fontsize=tamletra, labelpad = 0)
plt.plot(list_alpha_parallel,list_y_rpTM_imag,'-',color = 'blue',ms = lw, label = 'TM Pole')
plt.plot(list_alpha_parallel,list_y_rpTM_imag_Fresnel,'-.',color = 'purple',ms = lw, label = 'TM Fresnel')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
plt.tick_params(labelsize = tamnum,pad = pad)
plt.plot(np.ones(10)*k_air/k1, ejey_rp_imag,color = 'grey')
plt.plot(np.ones(10)*k_medium/k1, ejey_rp_imag,color = 'black')
plt.grid(1)
plt.tight_layout()
plt.savefig('Im_rp' + labelpng2  + '.png')



plt.figure(figsize=tamfig)    
plt.title(title,fontsize=tamtitle)
plt.xlabel(r'$\alpha_\parallel$',fontsize=tamletra, labelpad = 0)
plt.ylabel(r'Re($r_p$)',fontsize=tamletra, labelpad = 0)
plt.plot(list_alpha_parallel,list_y_rpTM_real,'-.',color = 'blue',ms = lw, label = 'TM Pole')
plt.plot(list_alpha_parallel,list_y_rpTM_real_Fresnel,'-.',color = 'purple',ms = lw, label = 'TM Fresnel')
plt.plot(np.ones(10)*k_air/(k1), ejey_rp_real,color = 'grey')
plt.plot(np.ones(10)*k_medium/(k1), ejey_rp_real,color = 'black')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
plt.tick_params(labelsize = tamnum,pad = pad)
#plt.yscale('log')
plt.grid(1)
plt.tight_layout()
plt.savefig('Re_rp' + labelpng2  + '.png')


#
plt.figure(figsize=tamfig)    
plt.title(title,fontsize=tamtitle)
plt.xlabel(r'$\alpha_\parallel$',fontsize=tamletra, labelpad = 0)
plt.ylabel(r'Im($r_s$)',fontsize=tamletra, labelpad = 0)
plt.plot(list_alpha_parallel,list_y_rsTE_imag,'-',color = 'red',ms = lw, label = 'TE Pole')
plt.plot(list_alpha_parallel,list_y_rsTE_imag_Fresnel,'-.',color = 'darkorange',ms = lw, label = 'TE Fresnel')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
plt.tick_params(labelsize = tamnum,pad = pad)
plt.plot(np.ones(10)*k_air/k1, ejey_rs_imag,color = 'grey')
plt.plot(np.ones(10)*k_medium/k1, ejey_rs_imag,color = 'black')
plt.grid(1)
plt.tight_layout()
# os.chdir(path_save)
plt.savefig('Im_rs' + labelpng2  + '.png')


plt.figure(figsize=tamfig)    
plt.title(title,fontsize=tamtitle)
plt.xlabel(r'$\alpha_\parallel$',fontsize=tamletra, labelpad = 0)
plt.ylabel(r'Re($r_s$)',fontsize=tamletra, labelpad = 0)
plt.plot(list_alpha_parallel,list_y_rsTE_real,'-.',color = 'red',ms = lw, label = 'TE Pole')
plt.plot(list_alpha_parallel,list_y_rsTE_real_Fresnel,'-.',color = 'darkorange',ms = lw, label = 'TE Fresnel')
plt.plot(np.ones(10)*k_air/k1, ejey_rs_real,color = 'grey')
plt.plot(np.ones(10)*k_medium/k1, ejey_rs_real,color = 'black')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
plt.tick_params(labelsize = tamnum,pad = pad)
plt.grid(1)
plt.tight_layout()
# os.chdir(path_save)
plt.savefig('Re_rs' + labelpng2  + '.png')

