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
from scipy.optimize import minimize   

sns.set()

plot_disp_relation = 1
find_zeros_imag_part = 0
pedir_2_condiciones_rel_disp = 0
plot_rp = 1

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
#path_ctes =  path_basic.replace('/' + 'Silica','')
#print('Importar modulos necesarios para este codigo')
path_save = path_basic + '/' + 'disp_relation_Silica'
try:
    sys.path.insert(1, path_basic)
    from constants import constantes
except ModuleNotFoundError:
    print('constants.py no se encuentra en ' + path_basic)

pi,hb,c,alfac,mu1,mu2,mu3 = constantes()
aux = c*hb

#%%

def k_parallel_WG_TM1(omega,epsilon2,d,k_parallel):
    k0 = omega/c
    chi_prima = d*np.sqrt(epsilon2*k0**2 - k_parallel**2)
    chi = d*np.sqrt(k_parallel**2 - k0**2)   
    return np.abs(np.tan(chi_prima/2) - epsilon2*chi/chi_prima)


def k_parallel_WG_TM2(omega,epsilon2,d,k_parallel):
    k0 = omega/c
    chi_prima = d*np.sqrt(epsilon2*k0**2 - k_parallel**2)
    chi = d*np.sqrt(k_parallel**2 - k0**2)
    return np.abs(1/np.tan(chi_prima/2) + epsilon2*chi/chi_prima)

def k_parallel_WG_TE1(omega,d,k_parallel):
    k0 = omega/c
    chi_prima = d*np.sqrt(epsilon2*k0**2 - k_parallel**2)
    chi = d*np.sqrt(k_parallel**2 - k0**2)
    return np.abs(np.tan(chi_prima/2) - chi/chi_prima)


def k_parallel_WG_TE2(omega,d,k_parallel):
    k0 = omega/c
    chi_prima = d*np.sqrt(epsilon2*k0**2 - k_parallel**2)
    chi = d*np.sqrt(k_parallel**2 - k0**2)
    return np.abs(1/np.tan(chi_prima/2) + chi/chi_prima)

def eq_disp_relation_TE(omega,epsilon2,d,k_parallel):
    
    k0 = omega/c
    
    if k_parallel**2 <= epsilon2*(k0**2):
    
        kz2 = np.sqrt(epsilon2*(k0**2) - k_parallel**2)
    
    else:
        kz2 = 1j*np.sqrt(k_parallel**2 - epsilon2*(k0**2))
    
    if k_parallel**2 <= k0**2 :
    
        kz1 = np.sqrt(k0**2 - k_parallel**2)
    
    else:

        kz1 = 1j*np.sqrt(k_parallel**2 - k0**2)
        
    kz3 = kz1
    
    r21 = (kz2 - kz1)/(kz2 + kz1)

    r23 = (kz2 - kz3)/(kz2 + kz3)
    
    
    return 1 - r21*r23*np.exp(1j*2*kz2*d)


def eq_disp_relation_TM(omega,epsilon2,d,k_parallel):
    
    k0 = omega/c
    
    if k_parallel**2 <= epsilon2*(k0**2):
    
        kz2 = np.sqrt(epsilon2*(k0**2) - k_parallel**2)
    
    else:
        kz2 = 1j*np.sqrt(k_parallel**2 - epsilon2*(k0**2))
    
    if k_parallel**2 <= k0**2 :
    
        kz1 = np.sqrt(k0**2 - k_parallel**2)
    
    else:

        kz1 = 1j*np.sqrt(k_parallel**2 - k0**2)
    
    kz3 = kz1
    
    r21 = (kz2*epsilon2 - kz1)/(kz2*epsilon2 + kz1)

    r23 = (kz2*epsilon2 - kz3)/(kz2*epsilon2 + kz3)
    
    return 1 - r21*r23*np.exp(1j*2*kz2*d)


def k_parallel_air(omega):

    return omega/c

def k_parallel_medium(omega,epsilon2):

    return omega*np.sqrt(epsilon2)/c

#%%
    
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
 
epsilon2 = 12
d = 15*1e-3

hb = 6.582118989999999e-16   #electron volts/seg
alfac = 1/137
c = 299792458000000 #microns/seg

list_colours = ['black', 'gold', 'orange', 'darksalmon', 'brown', 'darkred']

title = '$\epsilon_2$ = %i, d = %inm' %(epsilon2,d*1e3)

if pedir_2_condiciones_rel_disp == 1:
    N = 800
else:
    N = 400
list_omega_omegaWG = np.linspace(0.005,4,N)
omegaWG = (np.pi*c/d)/np.sqrt(epsilon2-1) 
cte_new = c/omegaWG
        
list_k_parallel = np.linspace(0.01,500,N)

list_y_air = []
list_y_medium = []
for omega_omegaWG in list_omega_omegaWG:
    omega = omega_omegaWG*omegaWG
    k_air = k_parallel_air(omega)
    k_medium = k_parallel_medium(omega,epsilon2)    
    list_y_air.append(k_air*cte_new)
    list_y_medium.append(k_medium*cte_new)

zeros_eqTM_re_tot = []
zeros_eqTE_re_tot = []

zeros_eqTM_im_tot = []
zeros_eqTE_im_tot = []

plt.figure(figsize=tamfig)    
plt.title(title,fontsize=tamtitle)
for omega_omegaWG in list_omega_omegaWG : 

    omega = omega_omegaWG*omegaWG

    kz1 = omega/c
    kz2 = np.sqrt(epsilon2)*omega/c
   
    list_y_WG_TM_re = []
    list_y_WG_TE_re = []
    
    zeros_eqTM_re = []
    zeros_eqTE_re = []

    for k_parallel in list_k_parallel:  
        
        eq1 = eq_disp_relation_TM(omega,epsilon2,d,k_parallel)
        eq2 = eq_disp_relation_TE(omega,epsilon2,d,k_parallel)
        if pedir_2_condiciones_rel_disp == 1:
            if np.real(eq1) < 1e-3 and np.abs(np.imag(eq1)) <1e-3:
                if kz1 < k_parallel < kz2:
                    zeros_eqTM_re.append(k_parallel)
                
            if np.real(eq2) < 1e-3 and np.abs(np.imag(eq2)) <1e-3:
                if kz1 < k_parallel < kz2:
                    zeros_eqTE_re.append(k_parallel)

        else:
            
            if np.real(eq1) < 1e-3:
                if kz1 < k_parallel < kz2:
                    zeros_eqTM_re.append(k_parallel)
                    print('TM1:', k_parallel_WG_TM1(omega,epsilon2,d,k_parallel))
                    print('TM2:',k_parallel_WG_TM2(omega,epsilon2,d,k_parallel))
                
            if np.real(eq2) < 1e-3:
                if kz1 < k_parallel < kz2:
                    zeros_eqTE_re.append(k_parallel)
                    print('TE1:',k_parallel_WG_TE1(omega,d,k_parallel))
                    print('TE2:',k_parallel_WG_TE2(omega,d,k_parallel))

                
        list_y_WG_TM_re.append(np.real(eq1))
        list_y_WG_TE_re.append(np.real(eq2))
        
    zeros_eqTM_re_tot.append(zeros_eqTM_re)
    zeros_eqTE_re_tot.append(zeros_eqTE_re)
    
    plt.plot(np.array(zeros_eqTE_re)*cte_new,omega_omegaWG*np.array(np.ones(len(zeros_eqTE_re))),'.',color = 'red',ms = 3)
    plt.plot(np.array(zeros_eqTM_re)*cte_new,omega_omegaWG*np.array(np.ones(len(zeros_eqTM_re))),'.',color = 'blue',ms = 3)

plt.xlabel(r'$k_\parallel$c/$\omega_{WG}$',fontsize=tamletra, labelpad = 0)
plt.ylabel(r'$\omega/\omega_{WG}$',fontsize=tamletra, labelpad = 0)
plt.plot(list_y_air,list_omega_omegaWG,'--',color = 'grey',ms = lw, label = 'air')
plt.plot(list_y_medium,list_omega_omegaWG,'--',color = 'black',ms = lw, label = 'medium')
plt.xlim([-0.5,8.5])
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
plt.tick_params(labelsize = tamnum,pad = pad)
plt.grid(1)
plt.tight_layout()
# os.chdir(path_save)
plt.savefig('disp_relation_Silica_re_tot_d%inm.png' %(d*1e3))
   
#list_mEv = np.array(list_omega*1e12*1e3*hb)
#%%
#
def rp_TM(omega,list_polos,d,k_parallel):
    k0 = omega/c
    kz1 = omega/c
    kz2 = np.sqrt(epsilon2)*omega/c    
    eta = epsilon2
    
    q = k_parallel/k0


    if kz1 < k_parallel < kz2:
        
        chi = d*np.sqrt(k_parallel**2 - k0**2)
        chi_prima = d*np.sqrt(epsilon2*k0**2 - k_parallel**2)
        
        term0 = 0
        for j in range(len(list_polos)):
            polo_j = list_polos[j]
            
            den =  (polo_j*d)**2
            
            Rj_aux = 2*chi_prima*chi/den
            
            term1 = (chi**2 + chi_prima**2)/(chi*chi_prima)
            term2 =  (eta*chi)**2 + chi_prima**2/(eta*chi_prima)
            
            Rj_term_final = 1/(2*term1 + term2)
            
            Rj = Rj_aux*Rj_term_final
                    
            qj_polo = polo_j/k0
            
            if q - qj_polo != 0:
            
                term0 = term0 + Rj*polo_j/(q - qj_polo)  
            
        return term0

    else:
        return 0

#%%
        
def rp_TE(omega,list_polos,d,k_parallel):
    k0 = omega/c
    kz1 = omega/c
    kz2 = np.sqrt(epsilon2)*omega/c    
    eta = 1
    
    q = k_parallel/k0


    if kz1 < k_parallel < kz2:
        
        chi = d*np.sqrt(k_parallel**2 - k0**2)
        chi_prima = d*np.sqrt(epsilon2*k0**2 - k_parallel**2)
        
        term0 = 0
        for j in range(len(list_polos)):
            polo_j = list_polos[j]
            
            den =  (polo_j*d)**2
            
            Rj_aux = 2*chi_prima*chi/den
            
            term1 = (chi**2 + chi_prima**2)/(chi*chi_prima)
            term2 =  (eta*chi)**2 + chi_prima**2/(eta*chi_prima)
            
            Rj_term_final = 1/(2*term1 + term2)
            
            Rj = Rj_aux*Rj_term_final
            
            qj_polo = polo_j/k0
            if q - qj_polo != 0:  
                      
                term0 = term0 + Rj*polo_j/(q - qj_polo)   
            
        return term0

    else:
        return 0

#%%
        
index = int(N*0.8)
omega_omegaWG = list_omega_omegaWG[index]
omega = omega_omegaWG*omegaWG
list_polos_TM_re = zeros_eqTM_re_tot[index]
#list_polos_TM_im = zeros_eqTM_im_tot[index]
title2 = title + r', $\omega/\omega_{WG}$ = %i' %(omega_omegaWG)
kz1 = omega/c
kz2 = np.sqrt(epsilon2)*omega/c

list_y_rpTM_imag = []
list_y_rpTM_real = []
for k_parallel in list_k_parallel:  
    
    rp_valueTM = rp_TM(omega,list_polos_TM_re,d,k_parallel) 
    list_y_rpTM_real.append(np.real(rp_valueTM))
    list_y_rpTM_imag.append(np.imag(rp_valueTM))

#plt.figure(figsize=tamfig)    
#plt.title(title,fontsize=tamtitle)
#plt.xlabel(r'$k_\parallel$ [1/$\mu$m]',fontsize=tamletra, labelpad = 0)
#plt.ylabel(r'Im($r_p$)',fontsize=tamletra, labelpad = 0)
#plt.plot(list_k_parallel,list_y_rpTM_imag,'.',color = 'blue',ms = lw, label = 'TM')
#plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
#plt.tick_params(labelsize = tamnum,pad = pad)
#plt.grid(1)
#plt.tight_layout()
## os.chdir(path_save)
#plt.savefig('Im_rp_vs_kparallel_TM_d%inm.png' %(d*1e3))


plt.figure(figsize=tamfig)    
plt.title(title2,fontsize=tamtitle)
plt.xlabel(r'$k_\parallel$ [1/$\mu$m]',fontsize=tamletra, labelpad = 0)
plt.ylabel(r'Re($r_p$)',fontsize=tamletra, labelpad = 0)
plt.plot(list_k_parallel,list_y_rpTM_real,'-.',color = 'blue',ms = lw, label = 'TM')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
plt.tick_params(labelsize = tamnum,pad = pad)
plt.grid(1)
plt.tight_layout()
# os.chdir(path_save)
plt.savefig('Re_rp_vs_kparallel_TM_d%inm.png' %(d*1e3))

#%%

index =  int(N*0.8)
omega_omegaWG = list_omega_omegaWG[index]
omega = omega_omegaWG*omegaWG
list_polos_TE = zeros_eqTE_re_tot[index]
kz1 = omega/c
kz2 = np.sqrt(epsilon2)*omega/c

list_y_rpTE_imag = []
list_y_rpTE_real = []
for k_parallel in list_k_parallel:  
    
    rp_valueTE = rp_TE(omega,list_polos_TE,d,k_parallel) 
    list_y_rpTE_real.append(np.real(rp_valueTE))
    list_y_rpTE_imag.append(np.imag(rp_valueTE))

#plt.figure(figsize=tamfig)    
#plt.title(title,fontsize=tamtitle)
#plt.xlabel(r'$k_\parallel$ [1/$\mu$m]',fontsize=tamletra, labelpad = 0)
#plt.ylabel(r'Im($r_s$)',fontsize=tamletra, labelpad = 0)
#plt.plot(list_k_parallel,list_y_rpTE_imag,'-.',color = 'red',ms = lw, label = 'TE')
#plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
#plt.tick_params(labelsize = tamnum,pad = pad)
#plt.grid(1)
#plt.tight_layout()
## os.chdir(path_save)
#plt.savefig('Im_rs_vs_kparallel_TM_d%inm.png' %(d*1e3))


plt.figure(figsize=tamfig)    
plt.title(title2,fontsize=tamtitle)
plt.xlabel(r'$k_\parallel$ [1/$\mu$m]',fontsize=tamletra, labelpad = 0)
plt.ylabel(r'Re($r_s$)',fontsize=tamletra, labelpad = 0)
plt.plot(list_k_parallel,list_y_rpTE_real,'-.',color = 'red',ms = lw, label = 'TE')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
plt.tick_params(labelsize = tamnum,pad = pad)
plt.grid(1)
plt.tight_layout()
# os.chdir(path_save)
plt.savefig('Re_rs_vs_kparallel_TM_d%inm.png' %(d*1e3))

#%%

if find_zeros_imag_part == 1:
    list_omega_omegaWG = np.linspace(0.005,4,100)
    omegaWG = (np.pi*c/d)/np.sqrt(epsilon2-1) 
    cte_new = c/omegaWG
            
    #omega = 2*omegaWG
    #k0 = omega/c
    N = 100
    list_k_parallel = np.linspace(0,500,N)
    
    
    list_y_air = []
    list_y_medium = []
    for omega_omegaWG in list_omega_omegaWG:
        omega = omega_omegaWG*omegaWG
        k_air = k_parallel_air(omega)
        k_medium = k_parallel_medium(omega,epsilon2)    
        list_y_air.append(k_air*cte_new)
        list_y_medium.append(k_medium*cte_new)
    
    zeros_eqTM_im_tot = []
    zeros_eqTE_im_tot = []
    
    plt.figure(figsize=tamfig)    
    plt.title(title,fontsize=tamtitle)
    for omega_omegaWG in list_omega_omegaWG : 
    
        omega = omega_omegaWG*omegaWG
    
        kz1 = omega/c
        kz2 = np.sqrt(epsilon2)*omega/c
        
        
        list_y_WG_TM_im = []
        list_y_WG_TE_im = []
        
        zeros_eqTM_im = []
        zeros_eqTE_im = []
            
        for k_parallel in list_k_parallel:  
            
            eq1 = eq_disp_relation_TM(omega,epsilon2,d,k_parallel)
            eq2 = eq_disp_relation_TE(omega,epsilon2,d,k_parallel)
        
            if np.imag(eq1) < 1e-15:
                if kz1 < k_parallel < kz2:
                    zeros_eqTM_im.append(k_parallel)
                
            if np.imag(eq2) < 1e-15:
                if kz1 < k_parallel < kz2:
                    zeros_eqTE_im.append(k_parallel)
    
        
            list_y_WG_TM_im.append(np.imag(eq1))
            list_y_WG_TE_im.append(np.imag(eq2))
            
        zeros_eqTM_re_tot.append(zeros_eqTM_im)
        zeros_eqTE_re_tot.append(zeros_eqTE_im)
        
        plt.plot(np.array(zeros_eqTE_im)*cte_new,omega_omegaWG*np.array(np.ones(len(zeros_eqTE_im))),'.',color = 'red',ms = 3)
        plt.plot(np.array(zeros_eqTM_im)*cte_new,omega_omegaWG*np.array(np.ones(len(zeros_eqTM_im))),'.',color = 'blue',ms = 3)
    
    plt.xlabel(r'$k_\parallel$c/$\omega_{WG}$',fontsize=tamletra, labelpad = 0)
    plt.ylabel(r'$\omega/\omega_{WG}$',fontsize=tamletra, labelpad = 0)
    plt.plot(list_y_air,list_omega_omegaWG,'--',color = 'grey',ms = lw, label = 'air')
    plt.plot(list_y_medium,list_omega_omegaWG,'--',color = 'black',ms = lw, label = 'medium')
    plt.xlim([-0.5,8.5])
    plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
    plt.tick_params(labelsize = tamnum,pad = pad)
    plt.grid(1)
    plt.tight_layout()
    # os.chdir(path_save)
    plt.savefig('disp_relation_Silica_im_tot_d%inm.png' %(d*1e3))

#%%

if plot_disp_relation == 1:
    omega_omegaWG = 2
    title2 = title + r', $\omega/\omega_{WG}$ = %i' %(omega_omegaWG)
    omega = omega_omegaWG*omegaWG
    
    list_y_WG_TM_re = []
    list_y_WG_TE_re = []
    
    
    list_y_WG_TM_im = []
    list_y_WG_TE_im = []
    
    zeros_eqTM_re = []
    zeros_eqTE_re = []
    
    kz1 = omega/c
    kz2 = np.sqrt(epsilon2)*omega/c
        
    for k_parallel in list_k_parallel:  
        
        eq1 = eq_disp_relation_TM(omega,epsilon2,d,k_parallel)
        eq2 = eq_disp_relation_TE(omega,epsilon2,d,k_parallel)
    
        if np.real(eq1) < 1e-3:
            zeros_eqTM_re.append(k_parallel)
            
        if np.real(eq2) < 1e-3:
            zeros_eqTE_re.append(k_parallel)
            
    
        list_y_WG_TM_re.append(np.real(eq1))
        list_y_WG_TE_re.append(np.real(eq2))
    
        list_y_WG_TM_im.append(np.imag(eq1))
        list_y_WG_TE_im.append(np.imag(eq2))
        
    maxi_re = np.max([np.max(list_y_WG_TE_re),np.max(list_y_WG_TM_re)])
    mini_re = np.min([np.min(list_y_WG_TE_re),np.min(list_y_WG_TM_re)])
    
    maxi_im = np.max([np.max(list_y_WG_TE_im),np.max(list_y_WG_TM_im)])
    mini_im = np.min([np.min(list_y_WG_TE_im),np.min(list_y_WG_TM_im)])
    
    plt.figure(figsize=tamfig)    
    plt.title(title2,fontsize=tamtitle)
    plt.plot(list_k_parallel,list_y_WG_TE_re,'.-',color = 'red',ms = 3, label = 'WG TE')
    plt.plot(zeros_eqTM_re,np.array(np.ones(len(zeros_eqTM_re)))*0,'.',color = 'hotpink',ms = 5)
    plt.plot(zeros_eqTE_re,np.array(np.ones(len(zeros_eqTE_re)))*0,'.',color = 'hotpink',ms = 5)
    plt.plot(np.array(np.ones(20)*kz1),np.linspace(mini_re,maxi_re,20),'-')
    plt.plot(np.array(np.ones(20)*kz2),np.linspace(mini_re,maxi_re,20),'-')
    plt.plot(list_k_parallel,list_y_WG_TM_re,'.-',color = 'blue',ms = 3, label = 'WG TM')
    plt.xlabel(r'$k_\parallel$ [1/$\mu$m]',fontsize=tamletra, labelpad = 0)
    plt.ylabel(r'Re(Eq dispersion)',fontsize=tamletra, labelpad = 0)
    plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
    plt.tick_params(labelsize = tamnum,pad = pad)
    plt.grid(1)
    plt.tight_layout()
    # os.chdir(path_save)
    plt.savefig('disp_relation_Silica_real_vs_kparallel_d%inm.png' %(d*1e3))
    
    plt.figure(figsize=tamfig)    
    plt.title(title2,fontsize=tamtitle)
    plt.plot(list_k_parallel,list_y_WG_TE_im,'.-',color = 'red',ms = 3, label = 'WG TE')
    plt.plot(list_k_parallel,list_y_WG_TM_im,'.-',color = 'blue',ms = 3, label = 'WG TM')
    plt.plot(np.array(np.ones(20)*kz1),np.linspace(mini_im,maxi_im,20),'-')
    plt.plot(np.array(np.ones(20)*kz2),np.linspace(mini_im,maxi_im,20),'-')
    
    plt.xlabel(r'$k_\parallel$ [1/$\mu$m]',fontsize=tamletra, labelpad = 0)
    plt.ylabel(r'Im(Eq dispersion)',fontsize=tamletra, labelpad = 0)
    plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
    plt.tick_params(labelsize = tamnum,pad = pad)
    plt.grid(1)
    plt.tight_layout()
    # os.chdir(path_save)
    plt.savefig('disp_relation_Silica_imag_vs_kparallel_d%inm.png' %(d*1e3))
