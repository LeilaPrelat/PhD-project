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
from scipy.optimize import minimize,fsolve 
 

sns.set()

pedir_2_condiciones_rel_disp = 0
minimizar_formulas_Eduardo = 1

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_ctes =  path_basic.replace('/' + 'new','')
#print('Importar modulos necesarios para este codigo')
path_save = path_basic + '/' + 'disp_relation_Silica'
try:
    sys.path.insert(1, path_ctes)
    from constants import constantes
except ModuleNotFoundError:
    print('constants.py no se encuentra en ' + path_ctes)

pi,hb,c,alfac,mu1,mu2,mu3 = constantes()
aux = c*hb

 
epsilon2 = 12
d = 20*1e-3

#%%

def k_parallel_air(omega):

    return omega/c

def k_parallel_medium(omega,epsilon2):

    return omega*np.sqrt(epsilon2)/c

#%%
    
def k_parallel_WG_TM1(omega,epsilon2,d,k_parallel):
    k0 = omega/c
    
    if k_parallel**2 <= epsilon2*(k0**2):    
        chi_prima = d*np.sqrt(epsilon2*k0**2 - k_parallel**2)
    else:
        chi_prima = d*1j*np.sqrt(k_parallel**2 - epsilon2*k0**2)

    if k_parallel**2 <= k0**2 :
        chi = d*1j*np.sqrt(k0**2 - k_parallel**2)
    else:
        chi = d*np.sqrt(k_parallel**2 - k0**2)
    
    return np.tan(chi_prima/2) - epsilon2*chi/chi_prima


def k_parallel_WG_TM2(omega,epsilon2,d,k_parallel):
    k0 = omega/c
    
    if k_parallel**2 <= epsilon2*(k0**2):    
        chi_prima = d*np.sqrt(epsilon2*k0**2 - k_parallel**2)
    else:
        chi_prima = d*1j*np.sqrt(k_parallel**2 - epsilon2*k0**2)

    if k_parallel**2 <= k0**2 :
        chi = d*1j*np.sqrt(k0**2 - k_parallel**2)
    else:
        chi = d*np.sqrt(k_parallel**2 - k0**2)
    
    return 1/np.tan(chi_prima/2) + epsilon2*chi/chi_prima

def k_parallel_WG_TE1(omega,epsilon2,d,k_parallel):
    k0 = omega/c
    
    if k_parallel**2 <= epsilon2*(k0**2):    
        chi_prima = d*np.sqrt(epsilon2*k0**2 - k_parallel**2)
    else:
        chi_prima = d*1j*np.sqrt(k_parallel**2 - epsilon2*k0**2)

    if k_parallel**2 <= k0**2 :
        chi = d*1j*np.sqrt(k0**2 - k_parallel**2)
    else:
        chi = d*np.sqrt(k_parallel**2 - k0**2)
    
    
    return np.tan(chi_prima/2) - chi/chi_prima


def k_parallel_WG_TE2(omega,epsilon2,d,k_parallel):
    k0 = omega/c
    
    if k_parallel**2 <= epsilon2*(k0**2):    
        chi_prima = d*np.sqrt(epsilon2*k0**2 - k_parallel**2)
    else:
        chi_prima = d*1j*np.sqrt(k_parallel**2 - epsilon2*k0**2)

    if k_parallel**2 <= k0**2 :
        chi = d*1j*np.sqrt(k0**2 - k_parallel**2)
    else:
        chi = d*np.sqrt(k_parallel**2 - k0**2)
    
    return 1/np.tan(chi_prima/2) + chi/chi_prima

#%%

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

#%%
    
def rp_TM(omega,list_polos,d,k_parallel):
    k0 = omega/c
#    kz1 = omega/c
#    kz2 = np.sqrt(epsilon2)*omega/c    
    eta = epsilon2
    
    q = k_parallel/k0
        
    if k_parallel**2 <= epsilon2*(k0**2):    
        chi_prima = d*np.sqrt(epsilon2*k0**2 - k_parallel**2)
    else:
        chi_prima = d*1j*np.sqrt(k_parallel**2 - epsilon2*k0**2)

    if k_parallel**2 <= k0**2 :
        chi = d*1j*np.sqrt(k0**2 - k_parallel**2)
    else:
        chi = d*np.sqrt(k_parallel**2 - k0**2)
    
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
        
def rs_TE(omega,list_polos,d,k_parallel):
    k0 = omega/c
#    kz1 = omega/c
#    kz2 = np.sqrt(epsilon2)*omega/c    
    eta = 1
    
    q = k_parallel/k0
        
    if k_parallel**2 <= epsilon2*(k0**2):    
        chi_prima = d*np.sqrt(epsilon2*k0**2 - k_parallel**2)
    else:
        chi_prima = d*1j*np.sqrt(k_parallel**2 - epsilon2*k0**2)

    if k_parallel**2 <= k0**2 :
        chi = d*1j*np.sqrt(k0**2 - k_parallel**2)
    else:
        chi = d*np.sqrt(k_parallel**2 - k0**2)
    
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


title = '$\epsilon_2$ = %i, d = %inm' %(epsilon2,d*1e3)




if pedir_2_condiciones_rel_disp == 1:
    N = 800
else:
    N = 500
    
cota_inf_omega_omegaWG = 0.00001
cota_sup_omega_omegaWG = 4
list_omega_omegaWGA = np.linspace(cota_inf_omega_omegaWG,cota_sup_omega_omegaWG,N)
#list_omega_omegaWG = np.linspace(4,0.005,N)
omegaWG = (np.pi*c/d)/(np.sqrt(epsilon2-1)) 
cte_new = c/omegaWG
        
list_k_parallel = np.linspace(0.000001,800,N)

list_y_air = []
list_y_medium = []

zeros_eqTM1A = []
zeros_eqTE1A = []

zeros_eqTM2A = []
zeros_eqTE2A = []

init_condit_TM1A = []
init_condit_TM2A = []
init_condit_TE1A = []
init_condit_TE2A = []

list_omega_omegaWG_TM1A = []
list_omega_omegaWG_TM2A = []
list_omega_omegaWG_TE1A = []
list_omega_omegaWG_TE2A = []

k_TM1 = 0
k_TM2 = 0
k_TE1 = 0
k_TE2 = 0
for omega_omegaWG in list_omega_omegaWGA:
    omega = omega_omegaWG*omegaWG
    k_air = k_parallel_air(omega)
    k_medium = k_parallel_medium(omega,epsilon2)   
    
    print(k_air)
    
    list_y_air.append(k_air*cte_new)
    list_y_medium.append(k_medium*cte_new)

    init_condit = [k_air]
    if minimizar_formulas_Eduardo == 1:
        
        def k_parallel_WG_TE2_1var(k_parallel):   
            return np.abs(k_parallel_WG_TE2(omega,epsilon2,d,k_parallel))
    
        def k_parallel_WG_TE1_1var(k_parallel):   
            return  np.abs(k_parallel_WG_TE1(omega,epsilon2,d,k_parallel))
    
        def k_parallel_WG_TM2_1var(k_parallel):   
            return  np.abs(k_parallel_WG_TM2(omega,epsilon2,d,k_parallel))
    
        def k_parallel_WG_TM1_1var(k_parallel):   
            return  np.abs(k_parallel_WG_TM1(omega,epsilon2,d,k_parallel))

        if k_TM1 == 0 :
            resTM1 = fsolve(k_parallel_WG_TM1_1var, [k_air])
        else:
            resTM1 = fsolve(k_parallel_WG_TM1_1var, [init_condit_TM1A[k_TM1-1]], maxfev=1500 )         
        resTM1_v = np.float(resTM1)
        if k_air < resTM1_v < k_medium:
            zeros_eqTM1A.append(resTM1_v)
            list_omega_omegaWG_TM1A.append(omega_omegaWG)
            k_TM1 = k_TM1 + 1
#            if resTM1_v*cte_new < 7:
#                init_condit_TM1A.append(resTM1_v)
#            else:
#                init_condit_TM1A.append(k_medium)
            
            if omega_omegaWG < 2.5:
                init_condit_TM1A.append(resTM1_v)
            else:
                init_condit_TM1A.append(k_medium)

        if k_TM2 == 0 :
            resTM2 = fsolve(k_parallel_WG_TM2_1var, [k_air])
        else:
            resTM2 = fsolve(k_parallel_WG_TM2_1var, [init_condit_TM2A[k_TM2-1]])
        resTM2_v = np.float(resTM2)            
        if k_air < resTM2_v < k_medium :
            zeros_eqTM2A.append(resTM2_v)
            list_omega_omegaWG_TM2A.append(omega_omegaWG)
            k_TM2 = k_TM2 + 1
            if resTM2_v*cte_new < 14:
                init_condit_TM2A.append(resTM2_v)
            else:
                init_condit_TM2A.append(k_air)
        

        if k_TE1 == 0 :
            resTE1 = fsolve(k_parallel_WG_TE1_1var, [k_air])
        else:
            resTE1 = fsolve(k_parallel_WG_TE1_1var, [init_condit_TE1A[k_TE1-1]])
        resTE1_v = np.float(resTE1)            
        if k_air <= resTE1_v <= k_medium:
            zeros_eqTE1A.append(resTE1_v)
            list_omega_omegaWG_TE1A.append(omega_omegaWG)
            k_TE1 = k_TE1 + 1

            if resTE1_v*cte_new < 14:
                init_condit_TE1A.append(resTE1_v)
            else:
                init_condit_TE1A.append(k_air)



        if k_TE2 == 0 :
            resTE2 = fsolve(k_parallel_WG_TE2_1var, [k_air])
        else:
            resTE2 = fsolve(k_parallel_WG_TE2_1var, [init_condit_TE2A[k_TE2-1]])
        resTE2_v = np.float(resTE2)            
        if k_air <= resTE2_v <= k_medium:
            zeros_eqTE2A.append(resTE2_v)
            list_omega_omegaWG_TE2A.append(omega_omegaWG)
            k_TE2 = k_TE2 + 1

            if resTE2_v*cte_new < 14:
                init_condit_TE2A.append(resTE2_v)
            else:
                init_condit_TE2A.append(k_air)
    
#%%


number_of_sol = 2*int(omega_omegaWG)




#%%

plt.figure(figsize=tamfig)    
plt.title(title,fontsize=tamtitle)
plt.plot(np.array(zeros_eqTE2A)*cte_new,list_omega_omegaWG_TE2A,'.',color = 'red',ms = 4, label = 'TE 2')
plt.plot(np.array(zeros_eqTE1A)*cte_new,list_omega_omegaWG_TE1A,'.',color = 'orange',ms = 4, label = 'TE 1')
plt.plot(np.array(zeros_eqTM2A)*cte_new,list_omega_omegaWG_TM2A,'.',color = 'blue',ms = 4, label = 'TM 2')
plt.plot(np.array(zeros_eqTM1A)*cte_new,list_omega_omegaWG_TM1A,'.',color = 'lightblue',ms = 4, label = 'TM 1')

#plt.plot(np.array(zeros_eqTE2B)*cte_new,list_omega_omegaWG_TE2B,'.',color = 'red',ms = 4)
#plt.plot(np.array(zeros_eqTE1B)*cte_new,list_omega_omegaWG_TE1B,'.',color = 'red',ms = 4)
#plt.plot(np.array(zeros_eqTM2B)*cte_new,list_omega_omegaWG_TM2B,'.',color = 'blue',ms = 4)
#plt.plot(np.array(zeros_eqTM1B)*cte_new,list_omega_omegaWG_TM1B,'.',color = 'blue',ms = 4)

plt.xlabel(r'$k_\parallel$c/$\omega_{WG}$',fontsize=tamletra, labelpad = 0)
plt.ylabel(r'$\omega/\omega_{WG}$',fontsize=tamletra, labelpad = 0)
plt.plot(list_y_air,list_omega_omegaWGA,'--',color = 'grey',ms = lw, label = 'air')
plt.plot(list_y_medium,list_omega_omegaWGA,'--',color = 'black',ms = lw, label = 'medium')
plt.xlim([-0.5,8.5])
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
plt.tick_params(labelsize = tamnum,pad = pad)
plt.grid(1)
plt.tight_layout()
# os.chdir(path_save)
plt.savefig('disp_relation_Silica_Edu_solutions_d%inm.png' %(d*1e3))
   
#list_mEv = np.array(list_omega*1e12*1e3*hb)
#%%
list_k_parallel = np.linspace(200,400)

indexB = int(len(list_omega_omegaWGA)/2)
indexA = indexB - 1
omega_omegaWGA = list_omega_omegaWGA[indexA]
omega_omegaWGB = list_omega_omegaWGB[indexB]

omegaA = omega_omegaWGA*omegaWG
omegaB = omega_omegaWGB*omegaWG

k_air = k_parallel_air(omegaA)
k_medium = k_parallel_medium(omegaA,epsilon2)    
   
list_polos = []    
if k_air <= zeros_eqTM1A[indexA] <= k_medium:  
    list_polos.append(zeros_eqTM1A[indexA])

if k_air <= zeros_eqTM1B[indexB] <= k_medium:  
    list_polos.append(zeros_eqTM1B[indexB])


if k_air <= zeros_eqTM2A[indexA] <= k_medium:  
    list_polos.append(zeros_eqTM2A[indexA])
    
if k_air <= zeros_eqTM2B[indexB] <= k_medium:  
    list_polos.append(zeros_eqTM2B[indexB])
    
    
#    list_polos = [zeros_eqTM1A[indexA],zeros_eqTM2A[indexA]]    
#list_polos_TM_im = zeros_eqTM_im_tot[index]
title2 = title + r', $\omega/\omega_{WG}$ = %.2f' %(omega_omegaWGA)

list_y_rpTM_real = []
list_y_rpTM_imag = []
for k_parallel in list_k_parallel:  
    
    rp_valueTM = rp_TM(omega,list_polos,d,k_parallel) 
    list_y_rpTM_real.append(np.real(rp_valueTM))
    list_y_rpTM_imag.append(np.imag(rp_valueTM))        
            
    
#    list_y_rpTM_imag.append(np.imag(rp_valueTM))
maxi_re,mini_re = np.max(list_y_rpTM_real),np.min(list_y_rpTM_real)
plt.figure(figsize=tamfig)    
plt.title(title2,fontsize=tamtitle)
plt.xlabel(r'$k_\parallel$ [1/$\mu$m]',fontsize=tamletra, labelpad = 0)
plt.ylabel(r'Re($r_p$)',fontsize=tamletra, labelpad = 0)
for value in list_polos:
    plt.plot(value*np.array(np.ones(10)),np.linspace(mini_re,maxi_re,10))
plt.plot(list_k_parallel,list_y_rpTM_real,'-.',color = 'blue',ms = lw, label = 'TM')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
plt.tick_params(labelsize = tamnum,pad = pad)
plt.grid(1)
plt.tight_layout()
# os.chdir(path_save)
plt.savefig('Re_rp_vs_kparallel_TM_d%inm.png' %(d*1e3))

maxi_im,mini_im = np.max(list_y_rpTM_imag),np.min(list_y_rpTM_imag)
plt.figure(figsize=tamfig)    
plt.title(title2,fontsize=tamtitle)
plt.xlabel(r'$k_\parallel$ [1/$\mu$m]',fontsize=tamletra, labelpad = 0)
plt.ylabel(r'Im($r_p$)',fontsize=tamletra, labelpad = 0)
for value in list_polos:
    plt.plot(value*np.array(np.ones(10)),np.linspace(mini_im,maxi_im,10))
plt.plot(list_k_parallel,list_y_rpTM_imag,'-.',color = 'blue',ms = lw, label = 'TM')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
plt.tick_params(labelsize = tamnum,pad = pad)
plt.grid(1)
plt.tight_layout()
# os.chdir(path_save)
plt.savefig('Im_rp_vs_kparallel_TM_d%inm.png' %(d*1e3))


#%%

list_k_parallel = np.linspace(200,400)

indexB = int(len(list_omega_omegaWGA)/2)
indexA = indexB - 1
omega_omegaWGA = list_omega_omegaWGA[indexA]
omega_omegaWGB = list_omega_omegaWGB[indexB]

omegaA = omega_omegaWGA*omegaWG
omegaB = omega_omegaWGB*omegaWG

k_air = k_parallel_air(omegaA)
k_medium = k_parallel_medium(omegaA,epsilon2)    
   
list_polos = []    
if k_air <= zeros_eqTE1A[indexA] <= k_medium:  
    list_polos.append(zeros_eqTE1A[indexA])

if k_air <= zeros_eqTE1B[indexB] <= k_medium:  
    list_polos.append(zeros_eqTE1B[indexB])


if k_air <= zeros_eqTE2A[indexA] <= k_medium:  
    list_polos.append(zeros_eqTE2A[indexA])
    
if k_air <= zeros_eqTE2B[indexB] <= k_medium:  
    list_polos.append(zeros_eqTE2B[indexB])
    
    
#    list_polos = [zeros_eqTM1A[indexA],zeros_eqTM2A[indexA]]    
#list_polos_TM_im = zeros_eqTM_im_tot[index]
title2 = title + r', $\omega/\omega_{WG}$ = %i' %(omega_omegaWG)

list_y_rpTE_real = []
list_y_rpTE_imag = []
for k_parallel in list_k_parallel:  
    
    rs_valueTE = rs_TE(omega,list_polos,d,k_parallel) 
    list_y_rpTE_real.append(np.real(rs_valueTE))
    list_y_rpTE_imag.append(np.imag(rs_valueTE))        
            
    
#    list_y_rpTM_imag.append(np.imag(rp_valueTM))
maxi_re,mini_re = np.max(list_y_rpTE_real),np.min(list_y_rpTE_real)
plt.figure(figsize=tamfig)    
plt.title(title2,fontsize=tamtitle)
plt.xlabel(r'$k_\parallel$ [1/$\mu$m]',fontsize=tamletra, labelpad = 0)
plt.ylabel(r'Re($r_s$)',fontsize=tamletra, labelpad = 0)
for value in list_polos:
    plt.plot(value*np.array(np.ones(10)),np.linspace(mini_re,maxi_re,10))
plt.plot(list_k_parallel,list_y_rpTE_real,'-.',color = 'red',ms = lw, label = 'TE')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
plt.tick_params(labelsize = tamnum,pad = pad)
plt.grid(1)
plt.tight_layout()
# os.chdir(path_save)
plt.savefig('Re_rs_vs_kparallel_TM_d%inm.png' %(d*1e3))

maxi_im,mini_im = np.max(list_y_rpTE_imag),np.min(list_y_rpTE_imag)
plt.figure(figsize=tamfig)    
plt.title(title2,fontsize=tamtitle)
plt.xlabel(r'$k_\parallel$ [1/$\mu$m]',fontsize=tamletra, labelpad = 0)
plt.ylabel(r'Im($r_s$)',fontsize=tamletra, labelpad = 0)
for value in list_polos:
    plt.plot(value*np.array(np.ones(10)),np.linspace(mini_im,maxi_im,10))
plt.plot(list_k_parallel,list_y_rpTE_imag,'-.',color = 'red',ms = lw, label = 'TE')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
plt.tick_params(labelsize = tamnum,pad = pad)
plt.grid(1)
plt.tight_layout()
# os.chdir(path_save)
plt.savefig('Im_rs_vs_kparallel_TM_d%inm.png' %(d*1e3))



#%%
