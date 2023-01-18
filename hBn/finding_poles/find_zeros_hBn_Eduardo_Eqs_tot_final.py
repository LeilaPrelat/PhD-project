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
from scipy.optimize import fsolve 

sns.set()

save_plots = 1
save_txt = 0

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_ctes =  path_basic.replace('/' + 'finding_poles','')
#print('Importar modulos necesarios para este codigo')
path_save = path_basic + '/' + 'disp_relation_hBn'
try:
    sys.path.insert(1, path_ctes)
    from constants import constantes
except ModuleNotFoundError:
    print('constants.py no se encuentra en ' + path_ctes)

pi,hb,c,alfac,mu1,mu2 = constantes()
aux = c*hb



#%%

def epsilon_x(hbw):
    
    epsi_inf = 4.87
    hbgamma = 0.87*1e-3
    f_x = 1.83
    
    num = (170.1*1e-3)**2
    
    den = hbw*(hbw + 1j*hbgamma ) - num
    
    return epsi_inf - f_x*(num/den)


def epsilon_z(hbw):
    
    epsi_inf = 2.95
    hbgamma = 0.25*1e-3
    f_x = 0.61
    
    num = (92.5*1e-3)**2
    
    den = hbw*(hbw + 1j*hbgamma ) - num
    
    return epsi_inf - f_x*(num/den)



def k_parallel_air(hbw):
    
    omegac = hbw/(hb*c)

    return omegac

def k_parallel_medium_x(hbw):
    
    omegac = hbw/(hb*c)
    
    
    epsilon2 =  epsilon_x(hbw)

    return omegac*np.sqrt(epsilon2)


def k_parallel_medium_z(hbw):
    
    omegac = hbw/(hb*c)
    
    
    epsilon2 =  epsilon_z(hbw)

    return omegac*np.sqrt(epsilon2)


#%%

def k_parallel_WG_TM2B_x(hbw,d,k_parallel):
    omegac = hbw/(hb*c)
    k0 = omegac
    
    epsilon2 =  epsilon_x(hbw)    
    if k_parallel**2 <= epsilon2*(k0**2):    
        chi_prima = d*np.sqrt(epsilon2*k0**2 - k_parallel**2)
    else:
        chi_prima = d*1j*np.sqrt(k_parallel**2 - epsilon2*k0**2)

    if k_parallel**2 <= k0**2 :
        chi = d*1j*np.sqrt(k0**2 - k_parallel**2)
    else:
        chi = d*np.sqrt(k_parallel**2 - k0**2)
    
    return 1/np.tan(chi_prima/2) + epsilon2*chi/chi_prima


#def k_parallel_WG_TM2B_z(hbw,d,k_parallel):
#    omegac = hbw/(hb*c)
#    k0 = omegac
#    
#    epsilon2 =  epsilon_z(hbw)    
#    if k_parallel**2 <= epsilon2*(k0**2):    
#        chi_prima = d*np.sqrt(epsilon2*k0**2 - k_parallel**2)
#    else:
#        chi_prima = d*1j*np.sqrt(k_parallel**2 - epsilon2*k0**2)
#
#    if k_parallel**2 <= k0**2 :
#        chi = d*1j*np.sqrt(k0**2 - k_parallel**2)
#    else:
#        chi = d*np.sqrt(k_parallel**2 - k0**2)
#    
#    return 1/np.tan(chi_prima/2) + epsilon2*chi/chi_prima
#
#
def k_parallel_WG_TM1B_x(hbw,d,k_parallel):
    omegac = hbw/(hb*c)
    k0 = omegac
    
    epsilon2 =  epsilon_x(hbw)    
    if k_parallel**2 <= epsilon2*(k0**2):    
        chi_prima = d*np.sqrt(epsilon2*k0**2 - k_parallel**2)
    else:
        chi_prima = d*1j*np.sqrt(k_parallel**2 - epsilon2*k0**2)

    if k_parallel**2 <= k0**2 :
        chi = d*1j*np.sqrt(k0**2 - k_parallel**2)
    else:
        chi = d*np.sqrt(k_parallel**2 - k0**2)
    
    return np.tan(chi_prima/2) - epsilon2*chi/chi_prima
#
#def k_parallel_WG_TM2B_z(hbw,d,k_parallel):
#    omegac = hbw/(hb*c)
#    k0 = omegac
#    
#    epsilon2 =  epsilon_z(hbw)    
#    if k_parallel**2 <= epsilon2*(k0**2):    
#        chi_prima = d*np.sqrt(epsilon2*k0**2 - k_parallel**2)
#    else:
#        chi_prima = d*1j*np.sqrt(k_parallel**2 - epsilon2*k0**2)
#
#    if k_parallel**2 <= k0**2 :
#        chi = d*1j*np.sqrt(k0**2 - k_parallel**2)
#    else:
#        chi = d*np.sqrt(k_parallel**2 - k0**2)
#    
#    return np.tan(chi_prima/2) - epsilon2*chi/chi_prima

#%% ### modos s ###

def k_parallel_Silver_x(hbw,epsi1,epsi3,d_nano): 
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

    hb_omega_p = 170.1*1e-3 #eV

    d_micros = d_nano*1e-3
    
    A = d_micros*(hb_omega_p**2)/(epsi1 + epsi3)
    
    return hbw**2/A


def k_parallel_Silver_z(hbw,epsi1,epsi3,d_nano): 
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

    hb_omega_p = 92.5*1e-3 #eV

    d_micros = d_nano*1e-3
    
    A = d_micros*(hb_omega_p**2)/(epsi1 + epsi3)
    
    return hbw**2/A


def k_parallel_Silver_guess_x(hbw,epsi1,epsi3,d_nano): 
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

    omegac = hbw/(hb*c)
    d_micros = d_nano*1e-3


    epsilon2 =  epsilon_x(hbw)  
    
    r12 = (epsilon2*np.sqrt(epsi1) - epsi1*np.sqrt(epsilon2))/(epsilon2*np.sqrt(epsi1) + epsi1*np.sqrt(epsilon2))
    
    r23 = (epsi3*np.sqrt(epsilon2) - epsilon2*np.sqrt(epsi3))/(epsi3*np.sqrt(epsilon2) + epsilon2*np.sqrt(epsi3))
    
    
    A = -np.real(np.log(r12*r23)*4/(1j*d_micros*np.sqrt(epsilon2)))
    
    k_parallel = omegac/np.sqrt(A)
    
    
    k_parallel = np.sqrt(4/((d*epsilon2)**2) + omegac**2)
    
    return k_parallel



def k_parallel_Silver_guess_z(hbw,epsi1,epsi3,d_nano): 
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

    omegac = hbw/(hb*c)
    d_micros = d_nano*1e-3


    epsilon2 =  epsilon_z(hbw)  
    
    r12 = (epsilon2*np.sqrt(epsi1) - epsi1*np.sqrt(epsilon2))/(epsilon2*np.sqrt(epsi1) + epsi1*np.sqrt(epsilon2))
    
    r23 = (epsi3*np.sqrt(epsilon2) - epsilon2*np.sqrt(epsi3))/(epsi3*np.sqrt(epsilon2) + epsilon2*np.sqrt(epsi3))
    
    
    A = -np.real(np.log(r12*r23)*4/(1j*d_micros*np.sqrt(epsilon2)))
    
#    k_parallel = omegac/np.sqrt(A)
    
    
    k_parallel = np.sqrt(4/((d*epsilon2)**2) + omegac**2)
    
    return k_parallel

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

def graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad):
    plt.figure(figsize=tamfig)
    plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
    plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
    plt.tick_params(labelsize = tamnum, pad = pad)
    return   

#%%
    
ind = 0
d_nano_list = [1,10,20,100]
d = d_nano_list[ind]*1e-3  ## en micrones
epsi1, epsi3 = 1,1
d_nano = d*1e3

title = '$\epsilon_1 = %i$, $\epsilon_3 = %i$, d = %i nm, $\hbar\omega_{bulk}$ = 9.17 eV, $\hbar\gamma$ = 21 meV' %(epsi1,epsi3,d*1e3)

cota_inf_E = 0.001
cota_sup_E = 4
N = 800
list_E = np.linspace(cota_inf_E,cota_sup_E,N)

list_int_v1 = [1,1.25,1.5,1.75,2]
list_int_v2 = [1,2,3,4,5]
list_int_v3 = [1,2,5,10]
list_int_v4 = [1,10,25,35]

if d_nano == 100:
    list_int_v = list_int_v1
elif d_nano == 20:
    list_int_v = list_int_v2

elif d_nano == 10:
    list_int_v = list_int_v3

elif d_nano == 1:
    list_int_v = list_int_v4

#%%

list_y_air = []
list_y_medium_x = []
list_y_medium_z = []

list_y_Silver_x = []
list_y_Silver_z = []
list_y_Silver_aprox_x  = []
list_y_Silver_aprox_z  = []

for E in list_E:

    k_air = k_parallel_air(E)
    list_y_air.append(k_air)
    
    k_medium_x = k_parallel_medium_x(E)   
    list_y_medium_x.append(k_medium_x)
    
    k_silver_x = k_parallel_Silver_x(E,epsi1,epsi3,d_nano)
    list_y_Silver_x.append(k_silver_x)
    
    k_silver_approx_x = k_parallel_Silver_guess_x(E,epsi1,epsi3,d_nano)
    list_y_Silver_aprox_x.append(k_silver_approx_x)
    


    k_medium_z = k_parallel_medium_z(E)   
    list_y_medium_z.append(k_medium_z)
    
    k_silver_z = k_parallel_Silver_z(E,epsi1,epsi3,d_nano)
    list_y_Silver_z.append(k_silver_z)
    
    k_silver_approx_z = k_parallel_Silver_guess_z(E,epsi1,epsi3,d_nano)
    list_y_Silver_aprox_z.append(k_silver_approx_z)

    
    
#%%
        

#zeros_eqTM1B = []
#list_energy_TM1B = []
#
#k_TM1 = 0
#for E in list_E:
#
#    if k_TM1 == 0:
#        init_condit = [k_parallel_Silver(E,epsi1,epsi3,d)]
#    else:
#        init_condit = zeros_eqTM1B[k_TM1-1]
#        
#        
#    def k_parallel_WG_TM1_1var(k_parallel):   
#        return  np.abs(k_parallel_WG_TM1B(E,d,k_parallel))
#    
#
#    resTM1 = fsolve(k_parallel_WG_TM1_1var, init_condit,maxfev = 1000 )         
#    resTM1_v = np.float(resTM1)
#    if resTM1_v < 0 :
#        resTM1_v = -resTM1_v
#    
#    zeros_eqTM1B.append(resTM1_v)
#    list_energy_TM1B.append(E)
#    k_TM1 = k_TM1 + 1

zeros_eqTM2B_x = []
list_energy_TM2B_x = []
#list_omega_omegaWGB = np.linspace(2,cota_sup_omega_omegaWG,N)

k_TM2_x = 0
for E in list_E:

    if k_TM2_x == 0:
        init_condit = [k_parallel_Silver_x(E,epsi1,epsi3,d)]
    else:
        init_condit = zeros_eqTM2B_x[k_TM2_x-1]
            
    def k_parallel_WG_TM2_1var_x(k_parallel):   
        return  np.abs(k_parallel_WG_TM2B_x(E,d,k_parallel))
    

    resTM2_x = fsolve(k_parallel_WG_TM2_1var_x, init_condit,maxfev = 1000 )         
    resTM2_v = np.float(resTM2_x)

    if resTM2_v < 0 :
        resTM2_v = -resTM2_v

    zeros_eqTM2B_x.append(resTM2_v)
    list_energy_TM2B_x.append(E)
    k_TM2_x = k_TM2_x + 1
    
del resTM2_x,k_TM2_x


zeros_eqTM2B_z = []
list_energy_TM2B_z = []
#list_omega_omegaWGB = np.linspace(2,cota_sup_omega_omegaWG,N)

k_TM2_z = 0
for E in list_E:

    if k_TM2_z == 0:
        init_condit = [k_parallel_Silver_z(E,epsi1,epsi3,d)]
    else:
        init_condit = zeros_eqTM2B_x[k_TM2_z-1]
            
    def k_parallel_WG_TM2_1var_z(k_parallel):   
        return  np.abs(k_parallel_WG_TM1B_x(E,d,k_parallel))
    

    resTM2_z = fsolve(k_parallel_WG_TM2_1var_z, init_condit,maxfev = 1000 )         
    resTM2_v = np.float(resTM2_z)

    if resTM2_v < 0 :
        resTM2_v = -resTM2_v

    zeros_eqTM2B_z.append(resTM2_v)
    list_energy_TM2B_z.append(E)
    k_TM2_z = k_TM2_z + 1
    
del resTM2_z,k_TM2_z


#%%

#
plt.figure(figsize=tamfig)    
plt.title(title,fontsize=tamtitle)
#plt.plot(np.array(zeros_eqTM1B),list_energy_TM1B,'.',color = 'darkred',ms = 4)
plt.plot(np.array(zeros_eqTM2B_x),list_energy_TM2B_x,'.',color = 'blue',ms = 4,label = '$1-r_{21}r_{23}e^{2ik_{z,2}d} = 0$')
plt.plot(list_y_Silver_aprox_x,list_E,'.',color = 'darkred',ms = 4,label = 'analytical approx')
plt.xlabel('$k_\parallel$ [1/$\mu$m]',fontsize=tamletra, labelpad = 0)
plt.ylabel(r'$\hbar\omega$ [eV]',fontsize=tamletra, labelpad = 0)
#plt.plot(list_y_air,list_E,'--',color = 'grey',ms = lw, label = 'air')
#if d_nano != 100:
#    plt.plot(list_y_Silver_x,list_E,'--',color = 'black',ms = lw, label = '$\hbar\omega = \sqrt{\mathcal{A} k_\parallel}$')
list_colours = ['black', 'gold', 'orange', 'darksalmon', 'brown', 'darkred']
j = 0
maxis = []
for int_v in list_int_v:
    v = c/int_v
    
    list_y2_re = []
    list_y2_im = []
    for hb_omega in list_E:
        omega = hb_omega/hb        
        valuey_e = omega/v
        list_y2_re.append(valuey_e.real)
        list_y2_im.append(valuey_e.imag)  
    maxis.append(np.max(list_y2_re))
    ##  
    if int_v != 1:
        if d_nano != 100:
            plt.plot(list_y2_re,list_E,'-', color = list_colours[j],ms = lw, label = 'v= c/%i' %(int_v))
        else:
            plt.plot(list_y2_re,list_E,'-', color = list_colours[j],ms = lw, label = 'v= c/%.2f' %(int_v))
    else:
        plt.plot(list_y2_re,list_E,'-', color = list_colours[j], ms = lw)
        
    j = j + 1

if d_nano == 1:
    min1 = -15
elif d_nano == 10: 
    min1 = -5
elif d_nano == 20: 
    min1 = -4
elif d_nano == 100: 
    min1 = -0.5


#plt.xlim([min1,np.max([np.max(list_y_Silver_x), np.max(maxis)])+0.5])
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
plt.grid(1)
plt.tight_layout()
if save_plots == 1 :
    plt.savefig('disp_relation_hBn_x_real_vs_Ev_d%inm.png' %(d*1e3))




#%%


#
plt.figure(figsize=tamfig)    
plt.title(title,fontsize=tamtitle)
#plt.plot(np.array(zeros_eqTM1B),list_energy_TM1B,'.',color = 'darkred',ms = 4)
#plt.plot(np.array(zeros_eqTM2B_z),list_energy_TM2B_z,'.',color = 'blue',ms = 4,label = '$1-r_{21}r_{23}e^{2ik_{z,2}d} = 0$')
plt.plot(list_y_Silver_aprox_z,list_E,'.',color = 'darkred',ms = 4,label = 'analytical approx')
plt.xlabel('$k_\parallel$ [1/$\mu$m]',fontsize=tamletra, labelpad = 0)
plt.ylabel(r'$\hbar\omega$ [eV]',fontsize=tamletra, labelpad = 0)
#plt.plot(list_y_air,list_E,'--',color = 'grey',ms = lw, label = 'air')
#if d_nano != 100:
#    plt.plot(list_y_Silver_z,list_E,'--',color = 'black',ms = lw, label = '$\hbar\omega = \sqrt{\mathcal{A} k_\parallel}$')
list_colours = ['black', 'gold', 'orange', 'darksalmon', 'brown', 'darkred']
j = 0
maxis = []
for int_v in list_int_v:
    v = c/int_v
    
    list_y2_re = []
    list_y2_im = []
    for hb_omega in list_E:
        omega = hb_omega/hb        
        valuey_e = omega/v
        list_y2_re.append(valuey_e.real)
        list_y2_im.append(valuey_e.imag)  
    maxis.append(np.max(list_y2_re))
    ##  
    if int_v != 1:
        if d_nano != 100:
            plt.plot(list_y2_re,list_E,'-', color = list_colours[j],ms = lw, label = 'v= c/%i' %(int_v))
        else:
            plt.plot(list_y2_re,list_E,'-', color = list_colours[j],ms = lw, label = 'v= c/%.2f' %(int_v))
    else:
        plt.plot(list_y2_re,list_E,'-', color = list_colours[j], ms = lw)
        
    j = j + 1

if d_nano == 1:
    min1 = -15
elif d_nano == 10: 
    min1 = -5
elif d_nano == 20: 
    min1 = -4
elif d_nano == 100: 
    min1 = -0.5


#plt.xlim([min1,np.max([np.max(list_y_Silver_z), np.max(maxis)])+0.5])
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
plt.grid(1)
plt.tight_layout()
if save_plots == 1 :
    plt.savefig('disp_relation_hBn_z_real_vs_Ev_d%inm.png' %(d*1e3))




#%%

if save_txt == 1:
    
    header1 = 'Energy [eV]     k_parallel [1/micrones]'  + ', '  + title + ', ' + name_this_py
     
#    tabla_TM1B = np.array([list_energy_TM1B_x,zeros_eqTM1B])
#    tabla_TM1B = np.transpose(tabla_TM1B)
#    
    tabla_TM2B = np.array([list_energy_TM2B_x,zeros_eqTM2B_x])
    tabla_TM2B = np.transpose(tabla_TM2B)


#    np.savetxt('sol_TM1B_d%inm.txt' %(d*1e3), tabla_TM1B, fmt='%1.11e', delimiter='\t', header = header1 )
    np.savetxt('sol_TM2B_x_d%inm.txt' %(d*1e3), tabla_TM2B, fmt='%1.11e', delimiter='\t', header = header1 )


    tabla_TM2B = np.array([list_energy_TM2B_z,zeros_eqTM2B_z])
    tabla_TM2B = np.transpose(tabla_TM2B)


#    np.savetxt('sol_TM1B_d%inm.txt' %(d*1e3), tabla_TM1B, fmt='%1.11e', delimiter='\t', header = header1 )
    np.savetxt('sol_TM2B_z_d%inm.txt' %(d*1e3), tabla_TM2B, fmt='%1.11e', delimiter='\t', header = header1 )






#%%

