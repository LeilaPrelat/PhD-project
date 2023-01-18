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

graficar = 1

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_graphene =  path_basic.replace('/' + 'plane','')
#print('Importar modulos necesarios para este codigo')

try:
    sys.path.insert(1, path_basic)
    from constants import constantes
except ModuleNotFoundError:
    print('constants_plane.py no se encuentra en ' + path_basic)

pi,hb,c,alfac,mu1,mu2 = constantes()
aux = c*hb

#%%
    
def k_parallel(omega,omega_bulk,gamma_in,d,epsilon1,epsilon2,epsilon_b): 
    """
    Parameters
    ----------
    omega : frecuencia en Hz
    omega_bulk : frecuencia de bulk Hz
    gamma_in : frecuencia de colision en eV
    d : espesor del plano Ag en micrometros
    epsilon1 : permeabilidad electrica del medio de arriba
    epsilon2 : permeabilidad electrica del medio de abajo
    epsilon_b : permeabilidad electrica para Ag (aprox de DL)
    Returns
        relacion de dispersion para un plano de Ag
        (solucion analitica)
    -------
    """
    cte = (epsilon1 + epsilon2)/d
    
    den = omega*(omega + 1j*gamma_in/hb)
    f  = 1 - epsilon_b + omega_bulk**2/den
    
    return cte/f         

#%%

if graficar == 1:
    
    import matplotlib.pyplot as plt
    import seaborn as sns
    import numpy as np
    from scipy.optimize import minimize   
    
    sns.set()
    path_save = path_basic + '/' + 'disp_relation_Ag'
    colors = ['darkred','steelblue','coral','yellowgreen']
    
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
    
    epsilon1,epsilon2 = 12,2.4
    
    gamma_in = 0.0001      # collision frequency in eV
    Ebulk = 9.17
    omega_bulk = Ebulk/hb
    d = 0.0001 #micrones, 10 nm
    int_v = 300
    v = c/int_v
    epsilonb = 4
    
    n  = 150
    
    
    title1 = '$\epsilon_1 = %.1f$, $\epsilon_2 = %.1f$, $\epsilon_b = %i$, $\gamma_{in}$ = %.4feV' %(epsilon1,epsilon2,epsilonb,gamma_in)
    title2 = '$E_{bulk}$ = %.2feV, d = %.4f$\mu$m' %(Ebulk,d)

    def k_parallel_ev(omega):
        return k_parallel(omega,omega_bulk,gamma_in,d,epsilon1,epsilon2,epsilonb)
    
    list_E = np.linspace(0.01,1,n)
    list_omega = []
    list_y_re = []
    list_y_im = []
    for E in list_E:
        omega = E/hb
        valuey = k_parallel_ev(omega)
        list_omega.append(omega)
        list_y_re.append(valuey.real)
        list_y_im.append(valuey.imag)


    recta = []
    recta_luz = []
    for kparallel in list_y_re:
        rta = kparallel*v*hb
        recta.append(rta)
        rta2 = kparallel*c*hb
        recta_luz.append(rta2)
        
    plt.figure(figsize=tamfig)    
    plt.title(title1 + '\n' + title2,fontsize=tamtitle)
    plt.plot(list_y_re,list_E,'-',color = colors[0],ms = lw, label = 'SP')
    plt.plot(list_y_re,recta,'-',color = colors[1],ms = lw, label = 'v= c/%i' %(int_v))
    #plt.plot(list_y_re,recta_luz,'-',color = colors[2],ms = lw)
    plt.xlabel('Re($k_\parallel$) [1/$\mu$m]',fontsize=tamletra, labelpad = 0)
    plt.ylabel('energy [eV]',fontsize=tamletra, labelpad = 0)
    plt.tick_params(labelsize = tamnum,pad = pad)
    plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
    plt.grid(1)
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig('disp_relationAg_real.png')

    plt.figure(figsize=tamfig)      
    plt.title(title1 + '\n' + title2,fontsize=tamtitle)
    plt.plot(list_y_im,list_E,'-',color = colors[0],ms = lw)
    plt.xlabel('Im($k_\parallel$) [1/$\mu$m]',fontsize=tamletra, labelpad = 0)
    plt.ylabel('energy [eV]',fontsize=tamletra, labelpad = 0)
    plt.tick_params(labelsize = tamnum,pad = pad)
    plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
    plt.grid(1)
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig('disp_relationAg_imag.png')

##

    def electron(v,omega):
        return omega/v
    
    def resta_function(omega):
        elec = electron(v,omega)
        sp = k_parallel_ev(omega)
        return hb*np.abs(elec - sp)

    res = minimize(resta_function, 0.75*1e15, method='Nelder-Mead', tol=1e-13, 
                   options={'maxiter':1150})
    
    omega_c = res.x[0]*1e-15
    
    list_y2_re = []
    list_y2_im = []
    for E in list_E:
        omega = E/hb
        valuey_e = electron(v,omega)
        valuey_e = valuey_e*hb
        list_y2_re.append(valuey_e.real)
        list_y2_im.append(valuey_e.imag)    
    
    
    plt.figure(figsize=tamfig)    
    plt.title(title1 + '\n' + title2 + '\n' + 'v= c/%i' %(int_v) + ', $\omega_c$ = %.3f 10$^{15}$Hz'%(omega_c),fontsize=tamtitle)
    plt.plot(list_y_re,list_omega,'-',color = colors[0],ms = lw, label = 'Ag')
    plt.plot(np.array(list_omega)/v,list_omega,'-',color = colors[1],ms = lw, label = r'$\omega = v k_\parallel$')
    # plt.plot(np.array(list_omega)/c,list_omega,'-',color = colors[2],ms = lw, label = 'line light')    
    plt.xlabel('Re($k_\parallel$) [1/$\mu$m]',fontsize=tamletra, labelpad = 0.5)
    plt.ylabel('$\omega$ [Hz]',fontsize=tamletra, labelpad = 0.5)
    plt.tick_params(labelsize = tamnum,pad = pad)
    plt.legend(loc='best',markerscale=2,fontsize=int(tamlegend*0.8))
    plt.grid(1)
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig('disp_relationAg_real2.png')

    plt.figure(figsize=tamfig)      
    plt.title(title1 + '\n' + title2,fontsize=tamtitle)
    plt.plot(list_y_im,list_omega,'-',color = colors[0],ms = lw)
    # plt.plot(list_y2_im,list_omega,'-',color = colors[1],ms = lw)
    plt.xlabel('Im($k_\parallel$) [1/$\mu$m]',fontsize=tamletra, labelpad = 0.5)
    plt.ylabel('$\omega$ [Hz]',fontsize=tamletra, labelpad = 0.5)
    plt.tick_params(labelsize = tamnum,pad = pad)
    plt.grid(1)
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig('disp_relationAg_imag2.png')

#%%

