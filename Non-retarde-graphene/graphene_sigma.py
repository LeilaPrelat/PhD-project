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
graficar = 0

#%%

def sigma(hbw,hbmu,hbgama):
    """
    Parameters
    ----------
    hbw : energia = hbar*omega en eV
    hbmu : potencial quimico del grafeno mu (en eV)
    hbgama : collision frequency in eV 

    Returns
    -------
    devuelve la conductividad del grafeno
    en unidades de e^2/hbar 
	---> multiplicar despues por fine_structure*c (unidades gaussianas)

    """
    global Tk  ### variables globales ---> no pueden cambiarse
    akb = 8.6173324E-5           ### Boltzman constant in eV/K  
    Tk = 300 			### temperature in K ---> no cambia
    
    # TKmu = Tk*akb/hbmu
    
    # if TKmu==0:
    #     # intra = 1j*(1/np.pi)*abs(hbmu)/(hbw+1j*hbgama)
    #     # if hbw-2*abs(hbmu)>0:
    #     #     inter1 = 0.25
    #     # else: 
    #     #     inter1 = 0
         
    #     # inter = inter1 + 1j*np.log((hbw-2*abs(hbmu))**2/(hbw+2*abs(hbmu))**2)/(4*np.pi)
    #     raise ValueError('agregar formula para T = 0K en graphene_sigma.py')
    # else:

    aux2 = hbw + 1j*hbgama
    aux3 = (aux2 - 2*hbmu)**2 + (2*Tk*akb)**2
    aux4 = (aux2 + 2*hbmu)**2 
    aux5 = (aux2 - 2*hbmu)/(2*akb*Tk)   
    inter = 0.25*(0.5 + np.arctan(aux5)/np.pi - 1j*np.log(aux4/aux3)/(2*np.pi))   

    auxTkmu = hbmu/(2*Tk*akb)
    aux = np.log(np.exp(auxTkmu) + np.exp(-auxTkmu))
    intra_aux = aux/aux2
    intra = 2j*akb*Tk*intra_aux/np.pi  

    sigmatot = inter + intra
    return sigmatot  

#%%

def sigma_DL(hbw,hbmu,hbgama):
    """
    Parameters
    ----------
    hbw : energia = hbar*omega en eV
    hbmu : potencial quimico del grafeno mu (en eV)
    hbgama : collision frequency in eV 

    Returns
    -------
    devuelve la conductividad del grafeno
    del modelo drude lorentz (parte intra aproximada)
    en unidades de e^2/hbar 
	---> multiplicar despues por fine_structure*c (unidades gaussianas)
    

    """
    
#    c = 3*10**(14)                ### light velocity in micron/seg
#    alfac = 1/137.0359            ### fine structure
    
    den = hbw + 1j*hbgama
    num = 1j*hbmu/np.pi

    return num/den

#%%

if graficar == 1:

    import matplotlib.pyplot as plt
    import sys
    import seaborn as sns
    import os
    sns.set()
    name_this_py = os.path.basename(__file__)
    path = os.path.abspath(__file__) #path absoluto del .py actual
    path_basic = path.replace('/' + name_this_py,'')
    path_save = path_basic + '/' + 'sigma_graphene'
    
    try:
        sys.path.insert(1, path_basic)
        from constants import constantes
    except ModuleNotFoundError:
        print('constants.py no se encuentra en ' + path_basic)
    
    pi,hb,c,alfac,mu1,mu2 = constantes()
    
    aux = hb*c
    tamfig = (4.5,3.5)
    tamlegend = 10
    tamletra = 10
    tamtitle = 10
    tamnum = 9
    loc2 = [0,1]
    pad = -2
    lw = 1.5
     
    hbargama = 0.0001      # collision frequency in eV
    list_mu = np.linspace(0.3,0.9,601)
    colors = ['darkred','steelblue','coral','yellowgreen']
    list_hbar_omega = [0.1,0.3,0.6,0.9]
    
    plt.figure(figsize=tamfig)
    k = 0
    for hbar_omega in list_hbar_omega: 
        sigma_real = sigma(hbar_omega,list_mu,hbargama).real            
        plt.plot(list_mu,sigma_real,'-',color = colors[k],ms = lw,label = '$\hbar \omega$ = %.2f eV' %(hbar_omega))
        k = k + 1
    #plt.title('$\sigma$ con $\hbar \omega$ = %.2f eV, $\gamma$ = %.4f eV' %(hbar_omega,hbargama),fontsize=tamtitle)
    plt.ylabel('Re($\sigma$)',fontsize=tamletra, labelpad = -1)
    plt.xlabel('$\mu_c$ [eV]',fontsize=tamletra, labelpad = -1)
    
 #   plt.yscale('log')
    plt.tick_params(labelsize = tamnum,pad = pad)
    plt.legend(loc = loc2, ncol = 2,markerscale=1,fontsize=tamlegend,frameon=False,handletextpad=0.05)
    plt.grid(1)
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig('sigma_real.png')

    plt.figure(figsize=tamfig)
    k = 0
    for hbar_omega in list_hbar_omega:          
        sigma_imag = sigma(hbar_omega,list_mu,hbargama).imag
        plt.plot(list_mu,sigma_imag,'-',color = colors[k],ms = lw,label = '$\hbar \omega$ = %.2f eV' %(hbar_omega))
        k = k + 1
    #plt.title('$\sigma$ con $\hbar \omega$ = %.2f eV, $\gamma$ = %.4f eV' %(hbar_omega,hbargama),fontsize=tamtitle)
    plt.ylabel('Im($\sigma$)',fontsize=tamletra, labelpad = -1)
    plt.xlabel('$\mu_c$ [eV]',fontsize=tamletra, labelpad = -1)
    
 #   plt.yscale('log')
    plt.tick_params(labelsize = tamnum,pad = pad)
    plt.legend(loc = loc2, ncol = 2,markerscale=1,fontsize=tamlegend,frameon=False,handletextpad=0.05)
    plt.tight_layout()
    plt.grid(1)
    plt.savefig('sigma_imag.png')

    #
    plt.figure(figsize=tamfig)
    k = 0
    for hbar_omega in list_hbar_omega: 
        sigma_real = sigma_DL(hbar_omega,list_mu,hbargama).real            
        plt.plot(list_mu,sigma_real,'-',color = colors[k],ms = lw,label = '$\hbar \omega$ = %.2f eV' %(hbar_omega))
        k = k + 1
    #plt.title('$\sigma$ con $\hbar \omega$ = %.2f eV, $\gamma$ = %.4f eV' %(hbar_omega,hbargama),fontsize=tamtitle)
    plt.ylabel('Re($\sigma$)',fontsize=tamletra, labelpad = -1)
    plt.xlabel('$\mu_c$ [eV]',fontsize=tamletra, labelpad = -1)
    
 #   plt.yscale('log')
    plt.tick_params(labelsize = tamnum,pad = pad)
    plt.legend(loc = loc2, ncol = 2,markerscale=1,fontsize=tamlegend,frameon=False,handletextpad=0.05)
    plt.grid(1)
    plt.tight_layout()
    plt.savefig('sigma_realDL.png')

    plt.figure(figsize=tamfig)
    k = 0
    for hbar_omega in list_hbar_omega:          
        sigma_imag = sigma_DL(hbar_omega,list_mu,hbargama).imag
        plt.plot(list_mu,sigma_imag,'-',color = colors[k],ms = lw,label = '$\hbar \omega$ = %.2f eV' %(hbar_omega))
        k = k + 1
    #plt.title('$\sigma$ con $\hbar \omega$ = %.2f eV, $\gamma$ = %.4f eV' %(hbar_omega,hbargama),fontsize=tamtitle)
    plt.ylabel('Im($\sigma$)',fontsize=tamletra, labelpad = -1)
    plt.xlabel('$\mu_c$ [eV]',fontsize=tamletra, labelpad = -1)
    
 #   plt.yscale('log')
    plt.tick_params(labelsize = tamnum,pad = pad)
    plt.legend(loc = loc2, ncol = 2,markerscale=1,fontsize=tamlegend,frameon=False,handletextpad=0.05)
    plt.grid(1)
    plt.tight_layout()
    plt.savefig('sigma_imagDL.png')

    
#%%
