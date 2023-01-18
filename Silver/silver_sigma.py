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
import sys
import os 

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')

graficar = 0

try:
    sys.path.insert(1, path_basic)
    from Silver_PP import Silver_lambda_p
except ModuleNotFoundError:
    print('graphene_sigma.py no se encuentra en ' + path_basic)

#%%

def sigma_DL(hbw,epsi1,epsi3,d_nano):
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

    d_micros = d_nano*1e-3
    lambda_p = Silver_lambda_p(hbw,epsi1,epsi3)*d_micros
    
#    c = 3*10**(14)                ### light velocity in micron/seg
#    alfac = 1/137.0359            ### fine structure
    

    num1 = 1j*(epsi1 + epsi3)/(4*np.pi*alfac)

    num2 = np.real(lambda_p)/(2*np.pi)
    
    return num1*num2

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
        sigma_real = sigma(hbar_omega,epsi1,epsi3,d_nano).real            
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
        sigma_imag = sigma(hbar_omega,epsi1,epsi3,d_nano).imag
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
