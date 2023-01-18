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


name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')



err = 'decay_rate_film3.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from polarizability_disk import polarizability
except ModuleNotFoundError:
    print(err)

try:
    sys.path.insert(1, path_basic)
    from constants import constantes
except ModuleNotFoundError:
    print('constants.py no se encuentra en ' + path_basic)

pi,hb,c,alfac,mu1,mu2 = constantes()


def cross_sec(hbw,hbmu,hbgama,D,epsilon):

    alpha = polarizability(hbw,hbmu,hbgama,D,epsilon)
    omegac = hbw/(c*hb)

    cross_section = 4*np.pi*omegac*np.imag(alpha)/np.sqrt(epsilon)

    return cross_section

    
graficar = 1
    
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
        
    tamfig = (4.5,3.5)
    tamlegend = 10
    tamletra = 10
    tamtitle = 10
    tamnum = 9
    loc2 = [0,1]
    pad = -2
    lw = 1.5
     
    hbargama = 0.0001      # collision frequency in eV
    hbmu = 0.35
    epsilon = 1
    
    #
    
    list_E = np.linspace(0.05,0.2)
    D_nano = 120
    D = D_nano*1e-3
    
    listy = cross_sec(list_E,hbmu,hbargama,D,epsilon)

    title = '$\hbar\mu$ = %.2f eV, $\hbar\gamma$ = %.2f meV, D = %i nm, $\epsilon$ = %i' %(hbmu,hbargama*1e3,D_nano,epsilon)

    plt.figure(figsize=tamfig)
    plt.plot(list_E,listy,'-',color = 'darkred',ms = lw)

    plt.title(title,fontsize=tamtitle)
    plt.ylabel('Cross section',fontsize=tamletra, labelpad = -1)
    plt.xlabel('Photon Energy [eV]',fontsize=tamletra, labelpad = -1)
    
 #   plt.yscale('log')
    plt.tick_params(labelsize = tamnum,pad = pad)
    plt.legend(loc = loc2, ncol = 2,markerscale=1,fontsize=tamlegend,frameon=False,handletextpad=0.05)
    plt.grid(1)
    plt.tight_layout()
    plt.savefig('cross_sections.png')
    
#%%
