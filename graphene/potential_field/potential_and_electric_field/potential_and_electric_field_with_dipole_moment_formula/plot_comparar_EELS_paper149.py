
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: leila

graficar el campo externo directo con la convencion de z hacia abajo
en z = 0
graficar mapa de color x,y
"""

import numpy as np
import sys
import os 
import matplotlib.pyplot as plt
from scipy import special


#%%
plot_vs_E = 1
plot_vs_c = 0
plot_vs_zp = 0

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/potential_field/potential_and_electric_field/potential_and_electric_field_with_dipole_moment_formula','')
path_save = path_basic + '/' + 'decay_rate_film'
if not os.path.exists(path_save):
    print('Creating folder to save graphs')
    os.mkdir(path_save)


try:
    sys.path.insert(1, path_basic)
    from dipole_moment import dipole_moment_anav2_for_decay_rate_resonance_dir
except ModuleNotFoundError:
    print('dipole_moment.py no se encuentra en ' + path_basic)

try:
    sys.path.insert(1, path_constants)
    from constants import constantes
except ModuleNotFoundError:
    print('constants.py no se encuentra en ' + path_constants)

pi,hb,c,alfac,mu1,mu2 = constantes()
aux = c*hb
aux2 = 1e12/c

#%%

print('Definir parametros del problema')

#v = c/int_v
#omega = 0.7*1e12

#v = c/int_v
#omega = 0.7*1e12

#v = c/int_v

epsi1,epsi2 = 1,1
hbmu,hbgama = 0.3,0.0001
zp = 0.05
b = -0.01


 
#title1 = r'$\kappa$ = %.2f$\omega_0$, $\kappa_r$ = %.2f$\kappa$, $E_0$=%i meV' %(kappa_factor_omega0, kappa_r_factor, energy0_pol)     
#title2 = r'$\hbar\mu$ = %.2feV, $\hbar\gamma$ = %.4feV' %(hbmu,hbgama) 
#title3 = r'$z_p$=%inm, px=%i, py=%i, pz=%i' %(zp*1e3,px,py,pz)
title4 = r'b = %i nm' %(b*1e3)
labelp = r'_res' 

N = 200

def function_imag_ana(energy0,int_v,zp_nano):
    omegac0 = energy0*1e-3/aux 
    zp = zp_nano*1e-3

    cte = 24*np.pi*(omegac0**2)


    arg = np.abs(b)*omegac0*int_v
    K1 = special.kn(1,arg)
    K0 = special.kn(0,arg)
    
    rta_K = K0**2 + K1**2
    px_dir,py_dir,pz_dir = dipole_moment_anav2_for_decay_rate_resonance_dir(omegac0,int_v,b,zp)
    denominador = np.abs(px_dir)**2 +  np.abs(py_dir)**2 +  np.abs(pz_dir)**2
    
    rta1 = rta_K*cte/denominador
    
#    rta2 = EELS_dir_ana_f(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,b,zp)
    
    return rta1

if plot_vs_c == 1 :
    E0 = 44 # meV
    # z0 = 0.06*1e3
    zp0 = 0.05*1e3 
    
    labelx = r'v/c'   
    title4 = title4 + ', ' + r'$\hbar\omega$ = %i meV, $z_p$ = %i nm' %(E0,zp0)
    label1 = 'vs_v' + labelp
#    listx = np.linspace(0.0001,2,N)
    listx = np.linspace(0.005,0.12,N)

if plot_vs_E ==1 :
    # z0 = 0.06*1e3
    int_v0 = 10
    zp0 = 0.05*1e3 
    
    labelx = r'$\hbar\omega$ [meV]'   
    title4 = title4 + ', ' + r'v = c/%i, $z_p$ = %i nm' %(int_v0,zp0)
    label1 = 'vs_E' + labelp
#    listx = np.linspace(0.0001,2,N)
    listx = np.linspace(15,65,N)

if plot_vs_zp == 1 : 
    E0 = 45 # meV
    int_v0 = 10

    labelx = r'Surface-dipole distance, $z_{\rm 0}$/$\lambda_{\rm p}$'   
    title4 = title4 + ', ' + r'v = c/%i, $\hbar\omega$ = %i meV' %(int_v0,E0)
    label1 = 'vs_zp' + labelp + '_E%imeV' %(E0)
#    listx = np.linspace(0.0001,2,N)
    listx = np.linspace(10,4000,N)
    

title =  title4 


#%%
    
tamfig = [2.5, 2]
tamletra = 9
tamtitle  = 8
tamnum = 7
tamlegend = 8
labelpady = 2
labelpadx = 3
pad = 3
mk = 1
ms = 2
hp = 0.3
length_marker = 0
dpi = 500

def graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad):
    plt.figure(figsize=tamfig)
    plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
    plt.ylabel(labely,fontsize=tamletra,labelpad =labelpadx)
    plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in", pad = pad)
#    plt.title(title,fontsize=int(tamtitle*0.9))

    return   


#%%


listy_im_ana = []
list_ana_parallel = []


if plot_vs_E == 1: 

    for value in listx: 

        y_im_ana = function_imag_ana(value,int_v0,zp0)        
        listy_im_ana.append(y_im_ana)
        


elif plot_vs_c == 1:       

    for value in listx: 
        
        value2 = 1/value
#        print(value2)

        y_im_ana = function_imag_ana(E0,value2,zp0)        
        listy_im_ana.append(y_im_ana)



elif plot_vs_zp == 1:
    
    for value in listx: 

        y_im_ana = function_imag_ana(E0,int_v0,value)        
        listy_im_ana.append(y_im_ana)
        
     
    
#%%??

graph(title,labelx,r'$\Gamma^{\#149}_{0}/\Gamma_{\rm EELS}$',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#plt.title(title,fontsize=int(tamtitle*0.9))
plt.plot(listx,np.array(listy_im_ana),'-',ms = ms,color = 'purple')
#plt.plot(listx,list_ana_parallel,'.-',ms = ms,color = 'darkred',label = r'$\Gamma_{\parallel}$')
#plt.plot(np.ones(10)*zp_crit_lambda_p_value, np.array(listy_aux)*1e-4,'--k')

plt.tight_layout()
#if plot_vs_c == 1:
#    plt.yscale('log')
os.chdir(path_save)
plt.savefig( 'Gamma_compared_' + label1 + '.png', format='png',bbox_inches='tight',pad_inches = 0.01, dpi=dpi)   
