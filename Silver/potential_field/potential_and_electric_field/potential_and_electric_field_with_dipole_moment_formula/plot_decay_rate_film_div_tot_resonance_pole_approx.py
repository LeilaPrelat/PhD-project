
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

#%%
plot_vs_E = 0
plot_vs_c = 0
plot_vs_zp = 1

plot_color_maps = 1
paper = 0


if plot_color_maps == 0:
    import seaborn as sns
    sns.set()


#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/potential_field/potential_and_electric_field/potential_and_electric_field_with_dipole_moment_formula','')
path_save = path_basic + '/' + 'decay_rate_film_div_tot_pole_approx'
if not os.path.exists(path_save):
    print('Creating folder to save graphs')
    os.mkdir(path_save)

err = 'decay_rate_film3.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from decay_rate_film_resonance_div_tot_pole_approx import EELS_film_pole_aprox_f, EELS_dir_pole_approx_f
except ModuleNotFoundError:
    print(err)
    
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
#omega = 0.7*1e12

epsi1,epsi3 = 1,1

zp = 0.05
b = -0.01

d_nano = 100

#energy0_pol = 43
#omega0 = energy0_pol*1e-3/hb 
##R = 10 # 10 nm en unidades de micrometros
#kappa_factor_omega0 = 0.1
#kappa_r_factor= 0.5
 
#title1 = r'$\kappa$ = %.2f$\omega_0$, $\kappa_r$ = %.2f$\kappa$, $E_0$=%i meV' %(kappa_factor_omega0, kappa_r_factor, energy0_pol)     
#title2 = r'$\hbar\mu$ = %.2feV, $\hbar\gamma$ = %.4feV' %(hbmu,hbgama) 
#title3 = r'$z_p$=%inm, px=%i, py=%i, pz=%i' %(zp*1e3,px,py,pz)
title4 = r'b = %i nm, d = %i nm' %(b*1e3,d_nano)
labelp = r'_res' 

N = 100


#%%

def function_imag_ana(energy0,int_v,zp_nano):
    omegac0 = energy0/aux 
    zp = zp_nano*1e-3

    rta = EELS_film_pole_aprox_f(omegac0,epsi1,epsi3,d_nano,int_v,b,zp)
    
    return rta



def function_imag_ana_epsilon(energy0,int_v,zp_nano):
    omegac0 = energy0/aux 
    zp = zp_nano*1e-3

    rta = EELS_dir_pole_approx_f(omegac0,epsi1,epsi3,d_nano,int_v,b,zp)
    
    return rta

if plot_color_maps == 1:
    
    if d_nano == 100:
        int_v0 = 20
    elif d_nano == 1: 
        int_v0 = 400
    
    title4 = title4 + ', ' + r'v = c/%i' %(int_v0)

    def function_imag_ana_par_3D(energy0,zp_nano):
        omegac0 = energy0/aux 
        zp = zp_nano*1e-3

        rta_par, rta_perp = EELS_film_pole_aprox_f(omegac0,epsi1,epsi3,d_nano,int_v0,b,zp)
    
        return np.log10(rta_par)


    def function_imag_ana_perp_3D(energy0,zp_nano):
        omegac0 = energy0/aux 
        zp = zp_nano*1e-3

        rta_par, rta_perp = EELS_film_pole_aprox_f(omegac0,epsi1,epsi3,d_nano,int_v0,b,zp)
    
        return np.log10(rta_perp)


#
#    def function_imag_ana_epsilon_3D(energy0,zp_nano):
#        omegac0 = energy0/aux 
#        zp = zp_nano*1e-3
#
#        rta = EELS_dir_pole_approx_f(omegac0,epsi1,epsi3,d_nano,int_v0,b,zp)
#    
#        return rta
    
#%%
    

if plot_vs_c == 1 :
    E0 = 3 # meV
    # z0 = 0.06*1e3
    zp0 = 0.05*1e3 
    
    labelx = r'v/c'   
    title4 = title4 + ', ' + r'$\hbar\omega$ = %i eV, $z_p$ = %i nm' %(E0,zp0)
    label1 = 'vs_v' + labelp
#    listx = np.linspace(0.0001,2,N)
    listx = np.linspace(0.005,0.12,N)

if plot_vs_E ==1 :
    # z0 = 0.06*1e3
    int_v0 = 10
    zp0 = 0.05*1e3 
    
    labelx = r'$\hbar\omega$ [eV]'   
    title4 = title4 + ', ' + r'v = c/%i, $z_p$ = %i nm' %(int_v0,zp0)
    label1 = 'vs_E' + labelp
#    listx = np.linspace(0.0001,2,N)
    listx = np.linspace(0.1,3.15,N)

if plot_vs_zp == 1 : 
    E0 = 3 # eV
    if d_nano == 100:
        int_v0 = 20
    elif d_nano == 1: 
        int_v0 = 400


    labelx = r'$z_p$ [nm]'   
    title4 = title4 + ', ' + r'v = c/%i, $\hbar\omega$ = %.2f eV' %(int_v0,E0)
    label1 = 'E%.2feV_d%inm_vs_zp' %(E0,d_nano) + labelp
#    listx = np.linspace(0.0001,2,N)
    listx = np.linspace(0.4,50,N)
    
title =  title4 

#%%
    
tamfig = (4.5,3.5)
tamlegend = 12
tamletra = 11
tamtitle = 10
tamnum = 9
labelpady = -1.5
labelpadx = 2
pad = 0
mk = 2
ms = 4
hp = 0.3
length_marker = 1

def graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad):
    plt.figure(figsize=tamfig)
    plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
    plt.ylabel(labely,fontsize=tamletra,labelpad =labelpadx)
    plt.tick_params(labelsize = tamnum, pad = pad)
    plt.title(title,fontsize=int(tamtitle*0.9))

    return   

#%%

if plot_color_maps == 0:
    listy_im_ana_par = []
    listy_im_ana_perp = []
    
    listy_im_ana_par_2 = []
    listy_im_ana_perp_2 = []
    
    if plot_vs_E == 1: 
    
        for value in listx: 
    
            y_im_ana_par, y_im_ana_perp = function_imag_ana(value,int_v0,zp0)        
            y_im_ana_epsilon_par, y_im_ana_epsilon_perp = function_imag_ana_epsilon(value,int_v0,zp0)    
    
            listy_im_ana_par_2.append(y_im_ana_par/(y_im_ana_epsilon_par + y_im_ana_par))
            listy_im_ana_perp_2.append(y_im_ana_perp/(y_im_ana_epsilon_perp + y_im_ana_perp))
            
            listy_im_ana_par.append(y_im_ana_par)
            listy_im_ana_perp.append(y_im_ana_perp)   
                   
    elif plot_vs_c == 1:       
    
        
        for value in listx: 
            
            value2 = 1/value
    #        print(value2)
    
            y_im_ana_par, y_im_ana_perp = function_imag_ana(E0,value2,zp0)        
            y_im_ana_epsilon_par, y_im_ana_epsilon_perp  = function_imag_ana_epsilon(E0,value2,zp0)    
    
            listy_im_ana_par_2.append(y_im_ana_par/(y_im_ana_epsilon_par + y_im_ana_par))
            listy_im_ana_perp_2.append(y_im_ana_perp/(y_im_ana_epsilon_perp + y_im_ana_perp))
    
            listy_im_ana_par.append(y_im_ana_par)
            listy_im_ana_perp.append(y_im_ana_perp)   
            
    elif plot_vs_zp == 1:
    
        
        for value in listx: 
    
            y_im_ana_par, y_im_ana_perp  = function_imag_ana(E0,int_v0,value)        
            y_im_ana_epsilon_par, y_im_ana_epsilon_perp  = function_imag_ana_epsilon(E0,int_v0,value)  
    
            listy_im_ana_par_2.append(y_im_ana_par/(y_im_ana_epsilon_par + y_im_ana_par))
            listy_im_ana_perp_2.append(y_im_ana_perp/(y_im_ana_epsilon_perp + y_im_ana_perp))
            
            listy_im_ana_par.append(y_im_ana_par)
            listy_im_ana_perp.append(y_im_ana_perp)   
    
    
    graph(title + '\n' + '$\mathbf{p}$ $\parallel$ $\hat{x}$' ,labelx,'$\Gamma_{SP}/\Gamma$',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    plt.plot(listx,listy_im_ana_par_2,'.-',ms = ms,color = 'purple')
    plt.tight_layout()
    #if plot_vs_c == 1:
    #    plt.yscale('log')
    os.chdir(path_save)
    plt.savefig( 'EELS_div_tot_ana_par_' + label1 + '.png', format='png')   
    
    
    graph(title + '\n' + '$\mathbf{p}$ $\parallel$ $\hat{z}$',labelx,'$\Gamma_{SP}/\Gamma$',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    plt.plot(listx,listy_im_ana_perp_2,'.-',ms = ms,color = 'purple')
    plt.tight_layout()
    #if plot_vs_c == 1:
    #    plt.yscale('log')
    os.chdir(path_save)
    plt.savefig( 'EELS_div_tot_ana_perp_' + label1 + '.png', format='png')   
    
    
    graph(title + '\n' + '$\mathbf{p}$ $\parallel$ $\hat{x}$' ,labelx,'$\Gamma_{SP}$',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    plt.plot(listx,listy_im_ana_par,'.-',ms = ms,color = 'purple')
    plt.tight_layout()
    #if plot_vs_c == 1:
    #    plt.yscale('log')
    os.chdir(path_save)
    plt.savefig( 'EELS_ana_par_' + label1 + '.png', format='png')   
    
    
    graph(title + '\n' + '$\mathbf{p}$ $\parallel$ $\hat{z}$',labelx,'$\Gamma_{SP}$',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    plt.plot(listx,listy_im_ana_perp,'.-',ms = ms,color = 'purple')
    plt.tight_layout()
    #if plot_vs_c == 1:
    #    plt.yscale('log')
    os.chdir(path_save)
    plt.savefig( 'EELS_ana_perp_' + label1 + '.png', format='png')   


#%%

if plot_color_maps == 1:
    
    labelx = r'$\hbar\omega$ [eV]'
    labely =  r'$z_p$ [nm]'

    listx = np.linspace(0.1,3.15,N)
    listy = np.linspace(0.4,50,N)
    
    X, Y = np.meshgrid(listx, listy)
    
    f_pole_approx_par = np.vectorize(function_imag_ana_par_3D)
    Z_pole_approx_par = f_pole_approx_par(X, Y)

    limits = [np.min(listx) , np.max(listx), np.min(listy) , np.max(listy)]
    graph(title  + '\n' + '$\mathbf{p}$ $\parallel$ $\hat{x}$' ,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    #for x in [x1,x2,x3,x4]:
    #    plt.plot(ejey_aux1,x*np.ones(10),'--',color = 'grey' )
    im = plt.imshow(Z_pole_approx_par, extent = limits, cmap=plt.cm.hot, aspect='auto', interpolation = 'none') 

    cbar = plt.colorbar(im, extend='both', fraction=0.046, pad=0.04)
    cbar.ax.tick_params(labelsize = tamnum)
    if paper == 0:
        cbar.set_label('$\Gamma_{SP}$',fontsize=tamlegend,labelpad = 1)
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'EELS_ana_parp_3D_' +  label1 + '.png', format='png')   



    f_pole_approx_perp = np.vectorize(function_imag_ana_perp_3D)
    Z_pole_approx_perp = f_pole_approx_perp(X, Y)

    limits = [np.min(listx) , np.max(listx), np.min(listy) , np.max(listy)]
    graph(title  + '\n' + '$\mathbf{p}$ $\parallel$ $\hat{z}$' ,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    #for x in [x1,x2,x3,x4]:
    #    plt.plot(ejey_aux1,x*np.ones(10),'--',color = 'grey' )
    im = plt.imshow(Z_pole_approx_par, extent = limits, cmap=plt.cm.hot, aspect='auto', interpolation = 'none') 

    cbar = plt.colorbar(im, extend='both', fraction=0.046, pad=0.04)
    cbar.ax.tick_params(labelsize = tamnum)
    if paper == 0:
        cbar.set_label('$\Gamma_{SP}$',fontsize=tamlegend,labelpad = 1)
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'EELS_ana_perp_3D_' +  label1 + '.png', format='png')   


#%%
