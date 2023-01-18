
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: leila

graficar el campo externo directo con la convencion de z hacia abajo
en z = 0
graficar mapa de color x,y
"""
from scipy.signal import find_peaks
import numpy as np
import sys
import os 
import matplotlib.pyplot as plt
#import seaborn as sns

plot_num = 0
plot_color_map = 0

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/potential_field/infinite_dipoles','')
path_save = path_basic + '/' + 'many_dipoles_phi_ind3D'
if not os.path.exists(path_save):
    print('Creating folder to save graphs')
    os.mkdir(path_save)

err = 'potential_many_dipoles_ind.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from potential_many_dipoles_ind import phi_many_dipoles_ana_n
except ModuleNotFoundError:
    print(err)

try:
    sys.path.insert(1, path_constants)
    from graphene_sigma import sigma_DL
except ModuleNotFoundError:
    print('graphene_sigma.py no se encuentra en ' + path_basic)

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


epsi1, epsi2 = 1,1
hbmu, hbgama = 0.3,0.0001

zp = 0.05
b = - 0.01

z0 = zp 

#Ndipoles = 50
#list_xD = np.linspace(-0.01,0.01,Ndipoles)



list_n = [1,2,3,4,5]
list_a = [0.005,0.001,0.01,0.05,0.01]
list_int_v = [10, 5, 2.5, 1.67,1.25]



energy0_pol = 43
omega0 = energy0_pol*1e-3/hb 
#R = 10 # 10 nm en unidades de micrometros
kappa_factor_omega0 = 0.1
kappa_r_factor= 0.5



title2 = r'b = %i nm, $z_p$=%i nm, z=%i nm' %(b*1e3,zp*1e3,z0*1e3) 
title3 = r'$\hbar\omega_0$ = %i meV' %(energy0_pol)



labelx = r'$\theta_{plot}$ [degrees]'  
labely = r'$\theta_{formula}$ [degrees]'  
title = title2 + ', ' +  title3 

listy = np.linspace(50,5000,100)
listx = np.linspace(-1000,1000,100)
listy = np.linspace(0,2000,100)

#listy = np.linspace(-cota_x,cota_x,100)
#listx = listy

#%%

tamfig = (4.5,3.5)
tamlegend = 12
tamletra = 12
tamtitle = 11
tamnum = 10
labelpady = -1.5
labelpadx = 2
pad = 0
mk = 2
lw = 1.5
ms = 4
hp = 0.3
length_marker = 1

#%%
    
def graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad):
    plt.figure(figsize=tamfig)
    plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
    plt.ylabel(labely,fontsize=tamletra,labelpad =labelpadx)
    plt.tick_params(labelsize = tamnum, pad = pad)
    plt.title(title,fontsize=int(tamtitle*0.9))

    return   

#%%
    
def omega_n_THz(int_v,a_nm,Nmax):  ## omega puede valer esto o menos para que se generen plasmones
    
    a = a_nm*1e-3
    
    aux = alfac*int_v*hbmu/(hb)
    rta = aux + 0.5*np.sqrt((2*aux)**2 + 16*alfac*c*hbmu*Nmax*np.pi/(a*hb))
    
    return rta



def theta_from_formula(int_v,a_nm,Nmax):

    omega_n = omega_n_THz(int_v,a_nm,Nmax)
    omegac = omega_n/c
    
 #   omegac = omegac - 0.5 ## ver si funciona con una freq menor 

    E = omegac*aux    
    a = a_nm*1e-3
#    print('omega/c:', omegac)

    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    alfa_p = 1j*(epsi1 + epsi2)/(cond)
    kp = alfa_p*omegac    
    
    lambdda_p = 2*pi/kp
    

    theta0 = np.arccos((omegac*int_v + 2*np.pi*Nmax/a)/np.real(kp))
    

    return theta0*180/np.pi


def theta_from_plot(int_v,a_nm,Nmax):


    omega_n = omega_n_THz(int_v,a_nm,Nmax)
    omegac0 = omega_n/c
    E = omegac0*aux  

    a = a_nm*1e-3

    print('E en meV:', E*1e3)   

    def function_real_ana(x_nano,y_nano):
        y = y_nano*1e-3 
        x = x_nano*1e-3                 
        rta = phi_many_dipoles_ana_n(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,zp,x,y,z0,a,b,Nmax,omega0,kappa_factor_omega0,kappa_r_factor)
        # print(rta)
        return rta.real

    
    
    list_2D_re_x0 = []
    list_2D_re_y0 = []
    x0_fijo = 0
    y0_fijo = 1000
    for x_nano in listx:
        rta = function_real_ana(x_nano,y0_fijo)
        list_2D_re_x0.append(rta)
    
    for y_nano in listy:
        rta = function_real_ana(x0_fijo,y_nano)
        list_2D_re_y0.append(rta)
    

    
    peaks_x01, _ = find_peaks(list_2D_re_x0, height=0)
    peaks_x02, _ = find_peaks(-np.array(list_2D_re_x0), height=0)
    
    peaks_y01, _ = find_peaks(list_2D_re_y0, height=0)
    peaks_y02, _ = find_peaks(-np.array(list_2D_re_y0), height=0)
    
    list_2D_re_x0 = np.array(list_2D_re_x0)
    list_2D_re_y0 = np.array(list_2D_re_y0)
    
    
    delta_x_maximos = []    ## distancia entre maximos
    delta_x_minimos = []    ## distancia entre minimos
    for j in range(len(peaks_x01)-1):
        delta_x_maximos.append(listx[peaks_x01[j+1]] - listx[peaks_x01[j]])
    
    for j in range(len(peaks_x02)-1):
        delta_x_minimos.append(listx[peaks_x02[j+1]] - listx[peaks_x02[j]])
    
    
    delta_y_maximos = []     ## distancia entre maximos
    delta_y_minimos = []       ## distancia entre minimos
    
    for j in range(len(peaks_y01)-1):
        delta_y_maximos.append(listx[peaks_y01[j+1]] - listx[peaks_y01[j]])
    
    for j in range(len(peaks_y02)-1):
        delta_y_minimos.append(listx[peaks_y02[j+1]] - listx[peaks_y02[j]])
    
    
    if len(delta_x_maximos) + len(delta_x_minimos) >= 2 :
        delta_x_final = np.mean([np.mean(delta_x_maximos), np.mean(delta_x_minimos)])
    else:
        delta_x_final = np.abs((listx[peaks_x01[0]] - listx[peaks_x02[0]])*2)
    
    
    delta_y_final = np.mean([np.mean(delta_y_maximos), np.mean(delta_y_minimos)])
    
    theta_plot = np.arctan2(delta_y_final,delta_x_final)   # # # en radianes


    if np.isnan(theta_plot):
        raise TypeError('wrong value for theta')


    return theta_plot*180/np.pi

#%%#    
    
list_theta_plot_tot = []
list_theta_formula_tot = []

for j in range(len(list_n)):
    
    Nmax = list_n[j]
    a = list_a[j]
    a_nm = a*1e3
    int_v = list_int_v[j]
    
    
    theta_plot_value = theta_from_plot(int_v,a_nm,Nmax)
    theta_formula_value = theta_from_formula(int_v,a_nm,Nmax)
    
    list_theta_formula_tot.append(theta_formula_value)
    list_theta_plot_tot.append(theta_plot_value)




#%%#   

colors = ['darkred','steelblue','coral','yellowgreen']
 
graph(title ,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(list_theta_plot_tot,list_theta_formula_tot,'.-', color = colors[1], ms = ms)
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
os.chdir(path_save)
plt.grid(1)
plt.savefig( 'theta_formula_vs_theta_plot'  + '.png', format='png')       


graph(title ,'n',labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(list_n,list_theta_formula_tot,'.-', color = colors[1], ms = ms)
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
os.chdir(path_save)
plt.grid(1)
plt.savefig( 'theta_formula_vs_n'  + '.png', format='png')     



graph(title ,'a [nm]',labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(list_a,list_theta_formula_tot,'.-', color = colors[1], ms = ms)
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
os.chdir(path_save)
plt.grid(1)
plt.savefig( 'theta_formula_vs_a'  + '.png', format='png')   







graph(title ,'v/c',labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(list_int_v,list_theta_formula_tot,'.-', color = colors[1], ms = ms)
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
os.chdir(path_save)
plt.grid(1)
plt.savefig( 'theta_formula_vs_v'  + '.png', format='png')   


#%%# 