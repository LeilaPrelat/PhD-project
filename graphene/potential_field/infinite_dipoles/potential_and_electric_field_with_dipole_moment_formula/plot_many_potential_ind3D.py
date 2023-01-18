
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
#import seaborn as sns

plot_num = 0
plot_color_map = 1
plot_check_position_disp_relation = 0

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
    from potential_many_dipoles_ind import phi_many_dipoles_ana_n, phi_many_dipoles_pole_aprox_n
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

Nmax = 5


int_v0 = 20
int_v0 = 10


a = 5
a_nm = a*1e3


def omega_n_THz(int_v,a_nm,Nmax):  ## omega puede valer esto o menos para que se generen plasmones
    
    a = a_nm*1e-3
    
    aux = alfac*int_v*hbmu/(hb)
    rta = aux + 0.5*np.sqrt((2*aux)**2 + 16*alfac*c*hbmu*Nmax*np.pi/(a*hb))
    
    return rta*1.5 ## menos freq 



def theta(int_v,a_nm,Nmax):

    omega_n = omega_n_THz(int_v,a_nm,Nmax)
    omegac = omega_n/c

    
 #   omegac = omegac - 0.5 ## ver si funciona con una freq menor 

    E = omegac*aux    

    print('omega/c:', omegac)

    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    alfa_p = 1j*(epsi1 + epsi2)/(cond)
    kp = alfa_p*omegac    
    
    lambdda_p = 2*pi/kp
    
#    
#    den1 = omegac*int_v/(2*np.pi) - 1/lambdda_p
#    den2 = omegac*int_v/(2*np.pi) + 1/lambdda_p
#    
#    a_min = Nmax/den1
#    
#    a_max = Nmax/den2
#    
#    a = np.mean([a_min,a_max])
#    print('a:', a)
    
    theta0 = np.arccos(np.real(lambdda_p)*(omegac*int_v/(2*np.pi) - Nmax/a))
    
    return theta0 

#omega_n = omega_n_THz(int_v,a_nm,Nmax)
#lambda_SP = 2*pi*v/omega_n 
#omegac0 = omega_n/c
#E0 = omegac0*(c*hb)
#
#############################################
#E0 = 43
#def k_SP(E_mev,hbmu,hbgama):
#    num = 2/alfac 
#    omegac = E_mev*1e-3/aux
#    cte = E_mev*1e-3 + 1j*hbgama/(4*hbmu)
#    
#    return omegac*num*cte
#
#lambda_SP = 2*np.pi/k_SP(E0,hbmu,hbgama)
#
#a_min = np.real(lambda_SP)*Nmax/(int_v - 1)
#a_max = np.real(lambda_SP)*Nmax/(int_v + 1)
#
#a = np.mean([a_min,a_max])
#a = a_min

############################################
theta0 = theta(int_v0,a_nm,Nmax)
print('theta:', theta0*180/np.pi, 'degrees')
if np.isnan(theta0):
    raise TypeError('wrong value for theta')

energy0_pol = 43
omega0 = energy0_pol*1e-3/hb 
#R = 10 # 10 nm en unidades de micrometros
kappa_factor_omega0 = 0.1
kappa_r_factor= 0.5

omega_n = omega_n_THz(int_v0,a_nm,Nmax)
omegac0 = omega_n/c
E = omegac0*aux  

print('E en meV:', E*1e3)

title2 = r'v = c/%i, b = %i nm, $z_p$=%i nm, z=%i nm' %(int_v0,b*1e3,zp*1e3,z0*1e3) 
title3 = r'$\hbar\omega_0$ = %i meV, a = %i nm, $\theta$ = %.2f$^º$, n = %i' %(energy0_pol,a*1e3,theta0*180/np.pi,Nmax)

labelp = r'__a%inm_n%i' %(a*1e3,Nmax)

labelx = 'x [nm]'  
labely = 'y [nm]'

title = title2 + '\n' +  title3 + '\n'  + 'E = %.2f meV' %(E*1e3)
label1 = 'vs_y_theta0%.2f' %(theta0) + labelp
listy = np.linspace(50,5000,100)
Nlist = 100
#Nlist = 500
listx = np.linspace(-1000,1000,Nlist)
listy = np.linspace(0,2000,Nlist)

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

def disp_relation_omega(kx,a,Nmax,epsi1,epsi2,hbmu):

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

    cte = (epsi1 + epsi2)
    
    kn = kx + 2*np.pi*Nmax/a
    
    
    return np.sqrt(kn*4*hbmu*alfac*(c/hb)/cte)


def omega_n_corte(int_v,a_nm,Nmax):  ## omega puede valer esto o menos para que se generen plasmones
    
    a = a_nm*1e-3
    
    aux = alfac*int_v*hbmu/(hb)
    rta = aux + 0.5*np.sqrt((2*aux)**2 + 16*alfac*c*hbmu*Nmax*np.pi/(a*hb))
#    
#    rta2 = aux - 0.5*np.sqrt((2*aux)**2 + 16*alfac*c*hbmu*Nmax*np.pi/(a*hb))
    
#    print(rta*hb*1e3)
    
    return rta   

if plot_check_position_disp_relation == 1:
    list_colours = ['black', 'gold', 'orange', 'darksalmon', 'brown', 'darkred']
    title1_aux = '$\epsilon_1 = %.1f$, $\epsilon_2 = %.1f$, $\gamma_{in}$ = %.2f meV' %(epsi1,epsi2,hbgama*1e3)
    title2_aux = '$\mu_c$ = %.2f eV, a = %i nm, n = %i' %(hbmu,a*1e3,Nmax)
    
    
    Ntot = 200
    list_kx = np.linspace(0.001, 4 + Nmax/(a), Ntot) ##creo que va a estar en 1/micrones
    
    list_y_re = []
    list_y_im = []
    for kx in list_kx:
        
        valuey = disp_relation_omega(kx,a,Nmax,epsi1,epsi2,hbmu)
        valuey = hb*valuey*1e3
        list_y_re.append(valuey.real)
        list_y_im.append(valuey.imag)
        
    
    plt.figure(figsize=tamfig)    
    plt.title(title1_aux + '\n' + title2_aux ,fontsize=tamtitle)
    plt.plot(list_kx,list_y_re,'-',color = 'darkgreen',ms = lw, label = 'DR')
    list_hb_omega_meV = np.linspace(np.min(list_y_re),np.max(list_y_re),Ntot)
    j = 0
    list_cortes_omega = []
    list_cortes_k = []
    list_int_v = [1/int_v0]
    for int_v in list_int_v:
        v = c*int_v
        
        list_y2_re = []
        
        for hb_omega_meV in list_hb_omega_meV:
            omega = hb_omega_meV*1e-3/hb
            
            valuey_e = omega/v
            list_y2_re.append(valuey_e.real)
            
        
        hb_omega_corte = omega_n_corte(1/int_v,a_nm,Nmax)*hb*1e3
        
        k_corte = omega_n_corte(1/int_v,a_nm,Nmax)/v
        
        list_cortes_omega.append(hb_omega_corte)
        list_cortes_k.append(k_corte)
            
        ##  
        if int_v != 1:
            plt.plot(list_y2_re,list_hb_omega_meV,'-', color = list_colours[j],ms = lw, label = 'v= %.1fc' %(int_v))
            plt.plot(list_cortes_k,list_cortes_omega,'.', color = 'red',ms = '3')
        else:
            plt.plot(list_y2_re,list_hb_omega_meV,'-', color = list_colours[j], ms = lw, label = 'light')
            
        j = j + 1
         
    
    plt.xlabel('$k_\parallel$ [1/$\mu$m]',fontsize=tamletra, labelpad = 0)
    plt.ylabel('$\hbar\omega$ [meV]',fontsize=tamletra, labelpad = 0)
    plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
    plt.tick_params(labelsize = tamnum,pad = pad)
    plt.grid(1)
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig('disp_relation_for_potential' + label1 + '.png', format='png')

#%%

def function_real_ana(x_nano,y_nano):
    y = y_nano*1e-3 
    x = x_nano*1e-3                 
    rta = phi_many_dipoles_ana_n(omegac0,epsi1,epsi2,hbmu,hbgama,int_v0,zp,x,y,z0,a,b,Nmax,omega0,kappa_factor_omega0,kappa_r_factor)
    # print(rta)
    return rta.real

def function_imag_ana(x_nano,y_nano):
    y = y_nano*1e-3 
    x = x_nano*1e-3 
                                 
    rta = phi_many_dipoles_ana_n(omegac0,epsi1,epsi2,hbmu,hbgama,int_v0,zp,x,y,z0,a,b,Nmax,omega0,kappa_factor_omega0,kappa_r_factor)
    return rta.imag

#%%
    
def function_real_pole_aprox(x_nano,y_nano):
    y = y_nano*1e-3 
    x = x_nano*1e-3 
                                          
    rta = phi_many_dipoles_pole_aprox_n(omegac0,epsi1,epsi2,hbmu,hbgama,int_v0,zp,x,y,z0,a,b,Nmax,omega0,kappa_factor_omega0,kappa_r_factor)
    # print(rta)
    return rta.real

def function_imag_pole_aprox(x_nano,y_nano):
    
    y = y_nano*1e-3 
    x = x_nano*1e-3 
    
    rta = phi_many_dipoles_pole_aprox_n(omegac0,epsi1,epsi2,hbmu,hbgama,int_v0,zp,x,y,z0,a,b,Nmax,omega0,kappa_factor_omega0,kappa_r_factor)
    return rta.imag

#%%
    


def graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad):
    plt.figure(figsize=tamfig)
    plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
    plt.ylabel(labely,fontsize=tamletra,labelpad =labelpadx)
    plt.tick_params(labelsize = tamnum, pad = pad)
    plt.title(title,fontsize=int(tamtitle*0.9))

    return   

#%% analitico
#        
if plot_color_map == 1:
    
    N = 100
    X, Y = np.meshgrid(listx, listy)
    limits = [min(listx) , max(listx), min(listy) , max(listy)]
    
    aux_list = np.array(np.ones(N))
    
    f_real_ana = np.vectorize(function_real_ana)
    Z_real_ana = f_real_ana(X, Y)
    
    f_imag_ana = np.vectorize(function_imag_ana)
    Z_imag_ana = f_imag_ana(X, Y)
    
       
    print('Graficar el External field' + label1)
    
    # if plot_xz_3D == 1 and px == 1: 
    graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    vmin, vmax = np.min(Z_real_ana), np.max(Z_real_ana)
    im = plt.imshow(Z_real_ana, extent = limits, cmap=plt.cm.hot, interpolation='none',origin = 'lower') 
    cbar = plt.colorbar(im)
    #pcm = plt.pcolormesh(X, Y, Z_real_ana,
    #                      norm=colors.SymLogNorm(linthresh=0.03, linscale=0.03,
    #                                           vmin=int(vmin), vmax=int(vmax)),cmap=plt.cm.hot)
    #maxlog=int(np.ceil( np.log10( np.abs(vmax) )))
    #minlog=int(np.ceil( np.log10( np.abs(vmin) )))
    #if vmin < 0 :
    #      tick_locations = ( [-(10.0**x) for x in np.linspace(minlog,-1,8)] 
    #                        + [0] 
    #                        + [(10.0**x) for x in np.linspace(-1,maxlog,8)] )
    #else:
    #      tick_locations = ( [(10.0**x) for x in np.linspace(minlog,maxlog,maxlog + np.abs(minlog) + 1) ])    
    #    im = plt.imshow(Z1_real, extent = limits, cmap=plt.cm.hot, interpolation='bilinear')
    #cbar = plt.colorbar(pcm)
    cbar.ax.tick_params(labelsize = tamnum)
    #cbar.set_ticks(tick_locations)
    cbar.set_label(r'Re($\phi_n$)')
    os.chdir(path_save)
    plt.tight_layout()
    plt.savefig('Re_potential_ana' + label1 + '.png', format='png')
    
    
    
    
    graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    vmin, vmax = np.min(Z_imag_ana), np.max(Z_imag_ana)
    im = plt.imshow(Z_imag_ana, extent = limits, cmap=plt.cm.hot, interpolation='none',origin = 'lower') 
    cbar = plt.colorbar(im)
    #pcm = plt.pcolormesh(X, Y, Z_imag_ana,
    #                      norm=colors.SymLogNorm(linthresh=0.03, linscale=0.03,
    #                                           vmin=int(vmin), vmax=int(vmax)),cmap=plt.cm.hot)
    #maxlog=int(np.ceil( np.log10( np.abs(vmax) )))
    #minlog=int(np.ceil( np.log10( np.abs(vmin) )))
    #if vmin < 0 :
    #      tick_locations = ( [-(10.0**x) for x in np.linspace(minlog,-1,8)] 
    #                        + [0] 
    #                        + [(10.0**x) for x in np.linspace(-1,maxlog,8)] )
    #else:
    #      tick_locations = ( [(10.0**x) for x in np.linspace(minlog,maxlog,maxlog + np.abs(minlog) + 1) ])    
    #    im = plt.imshow(Z1_real, extent = limits, cmap=plt.cm.hot, interpolation='bilinear')
    #cbar = plt.colorbar(pcm)
    cbar.ax.tick_params(labelsize = tamnum)
    #cbar.set_ticks(tick_locations)
    cbar.set_label(r'Im($\phi_n$)')
    os.chdir(path_save)
    plt.tight_layout()
    plt.savefig('Im_potential_ana' + label1 + '.png', format='png')


#%%

# chequear si el angulo es el posta
    

list_2D_re_x0 = []
list_2D_re_y0 = []
x0_fijo = 0
y0_fijo = 1000
for x_nano in listx:
    rta = function_real_ana(x_nano,y0_fijo)
    list_2D_re_x0.append(rta)

listy_2 = np.linspace(-4000,4000,100)
for y_nano in listy:
    rta = function_real_ana(x0_fijo,y_nano)
    list_2D_re_y0.append(rta)

#%%
from scipy.signal import find_peaks
colors = ['darkred','steelblue','coral','yellowgreen']

peaks_x01, _ = find_peaks(list_2D_re_x0, height=0)
peaks_x02, _ = find_peaks(-np.array(list_2D_re_x0), height=0)

peaks_y01, _ = find_peaks(list_2D_re_y0, height=0)
peaks_y02, _ = find_peaks(-np.array(list_2D_re_y0), height=0)

list_2D_re_x0 = np.array(list_2D_re_x0)
list_2D_re_y0 = np.array(list_2D_re_y0)

#%%

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


#%%
if len(delta_x_maximos) + len(delta_x_minimos) >= 2 :
    delta_x_final = np.mean([np.mean(delta_x_maximos), np.mean(delta_x_minimos)])
else:
    delta_x_final = np.abs((listx[peaks_x01[0]] - listx[peaks_x02[0]])*2)


delta_y_final = np.mean([np.mean(delta_y_maximos), np.mean(delta_y_minimos)])

theta_plot = np.arctan2(delta_y_final,delta_x_final)   # # # en radianes

print('(theta grafico + theta formula)/pi:', (theta_plot + theta0 )/np.pi )


#%%
    
graph(title + r', y = %i nm, $\theta_{plot}$ = %.2f' %(y0_fijo,theta_plot) ,'x [nm]','Re{$\phi_n$}',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(listx,list_2D_re_x0,'.-', color = colors[1], ms = ms)
plt.plot(listx[peaks_x01], list_2D_re_x0[peaks_x01],'.', color = 'red')
plt.plot(listx[peaks_x02], list_2D_re_x0[peaks_x02],'.', color = 'red')
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
os.chdir(path_save)
plt.savefig( 'phi_yfijo' + label1 + '.png', format='png')   



graph(title  + r', x = %i nm, $\theta_{plot}$ = %.2f' %(x0_fijo,theta_plot) ,'y [nm]','Re{$\phi_n$}',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(listy,list_2D_re_y0,'.-', color = colors[1], ms = ms)
plt.plot(listy[peaks_y01], list_2D_re_y0[peaks_y01], '.', color = 'red')
plt.plot(listy[peaks_y02], list_2D_re_y0[peaks_y02], '.', color = 'red')
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
os.chdir(path_save)
plt.savefig( 'phi_xfijo' + label1 + '.png', format='png')   


#%%

#
if plot_num == 1:   
    
    f_real_pole_aprox = np.vectorize(function_real_pole_aprox)
    Z_real_pole_aprox = f_real_ana(X, Y)
    
    f_imag_pole_aprox = np.vectorize(function_imag_pole_aprox)
    Z_imag_pole_aprox = f_imag_pole_aprox(X, Y)
       
    print('Graficar el External field' + label1)
    
    # if plot_xz_3D == 1 and px == 1: 
    graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    vmin, vmax = np.min(Z_real_pole_aprox), np.max(Z_real_pole_aprox)
#    pcm = plt.pcolormesh(X, Y, Z_real_pole_aprox,
#                          norm=colors.SymLogNorm(linthresh=0.03, linscale=0.03,
#                                               vmin=int(vmin), vmax=int(vmax)),cmap=plt.cm.hot)
#    maxlog=int(np.ceil( np.log10( np.abs(vmax) )))
#    minlog=int(np.ceil( np.log10( np.abs(vmin) )))
#    if vmin < 0 :
#          tick_locations = ( [-(10.0**x) for x in np.linspace(minlog,-1,8)] 
#                            + [0] 
#                            + [(10.0**x) for x in np.linspace(-1,maxlog,8)] )
#    else:
#          tick_locations = ( [(10.0**x) for x in np.linspace(minlog,maxlog,maxlog + np.abs(minlog) + 1) ])    
    #    im = plt.imshow(Z1_real, extent = limits, cmap=plt.cm.hot, interpolation='bilinear')
#    cbar = plt.colorbar(pcm)

    im = plt.imshow(Z_real_pole_aprox, extent = limits, cmap=plt.cm.hot, interpolation='none')
    cbar = plt.colorbar(im)
    cbar.ax.tick_params(labelsize = tamnum)
#    cbar.set_ticks(tick_locations)
    cbar.set_label(r'Re($\phi_n$)')
    os.chdir(path_save)
    plt.tight_layout()
    plt.savefig('Re_potential_pole_aprox_log' + label1 + '.png', format='png')
    
    
    
    
    graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    vmin, vmax = np.min(Z_imag_pole_aprox), np.max(Z_imag_pole_aprox)
#    pcm = plt.pcolormesh(X, Y, Z_imag_pole_aprox,
#                          norm=colors.SymLogNorm(linthresh=0.03, linscale=0.03,
#                                               vmin=int(vmin), vmax=int(vmax)),cmap=plt.cm.hot)
#    maxlog=int(np.ceil( np.log10( np.abs(vmax) )))
#    minlog=int(np.ceil( np.log10( np.abs(vmin) )))
#    if vmin < 0 :
#          tick_locations = ( [-(10.0**x) for x in np.linspace(minlog,-1,8)] 
#                            + [0] 
#                            + [(10.0**x) for x in np.linspace(-1,maxlog,8)] )
#    else:
#          tick_locations = ( [(10.0**x) for x in np.linspace(minlog,maxlog,maxlog + np.abs(minlog) + 1) ])    
    #    im = plt.imshow(Z1_real, extent = limits, cmap=plt.cm.hot, interpolation='bilinear')
#    cbar = plt.colorbar(pcm)
    im = plt.imshow(Z_imag_pole_aprox, extent = limits, cmap=plt.cm.hot, interpolation='none')
    cbar = plt.colorbar(im)
    cbar.ax.tick_params(labelsize = tamnum)
#    cbar.set_ticks(tick_locations)
    cbar.set_label(r'Im($\phi_n$)')
    os.chdir(path_save)
    plt.tight_layout()
    plt.savefig('Im_potential_pole_aprox_log' + label1 + '.png', format='png')
#    