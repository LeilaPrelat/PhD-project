
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

#plot_color_map = 1
plot_check_position_disp_relation = 0
paper = 1

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/potential_field/infinite_dipoles','')
path_save = path_basic + '/' + 'find_lambda_p_many_potential'
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


energy0_pol = 43
omega0 = energy0_pol*1e-3/hb 
#R = 10 # 10 nm en unidades de micrometros
kappa_factor_omega0 = 0.1
kappa_r_factor= 0.5

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
    if paper == 0:
        plt.title(title,fontsize=int(tamtitle*0.9))

    return   

#%%

# chequear si el angulo es el posta
labelx = 'x [nm]'  
labely = 'y [nm]'
Nlist = 75
#Nlist = 500



#int_v0 = 20
int_v0 = 10
#
#
a0 = 1
a0_nm = a0*1e3
Nmax0 = 5


if Nmax0 == 0 and a0 > 0.05:
    listx = np.linspace(-3000,3000,Nlist)
    listy = np.linspace(0,6000,Nlist) 
elif Nmax0 == 0 and a0 <= 0.05:
    listx = np.linspace(-1000,1000,Nlist)
    listy = np.linspace(0,2000,Nlist) 

else:
    listx = np.linspace(-1000,1000,Nlist)
    listy = np.linspace(0,2000,Nlist) 
    


#Nmax = 3
#
#
#int_v0 = 20
#int_v0 = 10
#
#
#a = 5
#a_nm = a*1e3

#%%


def omega_n_THz(int_v,a_nm,Nmax):  ## omega puede valer esto o menos para que se generen plasmones
    
    a = a_nm*1e-3
    
    aux = alfac*int_v*hbmu/(hb)
    rta = aux + 0.5*np.sqrt((2*aux)**2 + 16*alfac*c*hbmu*Nmax*np.pi/(a*hb))
    
    return rta*1.5 ## menos freq 


def theta_formula_function(int_v,a_nm,Nmax):

    omega_n = omega_n_THz(int_v,a_nm,Nmax)
    omegac = omega_n/c

    
 #   omegac = omegac - 0.5 ## ver si funciona con una freq menor 

    E = omegac*aux    

#    print('omega/c:', omegac)

    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    alfa_p = 1j*(epsi1 + epsi2)/(cond)
    kp = alfa_p*omegac    
#    
#    lambdda_p = 2*pi/kp
    
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
    a = a_nm*1e-3
    theta0 = np.arccos((omegac*int_v + 2*np.pi*Nmax/a)/np.real(kp))
    
    return theta0

def theta_from_plot(int_v,a_nm,Nmax,plot_true_or_false):
    a = a_nm*1e-3
    
    theta0 = theta_formula_function(int_v,a_nm,Nmax)

#    print('theta:', theta0*180/np.pi, 'degrees')
#    if np.isnan(theta0):
#        raise TypeError('wrong value for theta')
    
    
    
    omega_n = omega_n_THz(int_v,a_nm,Nmax)
    omegac0 = omega_n/c
    E = omegac0*aux  
    
#    print('E en meV:', E*1e3)
    
    title2 = r'v = c/%i, b = %i nm, $z_p$=%i nm, z=%i nm' %(int_v,b*1e3,zp*1e3,z0*1e3) 
    title3 = r'$\hbar\omega_0$ = %i meV, a = %i nm, $\theta_{eq}$ = %.2f$^º$, n = %i' %(energy0_pol,a*1e3,theta0*180/np.pi,Nmax)
    
    title = title2 + '\n' +  title3 + '\n'  + 'E = %.2f meV' %(E*1e3)
    labelp = r'_a%inm_n%i' %(a*1e3,Nmax)
    label1 = '_theta0%.2f' %(theta0) + labelp

    
    
    from scipy.signal import find_peaks


#    print('E en meV:', E*1e3)   

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
    
    colors = ['darkred','steelblue','coral','yellowgreen']
    
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


#    if np.isnan(theta_plot):
#        raise TypeError('wrong value for theta')


#    print('(theta grafico + theta formula)/pi:', (theta_plot + theta0 )/np.pi )

    if plot_true_or_false == 1:
        graph(title + r', y = %i nm, $\theta_{plot}$ = %.2f$^º$' %(y0_fijo,theta_plot*180/np.pi) ,'x [nm]','Re{$\phi_n$}',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
        plt.plot(listx,list_2D_re_x0,'.-', color = colors[1], ms = ms)
        plt.plot(listx[peaks_x01], list_2D_re_x0[peaks_x01],'.', color = 'red')
        plt.plot(listx[peaks_x02], list_2D_re_x0[peaks_x02],'.', color = 'red')
#        plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
        plt.tight_layout()
        os.chdir(path_save)
        plt.savefig( 'phi_yfijo' + label1 + '.png', format='png')   
        
        
        
        graph(title  + r', x = %i nm, $\theta_{plot}$ = %.2f$^º$' %(x0_fijo,theta_plot*180/np.pi) ,'y [nm]','Re{$\phi_n$}',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
        plt.plot(listy,list_2D_re_y0,'.-', color = colors[1], ms = ms)
        plt.plot(listy[peaks_y01], list_2D_re_y0[peaks_y01], '.', color = 'red')
        plt.plot(listy[peaks_y02], list_2D_re_y0[peaks_y02], '.', color = 'red')
#        plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
        plt.tight_layout()
        os.chdir(path_save)
        plt.savefig( 'phi_xfijo' + label1 + '.png', format='png')   



    return theta_plot

def recta_theta_from_plot(int_v,a_nm,Nmax,x_nano,y0,plot_true_or_false):
    theta_plot = theta_from_plot(int_v,a_nm,Nmax,plot_true_or_false)
    
    pendiente = np.tan(theta_plot)
    
    y_nano = x_nano*pendiente + y0
    
    return y_nano  


def recta_for_find_lambda_p(int_v,a_nm,Nmax,x_nano,y0_nano,plot_true_or_false):
    theta_plot = theta_from_plot(int_v,a_nm,Nmax,plot_true_or_false)
    
    pendiente = np.tan(theta_plot)
    
    y_nano = -x_nano/pendiente + y0_nano
    
    return y_nano
    



def find_lambda_p_from_plot(int_v,a_nm,Nmax,y0_nano,plot_true_or_false):

    theta0 = theta_formula_function(int_v,a_nm,Nmax)
    
    omega_n = omega_n_THz(int_v,a_nm,Nmax)
    omegac0 = omega_n/c
    E = omegac0*aux  
#    
#    print('E en meV:', E*1e3)
    a = a_nm*1e-3
    title2 = r'v = c/%i, b = %i nm, $z_p$=%i nm, z=%i nm, $\hbar\omega_0$ = %i meV' %(int_v,b*1e3,zp*1e3,z0*1e3,energy0_pol) 
    title3 = r'a = %i nm, $\theta_{eq}$ = %.2f$^º$, n = %i' %(a*1e3,theta0*180/np.pi,Nmax)
    
    title = title2 + '\n' +  title3 + ', '  + 'E = %.2f meV' %(E*1e3)
    labelp = r'_a%inm_n%i' %(a*1e3,Nmax)
    label1 = '_theta0%.2f' %(theta0) + labelp

    
    
    from scipy.signal import find_peaks


#    print('E en meV:', E*1e3)   

    def function_real_ana(x_nano,y_nano):
        y = y_nano*1e-3 
        x = x_nano*1e-3                 
        rta = phi_many_dipoles_ana_n(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,zp,x,y,z0,a,b,Nmax,omega0,kappa_factor_omega0,kappa_r_factor)
        # print(rta)
        return rta.real

    list_2D_re_y0 = []
#    print(listx)
    listx2 = np.linspace(-8000,4500,100)
    for x_nano in listx2:
        
        y_value_recta = recta_for_find_lambda_p(int_v,a_nm,Nmax,x_nano,y0_nano,0)
#        print(y_value_recta)
        rta = function_real_ana(x_nano,y_value_recta)
        list_2D_re_y0.append(rta)
        
        
 
    colors = ['darkred','steelblue','coral','yellowgreen']
    
    theta_plot = theta_from_plot(int_v,a_nm,Nmax,plot_true_or_false)
#    
#    pendiente = np.tan(theta_plot)
#    
    peaks_y01, _ = find_peaks(list_2D_re_y0, height=0)
    peaks_y02, _ = find_peaks(-np.array(list_2D_re_y0), height=0)
    

    list_2D_re_y0 = np.array(list_2D_re_y0)
    

    delta_y_maximos = []     ## distancia entre maximos
    delta_y_minimos = []       ## distancia entre minimos
    
    for j in range(len(peaks_y01)-1):
        delta_y_maximos.append(listx2[peaks_y01[j+1]] - listx2[peaks_y01[j]])
    
    for j in range(len(peaks_y02)-1):
        delta_y_minimos.append(listx2[peaks_y02[j+1]] - listx2[peaks_y02[j]])
    
    
    lambda_p = np.mean([np.mean(delta_y_maximos), np.mean(delta_y_minimos)])  

#    if np.isnan(theta_plot):
#        raise TypeError('wrong value for theta')

    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    alfa_p = 1j*(epsi1 + epsi2)/(cond)
    kp = alfa_p*omegac0
    lambda_p_formula = np.real(2*np.pi/kp)
#    print('(theta grafico + theta formula)/pi:', (theta_plot + theta0 )/np.pi )
    if plot_true_or_false == 1:

        graph(title + '\n'  + r'$\lambda^{plot}_p$ = %.2f nm, $\lambda^{eq}_p$ = %.2f nm, $\theta^{plot}$ = %.2f$^º$' %(lambda_p,lambda_p_formula*1e3,theta_plot*180/np.pi) ,'x [nm]','Re{$\phi_n$}',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
        plt.plot(listx2,list_2D_re_y0,'.-', color = colors[1], ms = ms)
        plt.plot(listx2[peaks_y01], list_2D_re_y0[peaks_y01], '.', color = 'red')
        plt.plot(listx2[peaks_y02], list_2D_re_y0[peaks_y02], '.', color = 'red')
#        plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
        plt.tight_layout()
        os.chdir(path_save)
        plt.savefig( 'phi_recta_para_hallar_lambda_p' + label1 + '.png', format='png')       
        
    
    return lambda_p
    
#%%




theta_formula = theta_formula_function(int_v0,a0_nm,Nmax0)
theta_plot = theta_from_plot(int_v0,a0_nm,Nmax0,0)

print('theta formula:', theta_formula*180/np.pi, 'degrees')
if np.isnan(theta_formula):
    raise TypeError('wrong value for theta')



omega_n = omega_n_THz(int_v0,a0_nm,Nmax0)
omegac0 = omega_n/c
E = omegac0*aux  

print('E en meV:', E*1e3)


if Nmax0 == 0:
    y00 = 5000

    lambda_p_from_plot = find_lambda_p_from_plot(int_v0,a0_nm,Nmax0,y00,1)
else:
    y00 = 3000
    lambda_p_from_plot = find_lambda_p_from_plot(int_v0,a0_nm,Nmax0,y00,1)

print('lambda_p from plot:', lambda_p_from_plot, 'nm')
cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
Rp = 2*epsi1/(epsi1 + epsi2)
alfa_p = 1j*(epsi1 + epsi2)/(cond)
kp = alfa_p*omegac0
lambda_p_formula = np.real(2*np.pi/kp)
print('lambda_p from formula:', lambda_p_formula*1e3, 'nm')


title2 = r'v = c/%i, b = %i nm, $z_p$=%i nm, z=%i nm, $\hbar\omega_0$ = %i meV,' %(int_v0,b*1e3,zp*1e3,z0*1e3,energy0_pol) 
title3 = r' a = %i nm, $\theta_{eq}$ = %.2f$^º$, n = %i, E = %.2f meV' %(a0*1e3,theta_formula*180/np.pi,Nmax0,E*1e3)

title = title2 + '\n' +  title3 + '\n'  + r'$\lambda^{plot}_p$ = %.2f nm, $\lambda^{eq}_p$ = %.2f nm, $\theta^{plot}$ = %.2f$^º$' %(lambda_p_from_plot,lambda_p_formula*1e3,theta_plot*180/np.pi)
labelp = r'_a%inm_n%i' %(a0*1e3,Nmax0)
label1 = '_theta0%.2f' %(theta_formula) + labelp
  
    
def function_real_ana(x_nano,y_nano):
    y = y_nano*1e-3 
    x = x_nano*1e-3                 
    rta = phi_many_dipoles_ana_n(omegac0,epsi1,epsi2,hbmu,hbgama,int_v0,zp,x,y,z0,a0,b,Nmax0,omega0,kappa_factor_omega0,kappa_r_factor)
    # print(rta)
    return rta.real

def function_imag_ana(x_nano,y_nano):
    y = y_nano*1e-3 
    x = x_nano*1e-3 
                                 
    rta = phi_many_dipoles_ana_n(omegac0,epsi1,epsi2,hbmu,hbgama,int_v0,zp,x,y,z0,a0,b,Nmax0,omega0,kappa_factor_omega0,kappa_r_factor)
    return rta.imag



N = 100
X, Y = np.meshgrid(listx, listy)
limits = [min(listx) , max(listx), min(listy) , max(listy)]

aux_list = np.array(np.ones(N))

f_real_ana = np.vectorize(function_real_ana)
Z_real_ana = f_real_ana(X, Y)

f_imag_ana = np.vectorize(function_imag_ana)
Z_imag_ana = f_imag_ana(X, Y)

#%%




recta_para_lambda_p1 = recta_for_find_lambda_p(int_v0,a0_nm,Nmax0,listy,4000,0)   ### rectas perpendiculares
recta_para_lambda_p2 = recta_for_find_lambda_p(int_v0,a0_nm,Nmax0,listy,y00,0)
recta_para_lambda_p3 = recta_for_find_lambda_p(int_v0,a0_nm,Nmax0,listy,2000,0)


recta_from_plot_p1 = recta_theta_from_plot(int_v0,a0_nm,Nmax0,listy,500,0)     ### rectas paralelas : cortes para hallar theta desde el grafico
recta_from_plot_p2 = recta_theta_from_plot(int_v0,a0_nm,Nmax0,listy,250,0)
recta_from_plot_p3 = recta_theta_from_plot(int_v0,a0_nm,Nmax0,listy,100,0)

   
print('Graficar el External field' + label1)

# if plot_xz_3D == 1 and px == 1: 
graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
vmin, vmax = np.min(Z_real_ana), np.max(Z_real_ana)
im = plt.imshow(Z_real_ana, extent = limits, cmap=plt.cm.hot, interpolation='bilinear')
cbar = plt.colorbar(im)
#plt.plot(listx,recta_para_lambda_p1,'-',color = 'green')
if paper == 0:
    plt.plot(listx,recta_para_lambda_p2,'-',color = 'green')

#plt.plot(listx,recta_para_lambda_p3,'-',color = 'green')

#plt.plot(listx,recta_from_plot_p1,'-',color = 'blue')
    plt.plot(listx,recta_from_plot_p2,'-',color = 'blue')
#plt.plot(listx,recta_from_plot_p3,'-',color = 'blue')

cbar.ax.tick_params(labelsize = tamnum)
#cbar.set_ticks(tick_locations)
plt.ylim([listy[0],listy[-1]])
plt.xlim([listx[0],listx[-1]])
cbar.set_label(r'Re{$\phi_%i$}' %(Nmax0))
os.chdir(path_save)
plt.tight_layout()
plt.savefig('Re_potential_ana' + label1 + '.png', format='png')




graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
vmin, vmax = np.min(Z_imag_ana), np.max(Z_imag_ana)
im = plt.imshow(Z_imag_ana, extent = limits, cmap=plt.cm.hot, interpolation='bilinear')
cbar = plt.colorbar(im)
#plt.plot(listx,recta_para_lambda_p1,'-',color = 'green')
if paper == 0:

    plt.plot(listx,recta_para_lambda_p2,'-',color = 'green')
#plt.plot(listx,recta_para_lambda_p3,'-',color = 'green')

#plt.plot(listx,recta_from_plot_p1,'-',color = 'blue')
    plt.plot(listx,recta_from_plot_p2,'-',color = 'blue')
#plt.plot(listx,recta_from_plot_p3,'-',color = 'blue')

cbar.ax.tick_params(labelsize = tamnum)
#cbar.set_ticks(tick_locations)

plt.ylim([listy[0],listy[-1]])
plt.xlim([listx[0],listx[-1]])

cbar.set_label(r'Im{$\phi_%i$}' %(Nmax0))
os.chdir(path_save)
plt.tight_layout()
plt.savefig('Im_potential_ana' + label1 + '.png', format='png')

#%%
    

#def disp_relation_omega(kx,a,Nmax,epsi1,epsi2,hbmu):
#
#    """
#    Parameters
#    ----------
#    omega : frecuencia en Hz
#    mu_c : chemical potential of graphene in eV
#    gamma_in : frecuencia de colision en eV
#    epsilon1 : permeabilidad electrica del medio de arriba
#    epsilon2 : permeabilidad electrica del medio de abajo
#    Returns
#        relacion de dispersion para un plano de Ag
#        (solucion analitica)
#    -------
#    """
#
#    cte = (epsi1 + epsi2)
#    
#    kn = kx + 2*np.pi*Nmax/a
#    
#    
#    return np.sqrt(kn*4*hbmu*alfac*(c/hb)/cte)

#if plot_check_position_disp_relation == 1:
#    
#    list_colours = ['black', 'gold', 'orange', 'darksalmon', 'brown', 'darkred']
#    title1_aux = '$\epsilon_1 = %.1f$, $\epsilon_2 = %.1f$, $\gamma_{in}$ = %.2f meV' %(epsi1,epsi2,hbgama*1e3)
#    title2_aux = '$\mu_c$ = %.2f eV, a = %i nm, n = %i' %(hbmu,a0*1e3,Nmax0)
#    
#    
#    Ntot = 200
#    list_kx = np.linspace(0.001, 4 + Nmax0/(a0), Ntot) ##creo que va a estar en 1/micrones
#    
#    list_y_re = []
#    list_y_im = []
#    for kx in list_kx:
#        
#        valuey = disp_relation_omega(kx,a0,Nmax0,epsi1,epsi2,hbmu)
#        valuey = hb*valuey*1e3
#        list_y_re.append(valuey.real)
#        list_y_im.append(valuey.imag)
#        
#    
#    plt.figure(figsize=tamfig)    
#    plt.title(title1_aux + '\n' + title2_aux ,fontsize=tamtitle)
#    plt.plot(list_kx,list_y_re,'-',color = 'darkgreen',ms = lw, label = 'DR')
#    list_hb_omega_meV = np.linspace(np.min(list_y_re),np.max(list_y_re),Ntot)
#    j = 0
#    list_cortes_omega = []
#    list_cortes_k = []
#    list_int_v = [1/int_v0]
#    
#    for int_v in list_int_v:
#        v = c*int_v
#        
#        list_y2_re = []
#        
#        for hb_omega_meV in list_hb_omega_meV:
#            omega = hb_omega_meV*1e-3/hb
#            
#            valuey_e = omega/v
#            list_y2_re.append(valuey_e.real)
#            
#        
#        hb_omega_corte = omega_n_corte(1/int_v,a_nm,Nmax)*hb*1e3
#        
#        k_corte = omega_n_corte(1/int_v,a_nm,Nmax)/v
#        
#        list_cortes_omega.append(hb_omega_corte)
#        list_cortes_k.append(k_corte)
#            
#        ##  
#        if int_v != 1:
#            plt.plot(list_y2_re,list_hb_omega_meV,'-', color = list_colours[j],ms = lw, label = 'v= %.1fc' %(int_v))
#            plt.plot(list_cortes_k,list_cortes_omega,'.', color = 'red',ms = '3')
#        else:
#            plt.plot(list_y2_re,list_hb_omega_meV,'-', color = list_colours[j], ms = lw, label = 'light')
#            
#        j = j + 1
#         
#    
#    plt.xlabel('$k_\parallel$ [1/$\mu$m]',fontsize=tamletra, labelpad = 0)
#    plt.ylabel('$\hbar\omega$ [meV]',fontsize=tamletra, labelpad = 0)
#    plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
#    plt.tick_params(labelsize = tamnum,pad = pad)
#    plt.grid(1)
#    plt.tight_layout()
#    os.chdir(path_save)
#    plt.savefig('disp_relation_for_potential' + label1 + '.png', format='png')



#%%
