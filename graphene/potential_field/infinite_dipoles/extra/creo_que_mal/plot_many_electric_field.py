
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
from scipy.signal import find_peaks

#sys.setrecursionlimit(50000)

#%%

save_graphs = 1 #guardar los graficos 2D del campo

plot_Ex = 0
plot_Etot = 1

plot_vs_xy = 1
plot_vs_yz = 0

plot_vs_Energy = 0
plot_num = 1

periodicity = 1

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/potential_field/many_potential_Javier_formula','')


err = 'electric_field_many_dipoles.py no se encuentra en ' + path_basic
path_save = path_basic + '/' + 'many_electric_field_f'

if save_graphs==1:    
    if not os.path.exists(path_save):
        print('Creating folder to save graphs')
        os.mkdir(path_save)
try:
    sys.path.insert(1, path_basic)
    from electric_field_many_dipoles import phi_many_dipoles_num
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

int_v = 350
#v = c/int_v
#omega = 0.7*1e12
cota = 1200 #nanometros
epsi1,epsi2 = 1,1
hbmu,hbgama = 0.3,0.0001
zp = 0.05 #micrometros 

px = 1
py = 0
pz = 0

omega0THz = 0.7
omega0 = omega0THz*1e12
lambdda0 = c/(2*np.pi*omega0)
E0 =omega0THz*1e12*hb  #eV
omegac0 = E0/aux   

N = 10
a = 5 


title2 = r'v = c/%i $\mu$m/s, $\hbar\mu$ = %.2feV, $\hbar\gamma$ = %.4feV' %(int_v,hbmu,hbgama) 
title4 = r'$z_p$=%inm, px = %i, py = %i, pz = %i' %(zp*1e3,px,py,pz)
labelp = r'_px%i_py%i_pz%i' %(px,py,pz)
labelx = r'E [eV]'

list_xD = np.linspace(-1500,1500,11)
# list_xD = np.array([10])

if plot_num == 1:
    if plot_vs_yz == 1:
        
        x0 = 0
        label1 = r'_yz_3D_num' + labelp  + '_E%.4f' %(E0)
        labelx,labely = 'y [nm]', 'z [nm]'
        title1 = r'a = %inm, N = %i, E = %.2fmeV, x = %inm' %(a*1e3,N,E0*1e3,x0*1e3)     
        def function_Ex_real(y_nano,z_nano):
            
            y = y_nano*1e-3 
            z = z_nano*1e-3
            
            Ex,Ey,Ez = phi_many_dipoles_num(omegac0,epsi1,epsi2,hbmu,hbgama,x0,y,z,a,N,zp,int_v,px,py,pz)
            return Ex.real
    
    
        def function_Ex_imag(y_nano,z_nano):
            
            y = y_nano*1e-3 
            z = z_nano*1e-3
            
            Ex,Ey,Ez = phi_many_dipoles_num(omegac0,epsi1,epsi2,hbmu,hbgama,x0,y,z,a,N,zp,int_v,px,py,pz)
            return Ex.imag
        
        def function_Etot(y_nano,z_nano):
            y = y_nano*1e-3 
            z = z_nano*1e-3
                   
            
            Ex,Ey,Ez = phi_many_dipoles_num(omegac0,epsi1,epsi2,hbmu,hbgama,x0,y,z,a,N,zp,int_v,px,py,pz)
            rta = np.abs(Ex)**2 + np.abs(Ey)**2 + np.abs(Ez)**2 
            return rta
    
    elif plot_vs_xy == 1:
        
        z0 = zp
        label1 = r'_xy_3D_num' + labelp  + '_E%.4f' %(E0)
        labelx,labely = 'x [nm]', 'y [nm]'
        title1 = r'a = %inm, N = %i, E = %.2fmeV, z = %inm' %(a*1e3,N,E0*1e3,z0*1e3)    
        
        def function_Ex_real(x_nano,y_nano):
            
            x = x_nano*1e-3 
            y = y_nano*1e-3
            
            Ex,Ey,Ez = phi_many_dipoles_num(omegac0,epsi1,epsi2,hbmu,hbgama,x,y,z0,a,N,zp,int_v,px,py,pz)
            return Ex.real
    
    
        def function_Ex_imag(x_nano,y_nano):
            
            x = x_nano*1e-3 
            y = y_nano*1e-3
            
            Ex,Ey,Ez = phi_many_dipoles_num(omegac0,epsi1,epsi2,hbmu,hbgama,x,y,z0,a,N,zp,int_v,px,py,pz)
            return Ex.imag
        
        def function_Etot(x_nano,y_nano):
            x = x_nano*1e-3 
            y = y_nano*1e-3
                   
            
            Ex,Ey,Ez = phi_many_dipoles_num(omegac0,epsi1,epsi2,hbmu,hbgama,x,y,z0,a,N,zp,int_v,px,py,pz)
            rta = np.abs(Ex)**2 + np.abs(Ey)**2 + np.abs(Ez)**2 
            return rta
    
    
    
    elif plot_vs_Energy == 1:
        x0,y0,z0 = 0,0,0
        title1 = r'v = c/%i $\mu$m/s, a = %inm, N = %i, y = %i, z = %i' %(int_v,a*1e3,N,x0,y0,z0) 
        label1 = '_vs_E_num' + labelp
        labelx = 'E [eV]'
        list_OmegaTHz = np.linspace(0.01,2.01,251)    
        listE = np.linspace(0.001, 10, 101)
    
    
        def function_Ex_real(energy0):
            omegac0 = energy0/aux           
            Ex,Ey,Ez = phi_many_dipoles_num(omegac0,epsi1,epsi2,hbmu,hbgama,x0,y0,z0,a,N,zp,int_v,px,py,pz)
            return Ex.real
    
        def function_Ex_imag(energy0):
            omegac0 = energy0/aux           
            Ex,Ey,Ez = phi_many_dipoles_num(omegac0,epsi1,epsi2,hbmu,hbgama,x0,y0,z0,a,N,zp,int_v,px,py,pz)
            return Ex.imag
        
        def function_Etot(energy0):
            omegac0 = energy0/aux           
            Ex,Ey,Ez = phi_many_dipoles_num(omegac0,epsi1,epsi2,hbmu,hbgama,x0,y0,z0,a,N,zp,int_v,px,py,pz)
            rta = np.abs(Ex)**2 + np.abs(Ey)**2 + np.abs(Ez)**2 
            return rta 

title = title1 + '\n' + title2 + '\n' + title4 

#%%
    
tamfig = (4.5,3.5)
tamlegend = 10
tamletra = 11
tamtitle = 10
tamnum = 9
labelpady = -1.5
labelpadx = 2
pad = 0
mk = 2
ms = 4
hp = 0.3

def graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad):
    plt.figure(figsize=tamfig)
    plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
    plt.ylabel(labely,fontsize=tamletra,labelpad =labelpadx)
    plt.tick_params(labelsize = tamnum, pad = pad)
    plt.title(title,fontsize=int(tamtitle*0.9))

    return   

#%%
    
if plot_Ex == 1:
    if plot_vs_Energy == 0:
        Nlist = 100
        list_x = np.linspace(-cota,cota,Nlist)
        list_y = np.linspace(-cota,cota,Nlist)
        X, Y = np.meshgrid(list_x, list_y)
        limits = [min(list_x) , max(list_x), min(list_y) , max(list_y)]
        
        aux_list = np.array(np.ones(Nlist))
        
        f1_real = np.vectorize(function_Ex_real)
        Z1_real = f1_real(X, Y)
        
        f2_real = np.vectorize(function_Ex_imag)
        Z2_real = f2_real(X, Y)
           
        print('Graficar el External field' + label1)
        
        graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
        im = plt.imshow(Z1_real, extent = limits, cmap='RdBu_r', interpolation='bilinear')
        cbar = plt.colorbar(im)
        cbar.ax.tick_params(labelsize = tamnum)
        cbar.set_label(r'Re($E_x$)',fontsize=tamlegend,labelpad = 1)
        if save_graphs == True:
            os.chdir(path_save)
            plt.tight_layout()
            plt.savefig('Re_many_Ex' + label1 + '.png', format='png')
            
        
        graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
        im = plt.imshow(Z2_real, extent = limits, cmap='RdBu_r', interpolation='bilinear')
        cbar = plt.colorbar(im)
        cbar.ax.tick_params(labelsize = tamnum)
        
        cbar.set_label(r'Im($E_x$)',fontsize=tamlegend,labelpad = 1)
        if save_graphs == True:
            plt.tight_layout()
            plt.savefig('Im_many_Ex' + label1 + '.png', format='png')

 
elif plot_Etot == 1:
    if plot_vs_Energy == 0:
        Nlist = 100
        list_x = np.linspace(-cota,cota,Nlist)
        list_y = np.linspace(-cota,cota,Nlist)
        X, Y = np.meshgrid(list_x, list_y)
        limits = [min(list_x) , max(list_x), min(list_y) , max(list_y)]
        
        aux_list = np.array(np.ones(Nlist))
        
        f1_real = np.vectorize(function_Etot)
        Z1_real = f1_real(X, Y)
           
        print('Graficar el External field' + label1)
        
        graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
        im = plt.imshow(Z1_real, extent = limits, cmap='RdBu_r', interpolation='bilinear')
        cbar = plt.colorbar(im)
        cbar.ax.tick_params(labelsize = tamnum)
        cbar.set_label(r'$E_{tot}$',fontsize=tamlegend,labelpad = 1)

        yD, zD = 0,0
        for xD in list_xD:
            if plot_vs_xy == 1:
                plt.plot(xD*aux_list, list_y,'--',lw = 1,color = 'green')
                plt.plot(list_x, yD*aux_list,'--',lw = 1,color = 'green')
            elif plot_vs_yz == 1:    
                plt.plot(yD*aux_list, list_y,'--',lw = 1,color = 'green')
                plt.plot(list_x, zD*aux_list,'--',lw = 1,color = 'green')
        if save_graphs == True:
            os.chdir(path_save)
            plt.tight_layout()
            plt.savefig('Re_many_Etot' + label1 + '.png', format='png')
                 
        
        
        
#else:
#    import seaborn as sns
#    sns.set()
#    
#    listy_re = []
#    listy_im = []
#    j = 0
#    for E in listE:
#        rta_re = function_real(E)
#        rta_im = function_imag(E)
#        listy_re.append(np.log10(np.abs(rta_re)))
#        if pz == 1:
#            listy_im.append(np.log10(np.abs(rta_im)))       
#        else:
#            listy_im.append(rta_im) 
#        print(j)
#        j = j + 1
#        
#    graph(title,labelx,r'log Re($\phi$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#    plt.plot(listE,-np.array(listy_re),'.-',ms = ms,color = 'lightseagreen')
#    
#    # plt.yscale('log')
##    plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
#    if save_graphs==1:
#        plt.tight_layout()
#        os.chdir(path_save)
#        plt.savefig( 'Re_many_potential' + label1 + '.png', format='png')    
#
#    if pz == 1:
#        graph(title,labelx,r'log Im($\phi$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#        plt.plot(listE,-np.array(listy_im),'.-',ms = ms,color = 'lightseagreen')    
#    else:
#        graph(title,labelx,r'Im($\phi$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#        plt.plot(listE,listy_im,'.-',ms = ms,color = 'lightseagreen')
#        
#        # plt.yscale('log')
##    plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
#    if save_graphs==1:
#        plt.tight_layout()
#        os.chdir(path_save)
#        plt.savefig( 'Im_many_potential' + label1 + '.png', format='png')    

    
#%%