
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
import matplotlib.colors as colors
from scipy.signal import find_peaks
#sys.setrecursionlimit(50000)

#%%

save_graphs = 1 #guardar los graficos 2D del campo

plot_vs_xy = 0
plot_vs_yz = 0
plot_vs_Energy = 0

periodicity = 1

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/potential_field/many_potential','')
path_save = path_basic + '/' + 'many_potential_f'

if save_graphs==1:    
    if not os.path.exists(path_save):
        print('Creating folder to save graphs')
        os.mkdir(path_save)

err = 'many_potential.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from potential_many_dipoles import phi_many_dipoles_num
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

int_v = 600
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

N = 3
a = 0.5 
yD = 0
zD = 0

title2 = r'v = c/%i $\mu$m/s, $\hbar\mu$ = %.2feV, $\hbar\gamma$ = %.4feV' %(int_v,hbmu,hbgama) 
title4 = r'$z_p$=%inm, px = %i, py = %i, pz = %i' %(zp*1e3,px,py,pz)
labelp = r'_px%i_py%i_pz%i' %(px,py,pz)

list_xD = np.linspace(-a*1e3*(N-1),a*1e3*(N-1),2*N-1)
# list_xD = np.array([10])


if plot_vs_yz == 1:
    
    x0 = 0
    label1 = r'_yz_3D' + labelp  + '_E%.4f' %(E0)
    labelx,labely = 'y [nm]', 'z [nm]'
    title1 = r'a = %inm, N = %i, E = %.2fmeV, x = %inm' %(a*1e3,N,E0*1e3,x0*1e3)     
    def function_real(y_nano,z_nano):
        
        y = y_nano*1e-3 
        z = z_nano*1e-3
        
        rta = phi_many_dipoles_num(omegac0,epsi1,epsi2,hbmu,hbgama,x0,y,z,a,N,zp,int_v,px,py,pz)
        return rta.real
    
    def function_imag(y_nano,z_nano):
        y = y_nano*1e-3 
        z = z_nano*1e-3
               
        
        rta = phi_many_dipoles_num(omegac0,epsi1,epsi2,hbmu,hbgama,x0,y,z,a,N,zp,int_v,px,py,pz)
        return rta.imag

elif plot_vs_xy == 1:
    
    z0 = zp
    label1 = r'_xy_3D' + labelp  + '_E%.4f' %(E0)
    labelx,labely = 'x [nm]', 'y [nm]'
    title1 = r'a = %inm, N = %i, E = %.2fmeV, z = %inm' %(a*1e3,N,E0*1e3,z0*1e3)    
    
    def function_real(x_nano,y_nano):
        
        x = x_nano*1e-3 
        y = y_nano*1e-3
        
        rta = phi_many_dipoles_num(omegac0,epsi1,epsi2,hbmu,hbgama,x,y,z0,a,N,zp,int_v,px,py,pz)
        return rta.real
    
    def function_imag(x_nano,y_nano):
        x = x_nano*1e-3 
        y = y_nano*1e-3
               
        
        rta = phi_many_dipoles_num(omegac0,epsi1,epsi2,hbmu,hbgama,x,y,z0,a,N,zp,int_v,px,py,pz)
        return rta.imag



elif plot_vs_Energy == 1:
    x0,y0,z0 = 0,0,0
    title1 = r'v = c/%i $\mu$m/s, a = %inm, N = %i, y = %i, z = %i' %(int_v,a*1e3,N,x0,y0,z0) 
    label1 = '_vs_E' + labelp
    labelx = 'E [eV]'
    list_OmegaTHz = np.linspace(0.01,2.01,251)    
    listE = np.linspace(0.001, 10, 101)


    def function_real(energy0):
        omegac0 = energy0/aux           
        rta = phi_many_dipoles_num(omegac0,epsi1,epsi2,hbmu,hbgama,x0,y0,z0,a,N,zp,int_v,px,py,pz)
        return rta.real
    
    def function_imag(energy0):
        omegac0 = energy0/aux           
        rta = phi_many_dipoles_num(omegac0,epsi1,epsi2,hbmu,hbgama,x0,y0,z0,a,N,zp,int_v,px,py,pz)
        return rta.imag

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
if plot_vs_Energy == 0:
    Nlist = 100
    list_x = np.linspace(-cota,cota,Nlist)
    list_y = np.linspace(-cota,cota,Nlist)
    X, Y = np.meshgrid(list_x, list_y)
    limits = [min(list_x) , max(list_x), min(list_y) , max(list_y)]
    
    aux_list = np.array(np.ones(Nlist))
    
    f1_real = np.vectorize(function_real)
    Z1_real = f1_real(X, Y)
    
    f2_real = np.vectorize(function_imag)
    Z2_real = f2_real(X, Y)
       
    print('Graficar el External field' + label1)
    
    graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    im = plt.imshow(Z1_real, extent = limits, cmap='RdBu_r', interpolation='bilinear')
    cbar = plt.colorbar(im)
    cbar.ax.tick_params(labelsize = tamnum)
    cbar.set_label(r'Re($\phi$)',fontsize=tamlegend,labelpad = 1)
    
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
        plt.savefig('Re_many_potential' + label1 + '.png', format='png')
        
    ####################################################################################
#    graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#    vmin, vmax = np.min(Z1_real), np.max(Z1_real)
#    pcm = plt.pcolormesh(X, Y, Z1_real,
#                          norm=colors.SymLogNorm(linthresh=0.03, linscale=0.03,vmin=int(vmin), vmax=int(vmax)),cmap='RdBu_r')
#
#    maxlog=int(np.ceil( np.log10( np.abs(vmax) )))
#    minlog=int(np.ceil( np.log10( np.abs(vmin) )))
#    
#    if vmin < 0 :
#        tick_locations = ( [-(10.0**x) for x in np.linspace(minlog,-1,minlog+2)] 
#                            + [0] 
#                            + [(10.0**x) for x in np.linspace(-1,maxlog,maxlog+2)] )
#    else:
#        tick_locations = ( [(10.0**x) for x in np.linspace(minlog,maxlog,maxlog - np.abs(minlog) + 1) ]) 
#    cbar = plt.colorbar(pcm)
#     cbar.set_ticks(tick_locations)
#    cbar.ax.tick_params(labelsize = tamnum)
#    cbar.set_label(r'Re($\phi$)',fontsize=tamlegend,labelpad = 1)
#    if save_graphs == True:
#        os.chdir(path_save)
#        plt.tight_layout()
#        plt.savefig('Re_many_potential_log' + label1 + '.png', format='png')
    ###################################################################################
    
    graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    im = plt.imshow(Z2_real, extent = limits, cmap='RdBu_r', interpolation='bilinear')
    cbar = plt.colorbar(im)
    cbar.ax.tick_params(labelsize = tamnum)
    for xD in list_xD:
        if plot_vs_xy == 1:
            plt.plot(xD*aux_list, list_y,'--',lw = 1,color = 'green')
            plt.plot(list_x, yD*aux_list,'--',lw = 1,color = 'green')
        elif plot_vs_yz == 1:    
            plt.plot(yD*aux_list, list_y,'--',lw = 1,color = 'green')
            plt.plot(list_x, zD*aux_list,'--',lw = 1,color = 'green')    
    cbar.set_label(r'Im($\phi$)',fontsize=tamlegend,labelpad = 1)
    if save_graphs == True:
        plt.tight_layout()
        plt.savefig('Im_many_potential' + label1 + '.png', format='png')

    ###################################################################################
#    graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#    vmin, vmax = np.min(Z2_real), np.max(Z2_real)
#    pcm = plt.pcolormesh(X, Y, Z2_real,
#                          norm=colors.SymLogNorm(linthresh=0.03, linscale=0.03,vmin=int(vmin), vmax=int(vmax)),cmap='RdBu_r')
#
#    maxlog=int(np.ceil( np.log10( np.abs(vmax) )))
#    minlog=int(np.ceil( np.log10( np.abs(vmin) )))
#    
#    if vmin < 0 :
#        tick_locations = ( [-(10.0**x) for x in np.linspace(minlog,-1,minlog+2)] 
#                            + [0] 
#                            + [(10.0**x) for x in np.linspace(-1,maxlog,maxlog+2)] )
#    else:
#        tick_locations = ( [(10.0**x) for x in np.linspace(minlog,maxlog,maxlog - np.abs(minlog) + 1) ])                                               
#    cbar = plt.colorbar(pcm)
#    cbar.set_ticks(tick_locations)
#    cbar.ax.tick_params(labelsize = tamnum)
#    cbar.set_label(r'Im($\phi$)',fontsize=tamlegend,labelpad = 1)
#    if save_graphs == True:
#        os.chdir(path_save)
#        plt.tight_layout()
#        plt.savefig('Im_many_potential_log' + label1 + '.png', format='png')
    
else:
    import seaborn as sns
    sns.set()
    
    listy_re = []
    listy_im = []
    j = 0
    for E in listE:
        rta_re = function_real(E)
        rta_im = function_imag(E)
        listy_re.append(np.log10(np.abs(rta_re)))
        if pz == 1:
            listy_im.append(np.log10(np.abs(rta_im)))       
        else:
            listy_im.append(rta_im) 
        print(j)
        j = j + 1
        
    graph(title,labelx,r'log Re($\phi$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    plt.plot(listE,-np.array(listy_re),'.-',ms = ms,color = 'lightseagreen')
    
    # plt.yscale('log')
#    plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
    if save_graphs==1:
        plt.tight_layout()
        os.chdir(path_save)
        plt.savefig( 'Re_many_potential' + label1 + '.png', format='png')    

    if pz == 1:
        graph(title,labelx,r'log Im($\phi$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
        plt.plot(listE,-np.array(listy_im),'.-',ms = ms,color = 'lightseagreen')    
    else:
        graph(title,labelx,r'Im($\phi$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
        plt.plot(listE,listy_im,'.-',ms = ms,color = 'lightseagreen')
        
        # plt.yscale('log')
#    plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
    if save_graphs==1:
        plt.tight_layout()
        os.chdir(path_save)
        plt.savefig( 'Im_many_potential' + label1 + '.png', format='png')    



if periodicity == 1:
    Nlist = 100
    list_x = np.linspace(-cota,cota,Nlist)
    
    def function_Etot1D_re(x_nano):
        x = x_nano*1e-3 
        y = 0
        
        rta = phi_many_dipoles_num(omegac0,epsi1,epsi2,hbmu,hbgama,x,y,z0,a,N,zp,int_v,px,py,pz)
        return rta.real

    def function_Etot1D_im(x_nano):
        x = x_nano*1e-3 
        y = 0
        
        rta = phi_many_dipoles_num(omegac0,epsi1,epsi2,hbmu,hbgama,x,y,z0,a,N,zp,int_v,px,py,pz)
        return rta.imag


    list_y_re = []
    for x in list_x:  
        y = function_Etot1D_re(x)
        list_y_re.append(x)
            
    graph(title,labelx,r'Re($\phi$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    plt.plot(list_x,list_y_re,'.-',ms = ms,color = 'lightseagreen')   
    
    peaks, _ = find_peaks(list_y_re, height=0)
    plt.plot(peaks, list_y[peaks], "x")
    
    list_minimus = []
    for j in range(len(peaks)-1):
        minimum = np.abs(peaks[j] - peaks[j + 1])
        list_minimus.append(minimum)
    
    final_minimum = np.mean(list_minimus)
    np.savetxt('Etot1D_re' + label1 + '.txt',[final_minimum], fmt='%.18e', delimiter='\n',header = title)
    
    if save_graphs==1:
        plt.tight_layout()
        os.chdir(path_save)
        plt.savefig( 'Etot1D_re' + label1 + '.png', format='png')    


    list_y_im = []
    for x in list_x:  
        y = function_Etot1D_im(x)
        list_y_im.append(x)

    graph(title,labelx,r'Im($\phi$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    plt.plot(list_x,list_y_im,'.-',ms = ms,color = 'lightseagreen')   
    
    peaks, _ = find_peaks(list_y_im, height=0)
    plt.plot(peaks, list_y[peaks], "x")
    
    list_minimus = []
    for j in range(len(peaks)-1):
        minimum = np.abs(peaks[j] - peaks[j + 1])
        list_minimus.append(minimum)
    
    final_minimum = np.mean(list_minimus)
    np.savetxt('Etot1D_im' + label1 + '.txt',[final_minimum], fmt='%.18e', delimiter='\n',header = title)
    
    if save_graphs==1:
        plt.tight_layout()
        os.chdir(path_save)
        plt.savefig( 'Etot1D_im' + label1 + '.png', format='png')    
       
#%%
