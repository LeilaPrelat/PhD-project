
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
import seaborn as sns

#%%

save_graphs = 1 #guardar los graficos 2D del campo
plot_sin_QE = 1

plot_xy_3D = 0 # mapa de color del campo vs x,y
plot_vs_omegaTHz = 1 # graficar vs omega THz

sns.set()    

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/External_Efield','')
path_save = path_basic + '/' + 'fieldE_direct'

if save_graphs==1:    
    if not os.path.exists(path_save):
        print('Creating folder to save graphs')
        os.mkdir(path_save)

err = 'fieldE_direct_numerical.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from fieldE_direct_numerical2_2 import Efield_NUM_QE,Efield_NUM
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

int_v = 400
#v = c/den
epsi1 = 1

b = -2 
z = 0

if plot_xy_3D ==1 :
    cota = 2
    omegaTHz = 0.01
    omegac0 = omegaTHz*aux2   
    #omega = 0.7*1e12
    title = r'$\epsilon_1$ = %i, v = c/%i $\mu$m/s, $\omega$ = %.2fTHz' %(epsi1,int_v,omegaTHz) 
    label1 = '_b%i_3D' %(b)
    labelx,labely = 'x [$\mu$m]', 'y [$\mu$m]'

elif plot_vs_omegaTHz == 1:
    x0 = 0
    y0 = 0
    list_OmegaTHz = np.linspace(0.01,2.01,251)
    list_OmegaTHz = np.linspace(0.001,0.2,200)
    label1 = '_b%i_vsOmega' %(b)
    title = r'$\epsilon_1$ = %i, v = c/%i $\mu$m/s, x = %i$\mu$m, y = %i$\mu$m, ' %(epsi1,int_v,x0,y0) 
    labelx = r'$\omega$ [THz]'
    
title = title + '\n' + r'b = %.2f$\mu$m, z = %i$\mu$m'%(b,z) + ', ' +  name_this_py 

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

# color_rgb1 = (1.0, 0.8, 0.6)
# color_rgb2 = (0.99, 0.74, 0.71)
# color_rgb3 = (0.529, 0.992, 0.945)
# color_rgb4 = (0.41, 0.21, 0.61)

color_rgb1 = (1.0, 0.01, 0.24)
color_rgb2 = (0.8, 0.58, 0.46)
color_rgb3 = (0.55, 0.71, 0)
color_rgb4 = (0.6, 0.73, 0.89)

if plot_xy_3D ==1 : # mapa de color del campo vs x,y

    N = 100
    list_x = np.linspace(-cota,cota,N)
    list_y = np.linspace(-cota,cota,N)
    X, Y = np.meshgrid(list_x, list_y)
    limits = [min(list_x) , max(list_x), min(list_y) , max(list_y)]
    
    def functionQE_real(x,y):
        rta = Efield_NUM_QE(omegac0,epsi1,int_v,b,x,y,z)
        return rta.real
    
    def function_real(x,y):
        rta = Efield_NUM(omegac0,epsi1,int_v,b,x,y,z)
        return rta.real
    
    f1_real = np.vectorize(functionQE_real)
    Z1_real = f1_real(X, Y)
    
    f2_real = np.vectorize(function_real)
    Z2_real = f2_real(X, Y)
    
    def functionQE_imag(x,y):
        rta = Efield_NUM_QE(omegac0,epsi1,int_v,b,x,y,z)
        return rta.imag
    
    def function_imag(x,y):
        rta = Efield_NUM(omegac0,epsi1,int_v,b,x,y,z)
        return rta.imag
    
    f1_imag = np.vectorize(functionQE_imag)
    Z1_imag = f1_imag(X, Y)
    
    f2_imag = np.vectorize(function_imag)
    Z2_imag = f2_imag(X, Y)
    
    os.chdir(path_save)
    
    np.savetxt('X_cota%i' %(cota) + label1 + '.txt', X, fmt='%1.11e', delimiter='\t', header = title)
    np.savetxt('Y_cota%i' %(cota) + label1 + '.txt', X, fmt='%1.11e', delimiter='\t', header = title)
    
    if plot_sin_QE == 0:
        header1 = 'Re E_{dir,x} QE num' + ', ' +  title
        # tabla1 = np.array(Z1)
        np.savetxt('Re_Edir_QE_cota%i' %(cota) + label1 + '.txt', Z1_real, fmt='%1.11e', delimiter='\t', header = header1)
    else:
        header1 = 'Re E_{dir,x} QE num' + ', ' +  title
        # tabla1 = np.array(Z1)
        np.savetxt('Re_Edir_QE_cota%i' %(cota) + label1 + '.txt', Z1_real, fmt='%1.11e', delimiter='\t', header = header1)
        
        header2 = 'Re E_{dir,x} num' + ', ' +  title
        # tabla2 = np.array([Z2])    
        np.savetxt('Re_Edir_cota%i' %(cota) + label1 + '.txt', Z2_real, fmt='%1.11e', delimiter='\t', header = header2)    
     
    #####################################
    
    if plot_sin_QE == 0:
        header1 = 'Im E_{dir,x} QE num' + ', ' +  title
        # tabla1 = np.array(Z1)
        np.savetxt('Im_Edir_QE_cota%i' %(cota) + label1 + '.txt', Z1_imag, fmt='%1.11e', delimiter='\t', header = header1)
    else:
        header1 = 'Im E_{dir,x} QE num' + ', ' +  title
        # tabla1 = np.array(Z1)
        np.savetxt('Im_Edir_QE_cota%i' %(cota) + label1 + '.txt', Z1_imag, fmt='%1.11e', delimiter='\t', header = header1)
        
        header2 = 'Im E_{dir,x} num' + ', ' +  title
        # tabla2 = np.array([Z2])    
        np.savetxt('Im_Edir_cota%i' %(cota) + label1 + '.txt', Z2_imag, fmt='%1.11e', delimiter='\t', header = header2)   
       
    print('Graficar el External field' + label1)
    
    graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    im = plt.imshow(Z1_real, extent = limits, cmap=plt.cm.hot, interpolation='bilinear')
    cbar = plt.colorbar(im)
    cbar.ax.tick_params(labelsize = tamnum)
    if save_graphs == True:
        plt.savefig('Re_Edir_QE_cota%i' %(cota) + label1 + '.png', format='png')
    
    
    graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    im = plt.imshow(Z2_real, extent = limits, cmap=plt.cm.hot, interpolation='bilinear')
    cbar = plt.colorbar(im)
    cbar.ax.tick_params(labelsize = tamnum)
    if save_graphs == True:
        plt.savefig('Re_Edir_cota%i' %(cota) + label1 + '.png', format='png')
        
    
    graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    im = plt.imshow(Z1_imag, extent = limits, cmap=plt.cm.hot, interpolation='bilinear')
    cbar = plt.colorbar(im)
    cbar.ax.tick_params(labelsize = tamnum)
    if save_graphs == True:
        plt.savefig('Imag_Edir_QE_cota%i' %(cota) + label1 + '.png', format='png')
    
    
    graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    im = plt.imshow(Z2_imag, extent = limits, cmap=plt.cm.hot, interpolation='bilinear')
    cbar = plt.colorbar(im)
    cbar.ax.tick_params(labelsize = tamnum)
    if save_graphs == True:
        plt.savefig('Imag_Edir_cota%i' %(cota) + label1 + '.png', format='png')    

#%%

elif plot_vs_omegaTHz == 1:
    
    listyQE_re = []
    listy_re = []
    listyQE_im = []
    listy_im = []
    
    j = 0
    for omegaTHz in list_OmegaTHz:
        omegaTHz = np.round(omegaTHz,8)
        omegac = omegaTHz*aux2  

        rta = Efield_NUM_QE(omegac,epsi1,int_v,b,x0,y0,z)
        listyQE_re.append(rta.real)
        listyQE_im.append(rta.imag)
        
        if plot_sin_QE == 1:
            rta = Efield_NUM(omegac,epsi1,int_v,b,x0,y0,z)
            listy_re.append(rta.real)
            listy_im.append(rta.imag)
        
        j = j + 1
        print(j)
    
    graph(title,labelx,r'Re$(E_{dir,x})$',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    plt.plot(list_OmegaTHz,listyQE_re,'.',ms=ms,color = 'lightseagreen',label = 'numerical QE')
    if plot_sin_QE == 1:
        plt.plot(list_OmegaTHz,listy_re,'--',lw = 2,color = 'darkviolet',label = 'numerical')
    # plt.yscale('log')
    plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
    if save_graphs==1:
        plt.tight_layout()
        os.chdir(path_save)
        plt.savefig( 'Re_Edir_' + label1 + '.png', format='png')
    
    graph(title,labelx,r'Im$(E_{dir,x})$',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    plt.plot(list_OmegaTHz,listyQE_im,'.',ms = ms,color = 'lightseagreen',label = 'numerical QE')
    if plot_sin_QE == 1:
        plt.plot(list_OmegaTHz,listy_im,'--',lw = 2,color = 'darkviolet',label = 'numerical')
    # plt.yscale('log')
    plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
    if save_graphs==1:
        plt.tight_layout()
        plt.savefig( 'Im_Edir_' + label1 + '.png', format='png')    

    
    np.savetxt('list_OmegaTHz' + label1 + '.txt', list_OmegaTHz, fmt='%1.11e', delimiter='\t', header = title)
    if plot_sin_QE == 0:
        header1 = 'E_{dir,x} QE num' + ', ' +  title
        np.savetxt('Re_Edir_QE' + label1 + '.txt', listyQE_re, fmt='%1.11e', delimiter='\t', header = header1)
        np.savetxt('Im_Edir_QE'+ label1 + '.txt', listyQE_im, fmt='%1.11e', delimiter='\t', header = header1)
        
    else:
        
        header1 = 'E_{dir,x} QE num' + ', ' +  title
        np.savetxt('Re_Edir_QE'  + label1 + '.txt', listyQE_re, fmt='%1.11e', delimiter='\t', header = header1)
        np.savetxt('Im_Edir_QE'  + label1 + '.txt', listyQE_im, fmt='%1.11e', delimiter='\t', header = header1)
        
        header2 = 'E_{dir,x} num' + ', ' +  title
        np.savetxt('Re_Edir'  + label1 + '.txt', listy_re, fmt='%1.11e', delimiter='\t', header = header2)  
        np.savetxt('Im_Edir'  + label1 + '.txt', listy_im, fmt='%1.11e', delimiter='\t', header = header2)         
        
#%%