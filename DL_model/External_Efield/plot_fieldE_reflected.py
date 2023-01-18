
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: leila

graficar el campo externo reflejado en el plano con la convencion de z hacia abajo
en z = 0
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
plot_vs_z =  0# graficar vs z (desde antes del electron hasta despues del plano)
sns.set()    

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/External_Efield','')
path_save = path_basic + '/' + 'fieldE_ref'

if save_graphs==1:    
    if not os.path.exists(path_save):
        print('Creating folder to save graphs')
        os.mkdir(path_save)

err = 'fieldE_ref_numerical.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from fieldE_ref_numerical import Efield_NUM_QE,Efield_NUM,Efield_NUM_SP,Efield_ana1,Efield_ana2,Efield_ana1b
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
v = c/int_v
b = -0.05 #electron position in z < 0 (arriba del plano)
#omega = 0.7*1e12
cota = 2
epsi1,epsi2 = 1,1
hbgama = 0.0001
zp = 0.05

epsilon_b = 4
d =  0.03   #micrones, 10 nm
Ebulk = 9.17
omega_bulk = Ebulk/hb

title4 = r'$\epsilon_b$ = %i, d = %.4f$\mu$m' %(epsilon_b,d)

if plot_xy_3D ==1 :
    omegaTHz0 = 0.01
    omegac0 = omegaTHz0*aux2  
    z0 = 1
    title1 = r'$\epsilon_1$ = %i, v = c/%i $\mu$m/s, $\omega$ = %.2fTHz, b = %.2f$\mu$m' %(epsi1,int_v,omegaTHz0,b) 
    label1 = '_b%i_3D' %(b)
    labelx,labely = 'x [$\mu$m]', 'y [$\mu$m]'
    title2A = r'$\epsilon_2$ = %i, $\hbar\mu$ = %.2feV, $\hbar\gamma$ = %.4feV, $z_p=%.2f\mu$m, z = %i$\mu$m' %(epsi2,Ebulk,hbgama,zp,z0) 
    title = title1 + '\n'  + title2A   + ', ' + title4

elif plot_vs_omegaTHz == 1:
    x0 = 0
    y0 = 0
    z0 = 0
    list_OmegaTHz = np.linspace(0.01,2.01,251)
    list_OmegaTHz = np.linspace(0.001,0.4001,51)
    list_E = np.array(list_OmegaTHz)*1e12*hb
    
    label1 = '_b%i_vsOmega' %(b)
    title1 = r'$\epsilon_1$ = %i, v = c/%i $\mu$m/s, x = %i$\mu$m, y = %i$\mu$m, b = %.2f$\mu$m' %(epsi1,int_v,x0,y0,b) 
    labelx = r'E [eV]'
    title2A = r'$\epsilon_2$ = %i, $\hbar\mu$ = %.2feV, $\hbar\gamma$ = %.4feV, $z_p=%.2f\mu$m, z = %i$\mu$m' %(epsi2,Ebulk,hbgama,zp,z0) 
    title = title1 + '\n'  + title2A   + ', ' + title4

elif plot_vs_z == 1:
    x0 = 10
    y0 = 10
    omegaTHz0 = 0.01
    omegac0 = omegaTHz0*aux2  
    
    list_z = np.linspace(b-2,zp,65)
    label1 = '_b%i_vs_z' %(b)
    title1 = r'$\epsilon_1$ = %i, v = c/%i $\mu$m/s, x = %i$\mu$m, y = %i$\mu$m, b = %i$\mu$m' %(epsi1,int_v,x0,y0,b) 
    labelx = r'z [$\mu$m]' 
    
    title2B = r'$\omega$ = %.2fTHz, $\epsilon_2$ = %i, $\hbar\mu$ = %.2feV, $\hbar\gamma$ = %.4feV' %(omegaTHz0,epsi2,Ebulk,hbgama) 
    title3  = r'$z_p=%.2f\mu$m' %(zp)
    
    title = title1 + '\n'  + title2B + '\n' + title3  + ', ' + title4


#%%

tamfig = (4.5,3.5)
tamlegend = 10
tamletra = 11
tamtitle = 10
tamnum = 9
labelpady = -1.5
labelpadx = 0.8
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

color_rgb1 = (1.0, 0.01, 0.24)
color_rgb2 = (0.8, 0.58, 0.46)
color_rgb3 = (0.55, 0.71, 0)
color_rgb4 = (0.6, 0.73, 0.89)

#%%

if plot_xy_3D ==1 : 
    N = 100
    list_x = np.linspace(-cota,cota,N)
    list_y = np.linspace(-cota,cota,N)
    X, Y = np.meshgrid(list_x, list_y)
    limits = [min(list_x) , max(list_x), min(list_y) , max(list_y)]

    def functionQE_real(x,y):
        rta = Efield_NUM_QE(omegac0,epsi1,epsi2,omega_bulk,hbgama,d,epsilon_b,int_v,b,zp,x,y,z0)
        return rta.real
    
    def function_real(x,y):
        rta = Efield_NUM(omegac0,epsi1,epsi2,omega_bulk,hbgama,d,epsilon_b,int_v,b,zp,x,y,z0)
        return rta.real
    
    f1_real = np.vectorize(functionQE_real)
    Z1_real = f1_real(X, Y)
    
    f2_real = np.vectorize(function_real)
    Z2_real = f2_real(X, Y)
    
    def functionQE_imag(x,y):
        rta = Efield_NUM_QE(omegac0,epsi1,epsi2,omega_bulk,hbgama,d,epsilon_b,int_v,b,zp,x,y,z0)
        return rta.imag
    
    def function_imag(x,y):
        rta = Efield_NUM(omegac0,epsi1,epsi2,omega_bulk,hbgama,d,epsilon_b,int_v,b,zp,x,y,z0)
        return rta.imag
    
    f1_imag = np.vectorize(functionQE_imag)
    Z1_imag = f1_imag(X, Y)
    
    f2_imag = np.vectorize(function_imag)
    Z2_imag = f2_imag(X, Y)
    
    os.chdir(path_save)
    
    np.savetxt('X_cota%i' %(cota) + label1 + '.txt', X, fmt='%1.11e', delimiter='\t', header = title)
    np.savetxt('Y_cota%i' %(cota) + label1 + '.txt', X, fmt='%1.11e', delimiter='\t', header = title)
    
    if plot_sin_QE == 0:
        header1 = 'Re E_{ref,x} QE num' + ', ' +  title
        # tabla1 = np.array(Z1)
        np.savetxt('Re_Eref_QE_cota%i' %(cota) + label1 + '.txt', Z1_real, fmt='%1.11e', delimiter='\t', header = header1)
    else:
        header1 = 'Re E_{ref,x} QE num' + ', ' +  title
        # tabla1 = np.array(Z1)
        np.savetxt('Re_Eref_QE_cota%i' %(cota) + label1 + '.txt', Z1_real, fmt='%1.11e', delimiter='\t', header = header1)
        
        header2 = 'Re E_{ref,x} num' + ', ' +  title
        # tabla2 = np.array([Z2])    
        np.savetxt('Re_Eref_cota%i' %(cota) + label1 + '.txt', Z2_real, fmt='%1.11e', delimiter='\t', header = header2)    
     
    #####################################
    
    if plot_sin_QE == 0:
        header1 = 'Im E_{ref,x} QE num' + ', ' +  title
        # tabla1 = np.array(Z1)
        np.savetxt('Im_Eref_QE_cota%i' %(cota) + label1 + '.txt', Z1_imag, fmt='%1.11e', delimiter='\t', header = header1)
    else:
        header1 = 'Im E_{ref,x} QE num' + ', ' +  title
        # tabla1 = np.array(Z1)
        np.savetxt('Im_Eref_QE_cota%i' %(cota) + label1 + '.txt', Z1_imag, fmt='%1.11e', delimiter='\t', header = header1)
        
        header2 = 'Im E_{ref,x} num' + ', ' +  title
        # tabla2 = np.array([Z2])    
        np.savetxt('Im_Eref_cota%i' %(cota) + label1 + '.txt', Z2_imag, fmt='%1.11e', delimiter='\t', header = header2)   
    
       
    print('Graficar el External field' + label1)
    
    graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    im = plt.imshow(Z1_real, extent = limits, cmap=plt.cm.hot, interpolation='bilinear')
    cbar = plt.colorbar(im)
    cbar.ax.tick_params(labelsize = tamnum)
    if save_graphs == True:
        plt.savefig('Re_Eref_QE_cota%i' %(cota) + label1 + '.png', format='png')
    
    
    graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    im = plt.imshow(Z2_real, extent = limits, cmap=plt.cm.hot, interpolation='bilinear')
    cbar = plt.colorbar(im)
    cbar.ax.tick_params(labelsize = tamnum)
    if save_graphs == True:
        plt.savefig('Re_Eref_cota%i' %(cota) + label1 + '.png', format='png')
        
    
    graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    im = plt.imshow(Z1_imag, extent = limits, cmap=plt.cm.hot, interpolation='bilinear')
    cbar = plt.colorbar(im)
    cbar.ax.tick_params(labelsize = tamnum)
    if save_graphs == True:
        plt.savefig('Imag_Eref_QE_cota%i' %(cota) + label1 + '.png', format='png')
    
    
    graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    im = plt.imshow(Z2_imag, extent = limits, cmap=plt.cm.hot, interpolation='bilinear')
    cbar = plt.colorbar(im)
    cbar.ax.tick_params(labelsize = tamnum)
    if save_graphs == True:
        plt.savefig('Imag_Eref_cota%i' %(cota) + label1 + '.png', format='png')    

#%%

elif plot_vs_omegaTHz == 1:
    
    listyQE_re = []
    listy_re = []
    listyQE_im = []
    listy_im = []
    
    listySP_re = []
    listySP_im = []
    
    listy_ana1_re = []
    listy_ana1_im = []
    
    listy_ana1b_re = []
    listy_ana1b_im = []
    
    
    listy_ana2_re = []
    listy_ana2_im = []
    
    
    j = 0
    for omegaTHz in list_OmegaTHz:
        omegaTHz = np.round(omegaTHz,8)
        omegac = omegaTHz*aux2  

        rta = Efield_NUM_QE(omegac,epsi1,epsi2,omega_bulk,hbgama,d,epsilon_b,int_v,b,zp,x0,y0,z0)
        listyQE_re.append(rta.real)
        listyQE_im.append(rta.imag)
        
        rta = Efield_NUM_SP(omegac,epsi1,epsi2,omega_bulk,hbgama,d,epsilon_b,int_v,b,zp,z0)
        listySP_re.append(rta.real)
        listySP_im.append(rta.imag)
    
        rta = Efield_ana1(omegac,epsi1,epsi2,omega_bulk,hbgama,d,epsilon_b,int_v,b,zp,z0)
        listy_ana1_re.append(rta.real)
        listy_ana1_im.append(rta.imag)
    
        rta = Efield_ana2(omegac,epsi1,epsi2,omega_bulk,hbgama,d,epsilon_b,int_v,b,zp,z0)
        listy_ana2_re.append(rta.real)
        listy_ana2_im.append(rta.imag)

        rta = Efield_ana1b(omegac,epsi1,epsi2,omega_bulk,hbgama,d,epsilon_b,int_v,b,zp,z0)
        listy_ana1b_re.append(rta.real)
        listy_ana1b_im.append(rta.imag)
        
        
        if plot_sin_QE == 1:
            rta = Efield_NUM(omegac,epsi1,epsi2,omega_bulk,hbgama,d,epsilon_b,int_v,b,zp,x0,y0,z0)
            listy_re.append(rta.real)
            listy_im.append(rta.imag)
        
        j = j + 1
        print(j)
    
    graph(title,labelx,r'Re$(E_{ref,x})$',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    plt.plot(list_E,listyQE_re,'.',ms = ms,color = 'lightseagreen',label = 'numerical QE')
    if plot_sin_QE == 1:
        plt.plot(list_E,listy_re,'--',lw = 2,color = 'darkviolet',label = 'numerical')
    # plt.yscale('log')
    plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
    if save_graphs==1:
        plt.tight_layout()
        os.chdir(path_save)
        plt.savefig( 'Re_Eref1_' + label1 + '.png', format='png')

    graph(title,labelx,r'Re$(E_{ref,x})$',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    plt.plot(list_E,listySP_re,'--',ms = ms,color = 'red',label = 'numerical SP')
    plt.plot(list_E,listy_ana1_re,'.',ms = ms,color = 'blue',label = 'analytical 1')
    plt.plot(list_E,listy_ana1b_re,'.',ms = ms,color = 'magenta',label = 'analytical 1b')
    plt.plot(list_E,listy_ana2_re,'--',ms = ms,color = 'green',label = 'analytical 2')
    plt.yscale('log')
    plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
    if save_graphs==1:
        plt.tight_layout()
        os.chdir(path_save)
        plt.savefig( 'Re_Eref2_' + label1 + '.png', format='png')
    
    graph(title,labelx,r'Im$(E_{ref,x})$',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    plt.plot(list_E,listyQE_im,'.',ms = ms,color = 'lightseagreen',label = 'numerical QE')
    if plot_sin_QE == 1:
        plt.plot(list_E,listy_im,'--',lw = 2,color = 'darkviolet',label = 'numerical')
    # plt.yscale('log')
    plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
    if save_graphs==1:
        plt.tight_layout()
        plt.savefig( 'Im_Eref1_' + label1 + '.png', format='png')    

    graph(title,labelx,r'|Im$(E_{ref,x})$|',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    plt.plot(list_E,np.abs(listySP_im),'--',ms = ms,color = 'red',label = 'numerical SP')
    plt.plot(list_E,np.abs(listy_ana1_im),'.',ms = ms,color = 'blue',label = 'analytical 1')
    plt.plot(list_E,np.abs(listy_ana1b_im),'.',ms = ms,color = 'magenta',label = 'analytical 1b')
    plt.plot(list_E,np.abs(listy_ana2_im),'--',ms = ms,color = 'green',label = 'analytical 2')
    plt.yscale('log')
    plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
    if save_graphs==1:
        plt.tight_layout()
        plt.savefig( 'Im_Eref2_' + label1 + '.png', format='png')      
                    
#%%
    

elif plot_vs_z == 1:
    
    listyQE_re = []
    listy_re = []
    listyQE_im = []
    listy_im = []
    
    j = 0
    for z in list_z:
        
        z = np.round(z,8)
        rta = Efield_NUM_QE(omegac0,epsi1,epsi2,omega_bulk,hbgama,d,epsilon_b,int_v,b,zp,x0,y0,z)
        listyQE_re.append(rta.real)
        listyQE_im.append(rta.imag)
        
        if plot_sin_QE == 1:
            rta = Efield_NUM(omegac0,epsi1,epsi2,omega_bulk,hbgama,d,epsilon_b,int_v,b,zp,x0,y0,z)
            listy_re.append(rta.real)
            listy_im.append(rta.imag)
        
        j = j + 1
        print(j)
    
    graph(title,labelx,r'Re$(E_{ref,x})$',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    plt.plot(list_z,listyQE_re,'.',ms = ms,color = 'lightseagreen',label = 'numerical QE')
    if plot_sin_QE == 1:
        plt.plot(list_z,listy_re,'--',lw = 2,color = 'darkviolet',label = 'numerical')
    # plt.yscale('log')
    plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
    if save_graphs==1:
        plt.tight_layout()
        os.chdir(path_save)
        plt.savefig( 'Re_Eref_' + label1 + '.png', format='png')
    
    graph(title,labelx,r'Im$(E_{ref,x})$',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    plt.plot(list_z,listyQE_im,'.',ms = ms,color = 'lightseagreen',label = 'numerical QE')
    if plot_sin_QE == 1:
        plt.plot(list_z,listy_im,'--',lw = 2,color = 'darkviolet',label = 'numerical')
    # plt.yscale('log')
    plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
    if save_graphs==1:
        plt.tight_layout()
        plt.savefig( 'Im_Eref_' + label1 + '.png', format='png')    

    np.savetxt('list_z' + label1 + '.txt', list_z, fmt='%1.11e', delimiter='\t', header = title)
    if plot_sin_QE == 0:
        header1 = 'E_{ref,x} QE num' + ', ' +  title
        np.savetxt('Re_Eref_QE'  + label1 + '.txt', listyQE_re, fmt='%1.11e', delimiter='\t', header = header1)
        np.savetxt('Im_Eref_QE'  + label1 + '.txt', listyQE_im, fmt='%1.11e', delimiter='\t', header = header1)
        
    else:
        
        header1 = 'E_{ref,x} QE num' + ', ' +  title
        np.savetxt('Re_Eref_QE'  + label1 + '.txt', listyQE_re, fmt='%1.11e', delimiter='\t', header = header1)
        np.savetxt('Im_Eref_QE'  + label1 + '.txt', listyQE_im, fmt='%1.11e', delimiter='\t', header = header1)
        
        header2 = 'E_{ref,x} num' + ', ' +  title
        np.savetxt('Re_Eref'  + label1 + '.txt', listy_re, fmt='%1.11e', delimiter='\t', header = header2)  
        np.savetxt('Im_Eref'  + label1 + '.txt', listy_im, fmt='%1.11e', delimiter='\t', header = header2)         
                
#%%