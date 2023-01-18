#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: leila

graficar green tensor Gnn con la misma convencion 
que el paper 370 (z hacia abajo, zplane>0, particula en r = 0, medio 1 arriba del plano, 
                  medio 2 abajo del plano)
"""

import numpy as np
import sys
import os 
import matplotlib.pyplot as plt
import seaborn as sns

#%%

list_type_plane = ['graphene','Ag'] ### revisar conductividad de Ag antes de hacer el green tensor reflected
type_plane = list_type_plane[0]

save_graphs = 1 #guardar los graficos 2D del campo


plot_vs_z = 0 #graficar vs z
plot_vs_omegaTHz = 1 #graficar vs omega en THz
plot_vs_xD = 0
sns.set()

if plot_vs_z + plot_vs_omegaTHz + plot_vs_xD > 1:
    raise TypeError('Choose one option')

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('GreenTensor/' + 'self','')

if type_plane == 'graphene':
    err = 'green_tensor_ref_xx_graphene.py no se encuentra en ' + path_basic
    try:
        sys.path.insert(1, path_basic)
        from green_tensor_self_xx_graphene import green_tensor_NUM_QE, green_tensor_NUM, green_tensor_NUM_PP_QE, green_tensor_ana1, green_tensor_ana2
    except ModuleNotFoundError:
        print(err)
    path_save = path_basic + '/' + 'green_tensor_self_xx_graphene'  


else:
    err = 'green_tensor_self_xx_Ag.py no se encuentra en ' + path_basic
    try:
        sys.path.insert(1, path_basic)
        from green_tensor_self_xx_Ag import green_tensor_NUM_QE, green_tensor_NUM
    except ModuleNotFoundError:
        print(err)
    path_save = path_basic + '/' + 'green_tensor_self_xx_Ag'

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

zD = 0
xD, yD = 0,0
x,y = 0,0
#omega = 0.7*1e12
epsi1 = 1
epsi2 = 1
zp = 0.01 #micrones # posicion del plano (z apunta hacia abajo)
title1 = r'$\epsilon_1$ = %i, $\epsilon_2$ = %i, x = %i$\mu$m, y = %i$\mu$m' %(epsi1,epsi2,x,y)

if type_plane == 'graphene':
    hbmu, hbgama = 0.4, 0.0001 #potencial quimico in ev, freq de collision in ev
    title3 = r'$\hbar \mu$ = %.2feV, $\hbar \gamma$ = %.4feV, zp = %.2f$\mu$m' %(hbmu,hbgama,zp)
    
    def functionQE(omegac,z,xD):
        return green_tensor_NUM_QE(omegac,epsi1,epsi2,hbmu,hbgama,x,y,z,xD,yD,zD,zp)
    
    def function(omegac,z,xD):    
        return green_tensor_NUM(omegac,epsi1,epsi2,hbmu,hbgama,x,y,z,xD,yD,zD,zp)

    def functionSP(omegac,z,xD):
        return green_tensor_NUM_PP_QE(omegac,epsi1,epsi2,hbmu,hbgama,z,xD,zD,zp)
    
    def functionANA1(omegac,z,xD):
        return green_tensor_ana1(omegac,epsi1,epsi2,hbmu,hbgama,z,xD,zD,zp)
    
    def functionANA2(omegac,z,xD):
        return green_tensor_ana2(omegac,epsi1,epsi2,hbmu,hbgama,z,xD,zD,zp)
    
    
    if save_graphs==1:    
        path_save = path_save + '/' + 'mu_%.2f' %(hbmu) 
        if not os.path.exists(path_save):
            print('Creating folder to save graphs')
            os.mkdir(path_save)
        
else:
    energy_bulk = 9.17
    omega_bulk, hbar_gamma_in, d = energy_bulk/hb, 0.0001, 0.09 
    title3 = r'$E_{bulk}$ = %.2feV, $\hbar \gamma$ = %.4feV, d = %.2f$\mu$m, zp = %.2f$\mu$m' %(energy_bulk,hbar_gamma_in,d,zp)

    def functionQE(omegac,z,xD):
        return green_tensor_NUM_QE(omegac,epsi1,epsi2,omega_bulk,hbar_gamma_in,d,x,y,z,xD,yD,zD,zp)
    
    def function(omegac,z,xD):    
        return green_tensor_NUM(omegac,epsi1,epsi2,omega_bulk,hbar_gamma_in,d,x,y,z,xD,yD,zD,zp)    
    
########################
if plot_vs_z == 1:
    omegaTHz = 0.7
    omegaTHz0 = 1.511
    omegaTHz0 = 0.01
    # omegaTHz = 1.663 - 0.152j
    
    omegac0 = omegaTHz*1e12/c 
    xD0 = 0    

    list_z = np.linspace(1,31,161)
    list_x = list_z
    labelx = '$z$ [$\mu$m]'
    label1 = '_vs_z'
    title2 = r'$\omega$ = %.2f THz, $x_D$ = %i$\mu$m' %(omegaTHz0,xD0)
    
########################
if plot_vs_omegaTHz == 1:
    z0 = 0
    xD0 = 0    
    list_OmegaTHz = np.linspace(0.01,2.01,251)
    list_OmegaTHz = np.linspace(0.001,0.201,101)
    listE = np.array(list_OmegaTHz)*1e12*hb
    list_x = listE 
    labelx = 'E [eV]'
    label1 = '_vs_OmegaTHz'
    title2 = r'z = %.2f $\mu$m, $x_D$ = %i$\mu$m' %(z0,xD0)
    
if plot_vs_xD == 1:
    omegaTHz0  = 0.70
    omegac0 = omegaTHz0*aux2
    z0 = 8
    
    list_xD = np.linspace(-80,80,321)
    list_xD  = np.delete(list_xD , 160) #sacar el cero porque diverge
    list_x = list_xD
    labelx = r'$x_D$ [$\mu$m]'
    label1 = '_vs_xD'
    title2 = r'$\omega$ = %.2f THz, z = %.2f$\mu$m' %(omegaTHz0,z0)


#%%

title1 = r'$\epsilon_1$ = %i, $\epsilon_2$ = %i, x = %i$\mu$m, y = %i$\mu$m, $y_D$ = %i$\mu$m' %(epsi1,epsi2,x,y,yD)
title3 = r'$z_D$ = %i$\mu$m, $\hbar \mu$ = %.2feV, $\hbar \gamma$ = %.4feV, zp = %.2f$\mu$m' %(zD,hbmu,hbgama,zp) 
title = title1 + '\n' + title2 + '\n ' + title3

tamfig = (4.5,3.5)
tamlegend = 10
tamletra = 11
tamtitle = 10
tamnum = 9
labelpady = -1.5
labelpadx = 0
pad = 0.5
mk = 2
ms = 4
hp = 0.3
hspace = 0
wspace = 0.01
loc1 = [0.3,0.935] 

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

listI05_re = []
listI25_re = []
listI06_re = []
listI26_re = []
total1_re = []

listI05_im = []
listI25_im = []
listI06_im = []
listI26_im = []
total1_im = []

listI05_v3_re = []
listI25_v3_re = []
listI06_v3_re = []
listI26_v3_re = []
total2_re = []

listI05_v3_im = []
listI25_v3_im = []
listI06_v3_im = []
listI26_v3_im = []
total2_im = []


sp_re = []
sp_im = []

ana1_re = []
ana1_im = []

ana2_re = []
ana2_im = []

if plot_vs_z == 1:
    
    # if omegaTHz.imag == 0:
    #     title2A = r'$\omega$ = %.2fTHz, $x_D$ = %i$\mu$m, $y_D$ = %i$\mu$m, $z_D$ = %i$\mu$m' %(omegaTHz.real,xD,yD,zD)
    # elif omegaTHz.imag <0:
    #     title2A = r'$\omega$ = (%.2f %.2f)THz, $x_D$ = %i$\mu$m, $y_D$ = %i$\mu$m, $z_D$ = %i$\mu$m' %(omegaTHz.real,omegaTHz.imag,xD,yD,zD)
    # else:
    #     title2A = r'$\omega$ = (%.2f + %.2f)THz, $x_D$ = %i$\mu$m, $y_D$ = %i$\mu$m, $z_D$ = %i$\mu$m' %(omegaTHz.real,omegaTHz.imag,xD,yD,zD)
 
    # title = title1 + '\n' + title2A + '\n' + title3


    for z in list_z: 
        
        z = np.round(z,8)
        
        I0_5, I2_5, I0_6, I2_6 = functionQE(omegac0,z,xD0)
        I0_5_v3, I2_5_v3, I0_6_v3, I2_6_v3 = function(omegac0,z,xD0)
        
        rta_SP = functionSP(omegac0,z,xD0)
        rta_ana1 = functionANA1(omegac0,z,xD0)
        rta_ana2 = functionANA2(omegac0,z,xD0)

        sp_re.append(rta_SP.real)
        sp_im.append(rta_SP.imag)
        
        ana1_re.append(rta_ana1.real)
        ana1_im.append(rta_ana1.imag)

        ana2_re.append(rta_ana2.real)
        ana2_im.append(rta_ana2.imag)

        
        listI05_re.append(I0_5.real)
        listI25_re.append(I2_5.real)
        listI06_re.append(I0_6.real)
        listI26_re.append(I2_6.real)
        
        listI05_im.append(I0_5.imag)
        listI25_im.append(I2_5.imag)
        listI06_im.append(I0_6.imag)
        listI26_im.append(I2_6.imag)
        
        tot1 = I0_5 + I2_5 + I0_6 + I2_6
        total1_re.append(tot1.real)
        total1_im.append(tot1.imag)
        
        listI05_v3_re.append(I0_5_v3.real)
        listI25_v3_re.append(I2_5_v3.real)
        listI06_v3_re.append(I0_6_v3.real)
        listI26_v3_re.append(I2_6_v3.real)
        
        listI05_v3_im.append(I0_5_v3.imag)
        listI25_v3_im.append(I2_5_v3.imag)
        listI06_v3_im.append(I0_6_v3.imag)
        listI26_v3_im.append(I2_6_v3.imag)

        tot2 = I0_5_v3 + I2_5_v3 + I0_6_v3 + I2_6_v3
        total2_re.append(tot2.real)
        total2_im.append(tot2.imag)
    del z
elif plot_vs_omegaTHz == 1 :


    for OmegaTHz in list_OmegaTHz: 
        
        OmegaTHz = np.round(OmegaTHz,8)
        Omegac = OmegaTHz*1e12/c
        
        I0_5, I2_5, I0_6, I2_6 = functionQE(Omegac,z0,xD0)
        I0_5_v3, I2_5_v3, I0_6_v3, I2_6_v3 = function(Omegac,z0,xD0)      

        rta_SP = functionSP(Omegac,z0,xD0)
        rta_ana1 = functionANA1(Omegac,z0,xD0)
        rta_ana2 = functionANA2(Omegac,z0,xD0)

        sp_re.append(rta_SP.real)
        sp_im.append(rta_SP.imag)
        
        ana1_re.append(rta_ana1.real)
        ana1_im.append(rta_ana1.imag)

        ana2_re.append(rta_ana2.real)
        ana2_im.append(rta_ana2.imag)
        
        listI05_re.append(I0_5.real)
        listI25_re.append(I2_5.real)
        listI06_re.append(I0_6.real)
        listI26_re.append(I2_6.real)
        
        listI05_im.append(I0_5.imag)
        listI25_im.append(I2_5.imag)
        listI06_im.append(I0_6.imag)
        listI26_im.append(I2_6.imag)

        tot1 = I0_5 + I2_5 + I0_6 + I2_6
        total1_re.append(tot1.real)
        total1_im.append(tot1.imag)
        
        listI05_v3_re.append(I0_5_v3.real)
        listI25_v3_re.append(I2_5_v3.real)
        listI06_v3_re.append(I0_6_v3.real)
        listI26_v3_re.append(I2_6_v3.real)
        
        listI05_v3_im.append(I0_5_v3.imag)
        listI25_v3_im.append(I2_5_v3.imag)
        listI06_v3_im.append(I0_6_v3.imag)
        listI26_v3_im.append(I2_6_v3.imag)

        tot2 = I0_5_v3 + I2_5_v3 + I0_6_v3 + I2_6_v3
        total2_re.append(tot2.real)
        total2_im.append(tot2.imag)

    del Omegac
elif plot_vs_xD == 1 :


    for xD in list_xD: 
        
        I0_5, I2_5, I0_6, I2_6 = functionQE(omegac0,z0,xD)
        I0_5_v3, I2_5_v3, I0_6_v3, I2_6_v3 = function(omegac0,z0,xD)


        rta_SP = functionSP(omegac0,z0,xD)
        rta_ana1 = functionANA1(omegac0,z0,xD)
        rta_ana2 = functionANA2(omegac0,z0,xD)

        sp_re.append(rta_SP.real)
        sp_im.append(rta_SP.imag)
        
        ana1_re.append(rta_ana1.real)
        ana1_im.append(rta_ana1.imag)

        ana2_re.append(rta_ana2.real)
        ana2_im.append(rta_ana2.imag)
        
        listI05_re.append(I0_5.real)
        listI25_re.append(I2_5.real)
        listI06_re.append(I0_6.real)
        listI26_re.append(I2_6.real)
        
        listI05_im.append(I0_5.imag)
        listI25_im.append(I2_5.imag)
        listI06_im.append(I0_6.imag)
        listI26_im.append(I2_6.imag)

        tot1 = I0_5 + I2_5 + I0_6 + I2_6
        total1_re.append(tot1.real)
        total1_im.append(tot1.imag)
        
        listI05_v3_re.append(I0_5_v3.real)
        listI25_v3_re.append(I2_5_v3.real)
        listI06_v3_re.append(I0_6_v3.real)
        listI26_v3_re.append(I2_6_v3.real)
        
        listI05_v3_im.append(I0_5_v3.imag)
        listI25_v3_im.append(I2_5_v3.imag)
        listI06_v3_im.append(I0_6_v3.imag)
        listI26_v3_im.append(I2_6_v3.imag)

        tot2 = I0_5_v3 + I2_5_v3 + I0_6_v3 + I2_6_v3
        total2_re.append(tot2.real)
        total2_im.append(tot2.imag)


#%%

print('Graficar el green tensor ' + label1)
    
# graph(title,labelx,r'Re $I^{(0)}_7$',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
# plt.plot(list_x,listI05_re,'.',ms = ms,color = color_rgb1,label = 'numerical QE')
# # plt.yscale('log')
# plt.plot(list_x,listI05_v3_re,'-',color = 'green',label = 'numerical')
# plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
# if save_graphs==1:
#     plt.tight_layout()
#     os.chdir(path_save)
#     plt.savefig( 'green_tensor_Re_I07' + label1 + '.png', format='png')

# graph(title,labelx,r'Im $I^{(0)}_7$',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
# plt.plot(list_x,listI05_im,'.',ms = ms,color = color_rgb1,label = 'numerical QE')
# # plt.yscale('log')
# plt.plot(list_x,listI05_v3_im,'-',color = 'green',label = 'numerical')
# plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
# if save_graphs==1:
#     plt.tight_layout()
#     plt.savefig( 'green_tensor_Im_I07' + label1 + '.png', format='png')

# #############################################################################################
# graph(title,labelx,r'Re $I^{(2)}_7$',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
# plt.plot(list_x,listI25_re,'.',ms = ms,color = color_rgb3,label = 'numerical QE')
# # plt.yscale('log')
# plt.plot(list_x,listI25_v3_re,'-',color = 'green',label = 'numerical')
# plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
# if save_graphs==1:
#     plt.tight_layout()
#     plt.savefig( 'green_tensor_Re_I27' + label1 + '.png', format='png')

# graph(title,labelx,r'Im $I^{(2)}_7$',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
# plt.plot(list_x,listI25_im,'.',ms = ms,color = color_rgb3,label = 'numerical QE')
# # plt.yscale('log')
# plt.plot(list_x,listI25_v3_im,'-',color = 'green',label = 'numerical')
# plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
# if save_graphs==1:
#     plt.tight_layout()
#     plt.savefig( 'green_tensor_Im_I27' + label1 + '.png', format='png')
# #############################################################################################
# # if graphs_vs_z == 1:
# graph(title,labelx,r'Re $I^{(0)}_8$',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
# plt.plot(list_x,listI06_re,'.',ms = ms,color = color_rgb2,label = 'numerical QE')
# plt.plot(list_x,listI06_v3_re,'-',color = 'green',label = 'numerical')
# plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
# if save_graphs==1:
#     plt.tight_layout()
#     plt.savefig( 'green_tensor_Re_I08'+ label1 + '.png', format='png')
# # elif graphs_vs_omegaTHz == 1 and x==0 and y==0:
# #     labely = r'Re $I^{(0)}_8$'
# #     fig,axs = plt.subplots(2,1, sharex=True, facecolor='w', figsize = tamfig)
# #     plt.subplots_adjust(hspace =hspace,wspace = wspace)
    
# #     axs[0].plot(list_x,listI06_re,'.',ms = ms,color = color_rgb2,label = 'numerical QE')
# #     axs[1].plot(list_x,listI06_v3_re,'-',color = 'green',label =  'numerical')
    
# #     # axs[0].set_ylabel(labely,fontsize = tamletra,labelpad =labelpady)
# #     # axs[1].set_ylabel(labely,fontsize = tamletra,labelpad =labelpady)
# #     axs[1].set_xlabel(labelx,fontsize = tamletra,labelpad =0.75) 
# #     for i in [0,1]:
# #         axs[i].minorticks_on()        
# #         axs[i].tick_params(labelsize = tamnum,pad = pad)
# #     fig.legend(loc = loc1,ncol = 2,markerscale=mk,fontsize=tamlegend, 
# #                       columnspacing = 2,frameon=0,handletextpad=0.2, handlelength=1)  

# #     if save_graphs==1:
# #         plt.tight_layout()
# #         plt.savefig( 'green_tensor_Re_I08' + label1 + '.png', format='png') 
 
# graph(title,labelx,r'Im $I^{(0)}_8$',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
# plt.plot(list_x,listI06_im,'.',ms = ms,color = color_rgb2,label = 'numerical QE')
# plt.plot(list_x,listI06_v3_im,'-',color = 'green',label = 'numerical')
# plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
# if save_graphs==1:
#     plt.tight_layout()
#     plt.savefig( 'green_tensor_Im_I08'+ label1 + '.png', format='png')
# ############################################################################################
# graph(title,labelx,r'Re $I^{(2)}_8$',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
# plt.plot(list_x,listI26_re,'.',ms = ms,color = color_rgb4,label = 'numerical QE')
# # plt.yscale('log')
# plt.plot(list_x,listI26_v3_re,'-',color = 'green',label = 'numerical')
# plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
# if save_graphs==1:
#     plt.tight_layout()
#     plt.savefig( 'green_tensor_Re_I28'+ label1 + '.png', format='png')

# graph(title,labelx,r'Im $I^{(2)}_8$',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
# plt.plot(list_x,listI26_im,'.',ms = ms,color = color_rgb4,label = 'numerical QE')
# # plt.yscale('log')
# plt.plot(list_x,listI26_v3_im,'-',color = 'green',label = 'numerical')
# plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
# if save_graphs==1:
#     plt.tight_layout()
#     plt.savefig( 'green_tensor_Im_I28'+ label1 + '.png', format='png')

############################################################################################

# if graphs_vs_z == 1:
graph(title,labelx,r'Re(Reflected Green Tensor)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(list_x,total1_re,'.',ms = ms,color = 'black',label = 'numerical QE')
plt.plot(list_x,total2_re,'-',color = 'green',label =  'numerical')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
# plt.yscale('log')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
if save_graphs==1:
    plt.tight_layout()
    plt.savefig( 'green_tensor_ref_Re1' + label1 + '.png', format='png')


graph(title,labelx,r'Im(Reflected Green Tensor)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(list_x,total1_im,'.',ms = ms,color = 'black',label = 'numerical QE')
plt.plot(list_x,total2_im,'-',color = 'green',label =  'numerical')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
# plt.yscale('log')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
if save_graphs==1:
    plt.tight_layout()
    plt.savefig( 'green_tensor_ref_Im1' + label1 + '.png', format='png')



graph(title,labelx,r'Re(Reflected Green Tensor)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(list_x,sp_re,'-',color = 'blue',label =  'numerical SP')
plt.plot(list_x,ana1_re,'.',color = 'green',label =  'ana 1')
plt.plot(list_x,ana2_re,'-',color = 'red',label =  'ana 2')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
# plt.yscale('log')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
if save_graphs==1:
    plt.tight_layout()
    plt.savefig( 'green_tensor_ref_Re1' + label1 + '.png', format='png')


graph(title,labelx,r'Im(Reflected Green Tensor)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(list_x,sp_im,'-',color = 'blue',label =  'numerical SP')
plt.plot(list_x,ana1_im,'.',color = 'green',label =  'ana 1')
plt.plot(list_x,ana2_im,'-',color = 'red',label =  'ana 2')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
# plt.yscale('log')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
if save_graphs==1:
    plt.tight_layout()
    plt.savefig( 'green_tensor_ref_Im1' + label1 + '.png', format='png')
    
#%%
