
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


plot_num = 0

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/potential_field/infinite_dipoles','')
path_save = path_basic + '/' + 'many_dipoles_phi_ind3D'
if not os.path.exists(path_save):
    print('Creating folder to save graphs')
    os.mkdir(path_save)

err = 'many_potential_integrals.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from potential_many_dipoles2 import phi_many_dipoles_ana_ind
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
cota = 75 #nanometros
epsi1,epsi2 = 1,1
hbmu,hbgama = 0.3,0.0001
zp = 0.05
b = -0.01

x0,z0 = 0, zp 

#Ndipoles = 50
#list_xD = np.linspace(-0.01,0.01,Ndipoles)

Nmax = 3

px = 1
py = 0
pz = 0

int_v = 20
int_v = 10

E0 = 30 # meV 
E0 = 43 # meV 
omegac0 = E0*1e-3/aux

def k_SP(E_mev,hbmu,hbgama):
    num = 2/alfac 
    omegac = E_mev*1e-3/aux
    cte = E_mev*1e-3 + 1j*hbgama/(4*hbmu)
    
    return omegac*num*cte

lambda_SP = 2*np.pi/k_SP(E0,hbmu,hbgama)
print(k_SP(E0,hbmu,hbgama))
a_min = np.real(lambda_SP)*Nmax/(int_v - 1)
a_max = np.real(lambda_SP)*Nmax/(int_v + 1)

lambda_SP_nm = np.real(lambda_SP)*1e3
cota_x = int(lambda_SP_nm*2)

a = 6.5
a = 0.05
a = np.mean([a_min,a_max])

theta0 = np.arccos(int_v-np.real(lambda_SP)*Nmax/a)

title2 = r'v = c/%i, b = %inm, $z_p$=%inm, z=%inm' %(int_v,b*1e3,zp*1e3,z0*1e3) 
title3 = r'px=%i, py=%i, pz=%i, a = %inm, x=%inm' %(px,py,pz,a*1e3,x0*1e3)
title4 = r'$\theta$ = %.2f, N = %i' %(theta0,Nmax)

labelp = r'_px%i_py%i_pz%i_a%inm_N%i' %(px,py,pz,a*1e3,Nmax)

labelx = 'x [nm]'  
labely = 'y [nm]'
title = title2 + '\n' +  title3 + '\n' + title4 + ', E = %.2f meV' %(E0)
label1 = 'vs_y_E0%i' %(E0) + labelp
listy = np.linspace(50,5000,100)
listx = np.linspace(-2500,2500,100)
listy =listx

listy = np.linspace(-cota_x,cota_x,100)
listx = listy

#%%
#
#def function_real_num(x_nano,y_nano):
#    y = y_nano*1e-3 
#    x = x_nano*1e-3 
#
#    rta = phi_many_dipoles_num(omegac0,epsi1,epsi2,hbmu,hbgama,x,y,z0,zp,int_v,px,py,pz,a,Nmax)
#    # print(rta)
#    return rta.real
#
#def function_imag_num(x_nano,y_nano):
#    y = y_nano*1e-3 
#    x = x_nano*1e-3 
#    
#    rta = phi_many_dipoles_num(omegac0,epsi1,epsi2,hbmu,hbgama,x,y,z0,zp,int_v,px,py,pz,a,Nmax)
#    return rta.imag

#%%
    
def function_real_ana(x_nano,y_nano):
    y = y_nano*1e-3 
    x = x_nano*1e-3 
    rta = phi_many_dipoles_ana_ind(omegac0,epsi1,epsi2,hbmu,hbgama,x,y,z0,zp,int_v,px,py,pz,a,Nmax)
    # print(rta)
    return rta.real

def function_imag_ana(x_nano,y_nano):
    y = y_nano*1e-3 
    x = x_nano*1e-3 
    
    rta = phi_many_dipoles_ana_ind(omegac0,epsi1,epsi2,hbmu,hbgama,x,y,z0,zp,int_v,px,py,pz,a,Nmax)
    return rta.imag

#%%
    
#def function_real_pole_aprox(x_nano,y_nano):
#    y = y_nano*1e-3 
#    x = x_nano*1e-3 
#
#    rta = phi_many_dipoles_pole_aprox(omegac0,epsi1,epsi2,hbmu,hbgama,x,y,z0,zp,int_v,px,py,pz,a,Nmax)
#    # print(rta)
#    return rta.real
#
#def function_imag_pole_aprox(x_nano,y_nano):
#    
#    y = y_nano*1e-3 
#    x = x_nano*1e-3 
#    
#    rta = phi_many_dipoles_pole_aprox(omegac0,epsi1,epsi2,hbmu,hbgama,x,y,z0,zp,int_v,px,py,pz,a,Nmax)
#    return rta.imag

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

#%% analitico
#        
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
im = plt.imshow(Z_real_ana, extent = limits, cmap=plt.cm.hot, interpolation='bilinear')
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
cbar.set_label(r'Re($\phi$)')
os.chdir(path_save)
plt.tight_layout()
plt.savefig('Re_potential_ana' + label1 + '.png', format='png')




graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
vmin, vmax = np.min(Z_imag_ana), np.max(Z_imag_ana)
im = plt.imshow(Z_imag_ana, extent = limits, cmap=plt.cm.hot, interpolation='bilinear')
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
cbar.set_label(r'Im($\phi$)')
os.chdir(path_save)
plt.tight_layout()
plt.savefig('Im_potential_ana' + label1 + '.png', format='png')
    

#del Z_real_ana, Z_imag_ana

#%% num
#
#if plot_num == 1:
#    f_real_num = np.vectorize(function_real_num)
#    Z_real_num = f_real_ana(X, Y)
#    
#    f_imag_num = np.vectorize(function_imag_num)
#    Z_imag_num = f_imag_ana(X, Y)
#       
#    print('Graficar el External field' + label1)
#    
#    # if plot_xz_3D == 1 and px == 1: 
#    graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#    vmin, vmax = np.min(Z_real_num), np.max(Z_real_num)
##    pcm = plt.pcolormesh(X, Y, Z_real_num,
##                          norm=colors.SymLogNorm(linthresh=0.03, linscale=0.03,
##                                               vmin=int(vmin), vmax=int(vmax)),cmap=plt.cm.hot)
##    maxlog=int(np.ceil( np.log10( np.abs(vmax) )))
##    minlog=int(np.ceil( np.log10( np.abs(vmin) )))
##    if vmin < 0 :
##          tick_locations = ( [-(10.0**x) for x in np.linspace(minlog,-1,8)] 
##                            + [0] 
##                            + [(10.0**x) for x in np.linspace(-1,maxlog,8)] )
##    else:
##          tick_locations = ( [(10.0**x) for x in np.linspace(minlog,maxlog,maxlog + np.abs(minlog) + 1) ])    
##    #    im = plt.imshow(Z1_real, extent = limits, cmap=plt.cm.hot, interpolation='bilinear')
##    cbar = plt.colorbar(pcm)
#    im = plt.imshow(Z_real_num, extent = limits, cmap=plt.cm.hot, interpolation='bilinear')
#    cbar = plt.colorbar(im)
#    cbar.ax.tick_params(labelsize = tamnum)
##    cbar.set_ticks(tick_locations)
#    cbar.set_label(r'Re($\phi$)')
#    os.chdir(path_save)
#    plt.tight_layout()
#    plt.savefig('Re_potential_num_log' + label1 + '.png', format='png')
#    
#    
#    
#    
#    graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#    vmin, vmax = np.min(Z_imag_num), np.max(Z_imag_num)
##    pcm = plt.pcolormesh(X, Y, Z_imag_num,
##                          norm=colors.SymLogNorm(linthresh=0.03, linscale=0.03,
##                                               vmin=int(vmin), vmax=int(vmax)),cmap=plt.cm.hot)
##    maxlog=int(np.ceil( np.log10( np.abs(vmax) )))
##    minlog=int(np.ceil( np.log10( np.abs(vmin) )))
##    if vmin < 0 :
##          tick_locations = ( [-(10.0**x) for x in np.linspace(minlog,-1,8)] 
##                            + [0] 
##                            + [(10.0**x) for x in np.linspace(-1,maxlog,8)] )
##    else:
##          tick_locations = ( [(10.0**x) for x in np.linspace(minlog,maxlog,maxlog + np.abs(minlog) + 1) ])    
#    #    im = plt.imshow(Z1_real, extent = limits, cmap=plt.cm.hot, interpolation='bilinear')
#    im = plt.imshow(Z_imag_num, extent = limits, cmap=plt.cm.hot, interpolation='bilinear')
#    cbar = plt.colorbar(im)
#
##    cbar = plt.colorbar(pcm)
#    cbar.ax.tick_params(labelsize = tamnum)
##    cbar.set_ticks(tick_locations)
#    cbar.set_label(r'Im($\phi$)')
#    os.chdir(path_save)
#    plt.tight_layout()
#    plt.savefig('Im_potential_num_log' + label1 + '.png', format='png')
#        
#    
#    del Z_real_num, Z_imag_num

#%%

#
#if plot_num == 1:    
#    f_real_pole_aprox = np.vectorize(function_real_pole_aprox)
#    Z_real_pole_aprox = f_real_ana(X, Y)
#    
#    f_imag_pole_aprox = np.vectorize(function_imag_pole_aprox)
#    Z_imag_pole_aprox = f_imag_pole_aprox(X, Y)
#       
#    print('Graficar el External field' + label1)
#    
#    # if plot_xz_3D == 1 and px == 1: 
#    graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#    vmin, vmax = np.min(Z_real_pole_aprox), np.max(Z_real_pole_aprox)
##    pcm = plt.pcolormesh(X, Y, Z_real_pole_aprox,
##                          norm=colors.SymLogNorm(linthresh=0.03, linscale=0.03,
##                                               vmin=int(vmin), vmax=int(vmax)),cmap=plt.cm.hot)
##    maxlog=int(np.ceil( np.log10( np.abs(vmax) )))
##    minlog=int(np.ceil( np.log10( np.abs(vmin) )))
##    if vmin < 0 :
##          tick_locations = ( [-(10.0**x) for x in np.linspace(minlog,-1,8)] 
##                            + [0] 
##                            + [(10.0**x) for x in np.linspace(-1,maxlog,8)] )
##    else:
##          tick_locations = ( [(10.0**x) for x in np.linspace(minlog,maxlog,maxlog + np.abs(minlog) + 1) ])    
#    #    im = plt.imshow(Z1_real, extent = limits, cmap=plt.cm.hot, interpolation='bilinear')
##    cbar = plt.colorbar(pcm)
#
#    im = plt.imshow(Z_real_pole_aprox, extent = limits, cmap=plt.cm.hot, interpolation='bilinear')
#    cbar = plt.colorbar(im)
#    cbar.ax.tick_params(labelsize = tamnum)
##    cbar.set_ticks(tick_locations)
#    cbar.set_label(r'Re($\phi$)')
#    os.chdir(path_save)
#    plt.tight_layout()
#    plt.savefig('Re_potential_pole_aprox_log' + label1 + '.png', format='png')
#    
#    
#    
#    
#    graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#    vmin, vmax = np.min(Z_imag_pole_aprox), np.max(Z_imag_pole_aprox)
##    pcm = plt.pcolormesh(X, Y, Z_imag_pole_aprox,
##                          norm=colors.SymLogNorm(linthresh=0.03, linscale=0.03,
##                                               vmin=int(vmin), vmax=int(vmax)),cmap=plt.cm.hot)
##    maxlog=int(np.ceil( np.log10( np.abs(vmax) )))
##    minlog=int(np.ceil( np.log10( np.abs(vmin) )))
##    if vmin < 0 :
##          tick_locations = ( [-(10.0**x) for x in np.linspace(minlog,-1,8)] 
##                            + [0] 
##                            + [(10.0**x) for x in np.linspace(-1,maxlog,8)] )
##    else:
##          tick_locations = ( [(10.0**x) for x in np.linspace(minlog,maxlog,maxlog + np.abs(minlog) + 1) ])    
#    #    im = plt.imshow(Z1_real, extent = limits, cmap=plt.cm.hot, interpolation='bilinear')
##    cbar = plt.colorbar(pcm)
#    im = plt.imshow(Z_imag_pole_aprox, extent = limits, cmap=plt.cm.hot, interpolation='bilinear')
#    cbar = plt.colorbar(im)
#    cbar.ax.tick_params(labelsize = tamnum)
##    cbar.set_ticks(tick_locations)
#    cbar.set_label(r'Im($\phi$)')
#    os.chdir(path_save)
#    plt.tight_layout()
#    plt.savefig('Im_potential_pole_aprox_log' + label1 + '.png', format='png')
#    