
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
bigger_than_Lp = 1

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/potential_field/potential_and_electric_field/sum_dipoles/potential_and_electric_field_simple_version','')
path_save = path_basic + '/' + 'sum_dipoles_phi3D'
if not os.path.exists(path_save):
    print('Creating folder to save graphs')
    os.mkdir(path_save)

err = 'many_potential_integrals.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from potential_final_version import potential_final_num, potential_final_ana, potential_final_pole_aprox
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
#b = -0.01

z0 = zp 

#Ndipoles = 50
#list_xD = np.linspace(-0.01,0.01,Ndipoles)

px = 0
py = 0
pz = 1


int_v, E0 = 10, 60# meV 
#int_v, E0 = 20, 21

#int_v = 30
#E0 = 30   #----------> sin ver plasmones
omegac0 = E0*1e-3/aux

def k_SP(E_mev,hbmu,hbgama):
    num = 2/alfac 
    omegac = E_mev*1e-3/aux
    cte = E_mev*1e-3 + 1j*hbgama/(4*hbmu)
    
    return omegac*num*cte

lambda_SP = 2*np.pi/np.real(k_SP(E0,hbmu,hbgama))
Lp = 1*1e3/(2*np.imag(k_SP(E0,hbmu,hbgama)))  ## nanometros 
print(k_SP(E0,hbmu,hbgama))
#a_min = np.real(lambda_SP)*Nmax/(int_v - 1)
#a_max = np.real(lambda_SP)*Nmax/(int_v + 1)

lambda_SP_nm = np.real(lambda_SP)*1e3
cota_x = int(lambda_SP_nm*10)
#cota_x =  int(Lp*1.5)

#a = 6.5
#a = 0.05
#a = np.mean([a_min,a_max])
a = 0.05     
#a = int(np.real(lambda_SP)*1.5)
#a = int(np.real(lambda_SP)*0.5)
#a = 0.05
N = 3

#a = a_min

#theta0 = np.arccos(int_v-np.real(lambda_SP)*Nmax/a)

title2 = r'v = c/%i, a  = %inm, N = %i' %(int_v,a*1e3,N) 
#title3 = r'px=%i, py=%i, pz=%i, a = %inm' %(px,py,pz,a*1e3)
#title4 = r'Ndip = %i' %(Ndip)

if py == 0 and pz ==0:
    title4 = r'z = $z_p$=%inm, px = %i, Re(kp) = %.1f 1/$\mu$m' %(zp*1e3,px,np.real(k_SP(E0,hbmu,hbgama)))
elif py == 0 and px ==0:
    title4 = r'z = $z_p$=%inm, pz = %i,Re(kp) = %.1f 1/$\mu$m' %(zp*1e3,pz,np.real(k_SP(E0,hbmu,hbgama)))

title2 = title2 + ', ' + 'Lp = %inm' %(Lp) 

labelp = r'_px%i_py%i_pz%i_a%inm_Ndip%i_intv%i_E%imeV' %(px,py,pz,a*1e3,N,int_v,E0)

labelx = r'x [$\mu$m]'  
labely = r'y [$\mu$m]'
title = title2 + '\n' + title4 + ', E = %i meV' %(E0)
label1 = 'vs_y_E0%i' %(E0) + labelp
listy = np.linspace(50,5000,100)
listx = np.linspace(-2500,2500,100)
listy =listx


if bigger_than_Lp == 1:
    listy = np.array(np.linspace(-Lp*1.5,Lp*1.5,100))
    listx = listy
    
else:
    listx = np.array(np.linspace(-3500,3500,100))  ## micrones
    listy = listx

#%%

def function_real_num(x_nano,y_nano):
    y = y_nano*1e-3 
    x = x_nano*1e-3 - a*j 
    R = np.sqrt(x**2 + y**2)
    phi  = np.arctan2(y,x)
    

    rta = potential_final_num(omegac0,epsi1,epsi2,hbmu,hbgama,phi,R,z0,zp,px,py,pz)
    rta0 = rta*np.exp(a*j*omegac0*int_v)
        
    # print(rta)
    return rta0.real

def function_imag_num(x_nano,y_nano):
    y = y_nano*1e-3 
    x = x_nano*1e-3 - a*j 
    R = np.sqrt(x**2 + y**2)
    phi  = np.atan2(y,x)

    rta = potential_final_num(omegac0,epsi1,epsi2,hbmu,hbgama,phi,R,z0,zp,px,py,pz)
    rta0 = rta*np.exp(a*j*omegac0*int_v)
    return rta0.imag

#%%
    
def function_real_ana(x_nano,y_nano):
    y = y_nano*1e-3 

    term0 = 0
    for j in range(0,N):


        x = x_nano*1e-3 - a*j  
        R = np.sqrt(x**2 + y**2)
        phi  = np.arctan2(y,x)


        rta = potential_final_ana(omegac0,epsi1,epsi2,hbmu,hbgama,phi,R,z0,zp,px,py,pz)
        rta0 = rta*np.exp(a*j*omegac0*int_v)
        term0 = rta0 + term0
        
    # print(rta)
    return term0.real

def function_imag_ana(x_nano,y_nano):
    y = y_nano*1e-3 
    
    term0 = 0
    for j in range(0,N):
        x = x_nano*1e-3 - a*j 

        R = np.sqrt(x**2 + y**2)
        phi  = np.arctan2(y,x)
   
        rta = potential_final_ana(omegac0,epsi1,epsi2,hbmu,hbgama,phi,R,z0,zp,px,py,pz)
        rta0 = rta*np.exp(a*j*omegac0*int_v)     
        term0 = rta0 + term0
    
    return term0.imag

#%%
    
def function_real_pole_aprox(x_nano,y_nano):
    y = y_nano*1e-3 
    x = x_nano*1e-3 - a*j 
    R = np.sqrt(x**2 + y**2)
    phi  = np.arctan2(y,x)

    rta = potential_final_pole_aprox(omegac0,epsi1,epsi2,hbmu,hbgama,phi,R,z0,zp,px,py,pz)
    rta0 = rta*np.exp(a*j*omegac0*int_v)      
    # print(rta)
    return rta0.real

def function_imag_pole_aprox(x_nano,y_nano):
    
    y = y_nano*1e-3 
    x = x_nano*1e-3 - a*j 
    R = np.sqrt(x**2 + y**2)
    phi  = np.arctan2(y,x)
    rta = potential_final_pole_aprox(omegac0,epsi1,epsi2,hbmu,hbgama,phi,R,z0,zp,px,py,pz)
    rta0 = rta*np.exp(a*j*omegac0*int_v)      
    return rta0.imag

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
X, Y = np.meshgrid(listx, listy)
limits = [min(listx) , max(listx), min(listy) , max(listy)]

#aux_list = np.array(np.ones(100))

f_real_ana = np.vectorize(function_real_ana)
Z_real_ana = f_real_ana(X, Y)

f_imag_ana = np.vectorize(function_imag_ana)
Z_imag_ana = f_imag_ana(X, Y)

print('Graficar el External field' + label1)

# if plot_xz_3D == 1 and px == 1: 
graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
vmin, vmax = np.min(Z_real_ana), np.max(Z_real_ana)
im = plt.imshow(Z_real_ana, extent = limits, cmap='coolwarm')
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
im = plt.imshow(Z_imag_ana, extent = limits, cmap='coolwarm')
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

if plot_num == 1:
    f_real_num = np.vectorize(function_real_num)
    Z_real_num = f_real_ana(X, Y)
    
    f_imag_num = np.vectorize(function_imag_num)
    Z_imag_num = f_imag_ana(X, Y)
       
    print('Graficar el External field' + label1)
    
    # if plot_xz_3D == 1 and px == 1: 
    graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    vmin, vmax = np.min(Z_real_num), np.max(Z_real_num)
#    pcm = plt.pcolormesh(X, Y, Z_real_num,
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
#    #    im = plt.imshow(Z1_real, extent = limits, cmap=plt.cm.hot, interpolation='bilinear')
#    cbar = plt.colorbar(pcm)
    im = plt.imshow(Z_real_num, extent = limits, cmap=plt.cm.hot, interpolation='bilinear')
    cbar = plt.colorbar(im)
    cbar.ax.tick_params(labelsize = tamnum)
#    cbar.set_ticks(tick_locations)
    cbar.set_label(r'Re($\phi$)')
    os.chdir(path_save)
    plt.tight_layout()
    plt.savefig('Re_potential_num_log' + label1 + '.png', format='png')
    
    
    
    
    graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    vmin, vmax = np.min(Z_imag_num), np.max(Z_imag_num)
#    pcm = plt.pcolormesh(X, Y, Z_imag_num,
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
    im = plt.imshow(Z_imag_num, extent = limits, cmap=plt.cm.hot, interpolation='bilinear')
    cbar = plt.colorbar(im)

#    cbar = plt.colorbar(pcm)
    cbar.ax.tick_params(labelsize = tamnum)
#    cbar.set_ticks(tick_locations)
    cbar.set_label(r'Im($\phi$)')
    os.chdir(path_save)
    plt.tight_layout()
    plt.savefig('Im_potential_num_log' + label1 + '.png', format='png')
        
    
#    del Z_real_num, Z_imag_num

#%%


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

    im = plt.imshow(Z_real_pole_aprox, extent = limits, cmap=plt.cm.hot, interpolation='bilinear')
    cbar = plt.colorbar(im)
    cbar.ax.tick_params(labelsize = tamnum)
#    cbar.set_ticks(tick_locations)
    cbar.set_label(r'Re($\phi$)')
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
    im = plt.imshow(Z_imag_pole_aprox, extent = limits, cmap=plt.cm.hot, interpolation='bilinear')
    cbar = plt.colorbar(im)
    cbar.ax.tick_params(labelsize = tamnum)
#    cbar.set_ticks(tick_locations)
    cbar.set_label(r'Im($\phi$)')
    os.chdir(path_save)
    plt.tight_layout()
    plt.savefig('Im_potential_pole_aprox_log' + label1 + '.png', format='png')
    