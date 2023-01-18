
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

#%%

save_graphs = 1 #guardar los graficos 2D del campo

plot_xy_3D = 1
plot_yz_3D = 0
plot_xz_3D = 0

plot_vs_Energy = 0

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/potential_field','')
path_save = path_basic + '/' + 'many_electric_field_x'

if save_graphs==1:    
    if not os.path.exists(path_save):
        print('Creating folder to save graphs')
        os.mkdir(path_save)

err = 'fieldE_direct_numerical.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from many_electric_field_x_formula import electric_field_function
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
#v = c/int_v
b = -0.01 #electron position in z < 0 (arriba del plano)
#omega = 0.7*1e12
cota = 1200 #nanometros
epsi1,epsi2 = 1,1
hbmu,hbgama = 0.3,0.0001
zp = 0.05

#xD = 0.01
yD = 0
zD = 0

px = 0
py = 0
pz = 1

E0 = 0.1 #eV 
omegac0 = E0/aux   

title1 = r'v = c/%i $\mu$m/s, b = %inm, E = %.2feV' %(int_v,b*1e3,E0) 
title2 = r'$\hbar\mu$ = %.2feV, $\hbar\gamma$ = %.4feV, $z_p$=%inm' %(hbmu,hbgama,zp*1e3) 
title4 = r'px = %i, py = %i, pz = %i' %(px,py,pz)
title5 = r'yD = %inm, zD = %inm' %(yD*1e3,zD*1e3)
labelp = r'_px%i_py%i_pz%i_E%.4f' %(px,py,pz,E0)
labelx = r'E [eV]'

list_xD = np.linspace(-1000,1000,11)
# list_xD = [10]

if plot_xy_3D == 1 :    #omega = 0.7*1e12
    z0 = 0
    title3A = r'z = %inm' %(z0*1e3)

    title = title1 + '\n' + title2 + '\n' + title3A + ', ' + title4 + '\n'  + title5
    label1 = '_xy_3D' + labelp 
    labelx,labely = 'x [nm]', 'y [nm]'

    def function_real(x_nano,y_nano):
        
        x = x_nano*1e-3 
        y = y_nano*1e-3
        
        rta = electric_field_function(omegac0,epsi1,epsi2,hbmu,hbgama,x,y,z0,list_xD,yD,zD,b,zp,px,py,pz)
        return rta.real
    
    def function_imag(x_nano,y_nano):
        
        x = x_nano*1e-3 
        y = y_nano*1e-3
        
        rta = electric_field_function(omegac0,epsi1,epsi2,hbmu,hbgama,x,y,z0,list_xD,yD,zD,b,zp,px,py,pz)
        return rta.imag

elif plot_yz_3D == 1 :    #omega = 0.7*1e12
    x0 = 0
    title3B = r'x = %inm' %(x0*1e3)

    title = title1 + '\n' + title2 + '\n' + title3B + ', ' + title4 + '\n'  + title5
    label1 = '_yz_3D' + labelp
    labelx,labely = 'y [nm]', 'z [nm]'

    def function_real(y_nano,z_nano):
        y = y_nano*1e-3 
        z = z_nano*1e-3        
        
        rta = electric_field_function(omegac0,epsi1,epsi2,hbmu,hbgama,x0,y,z,list_xD,yD,zD,b,zp,px,py,pz)
        return rta.real

    def function_imag(y_nano,z_nano):
        y = y_nano*1e-3 
        z = z_nano*1e-3        
        
        rta = electric_field_function(omegac0,epsi1,epsi2,hbmu,hbgama,x0,y,z,list_xD,yD,zD,b,zp,px,py,pz)
        return rta.imag
    
elif plot_xz_3D == 1 :    #omega = 0.7*1e12
    y0 = 0
    title3C = r'y = %inm' %(y0*1e3)

    title = title1 + '\n' + title2 + '\n' + title3C + ', ' + title4 + '\n'  + title5
    label1 = '_xz_3D' + labelp
    labelx,labely = 'x [nm]', 'z [nm]'

    def function_real(x_nano,z_nano):
        x = x_nano*1e-3 
        z = z_nano*1e-3
        
        rta = electric_field_function(omegac0,epsi1,epsi2,hbmu,hbgama,x,y0,z,list_xD,yD,zD,b,zp,px,py,pz)
        return rta.real

    def function_imag(x_nano,z_nano):
        x = x_nano*1e-3 
        z = z_nano*1e-3
        
        rta = electric_field_function(omegac0,epsi1,epsi2,hbmu,hbgama,x,y0,z,list_xD,yD,zD,b,zp,px,py,pz)
        return rta.imag

elif plot_vs_Energy == 1: 
    x0 = 0
    y0 = 0
    z0 = 0
    title3D = r'x = %inm, y = %inm, z = %inm' %(x0*1e3,y0*1e3,z0*1e3)

    title = title1 + '\n' + title2 + '\n' + title3D + ', ' + title4 + '\n'  + title5
    label1 = '_vs_E' + labelp
    labelx = 'E [eV]'
    labely = r'$|E|^2$'
    list_OmegaTHz = np.linspace(0.01,2.01,251)    
    listE = np.linspace(0.001, 10, 1001)

    def function_real(energy0):
 
        
        omegac0 = energy0/aux         
        
        rta = electric_field_function(omegac0,epsi1,epsi2,hbmu,hbgama,x0,y0,z0,list_xD,yD,zD,b,zp,px,py,pz)
        return rta    
    

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
    N = 100
    list_x = np.linspace(-cota,cota,N)
    list_y = np.linspace(-cota,cota,N)
    X, Y = np.meshgrid(list_x, list_y)
    limits = [min(list_x) , max(list_x), min(list_y) , max(list_y)]
    
    aux_list = np.array(np.ones(N))
    
    f1_real = np.vectorize(function_real)
    Z1_real = f1_real(X, Y)
   
    print('Graficar el External field' + label1)
    
    graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    vmin, vmax = np.min(Z1_real), np.max(Z1_real)
    pcm = plt.pcolormesh(X, Y, Z1_real,
                          norm=colors.SymLogNorm(linthresh=0.01, linscale=0.01,
                                                vmin=int(vmin), vmax=int(vmax)),cmap='RdBu_r')
    
#    im = plt.imshow(Z1_real, extent = limits, cmap=plt.cm.hot, interpolation='bilinear')
    try: 
        maxlog=int(np.ceil( np.log10( np.abs(vmax) )))
        minlog=int(np.ceil( np.log10( np.abs(vmin) )))
        
        if vmin < 0 :
              tick_locations = ( [-(10.0**x) for x in np.linspace(minlog,-1,minlog+2)] 
                                + [0] 
                                + [(10.0**x) for x in np.linspace(-1,maxlog,maxlog+2)] )
        else:
              tick_locations = ( [(10.0**x) for x in np.linspace(minlog,maxlog,maxlog + np.abs(minlog) + 1) ])    
            
    
        cbar = plt.colorbar(pcm)
        cbar.ax.tick_params(labelsize = tamnum)
        cbar.set_ticks(tick_locations)
        cbar.set_label('Re(Ex)',fontsize=tamlegend,labelpad = 1)
        
        if save_graphs == True:
            os.chdir(path_save)
            plt.tight_layout()
            plt.savefig('Re_many_Efield_log' + label1 + '.png', format='png')
    except OverflowError:
        print('no se puede graficar en escala logaritmica')

    graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    im = plt.imshow(Z1_real, extent = limits, cmap='RdBu_r', interpolation='bilinear')
    cbar = plt.colorbar(im)
    cbar.ax.tick_params(labelsize = tamnum)
    cbar.set_label('Re(Ex)',fontsize=tamlegend,labelpad = 1)
    
    if save_graphs == True:
        os.chdir(path_save)
        plt.tight_layout()
        plt.savefig('Re_many_Efield' + label1 + '.png', format='png')
    
    f2_imag = np.vectorize(function_imag)
    Z2_imag = f2_imag(X, Y)
       
    print('Graficar el External field' + label1)
    
    graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    im = plt.imshow(Z2_imag, extent = limits, cmap='RdBu_r', interpolation='bilinear')
    cbar = plt.colorbar(im)
    cbar.ax.tick_params(labelsize = tamnum)
    cbar.set_label('Im(Ex)',fontsize=tamlegend,labelpad = 1)
    
    for xD in list_xD:
        if plot_xy_3D == 1:
            plt.plot(xD*aux_list, list_y,'--',color = 'green')
            plt.plot(list_x, yD*aux_list,'--',color = 'green')
        elif plot_yz_3D == 1:    
            plt.plot(yD*aux_list, list_y,'--',color = 'green')
            plt.plot(list_x, zD*aux_list,'--',color = 'green')
        elif plot_xz_3D == 1:    
            plt.plot(xD*aux_list, list_y,'--',color = 'green')
            plt.plot(list_x, zD*aux_list,'--',color = 'green')
    
    
    if save_graphs == True:
        os.chdir(path_save)
        plt.tight_layout()
        plt.savefig('Im_many_Efield' + label1 + '.png', format='png')        
        
else:
    listy = function_real(listE)

    graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    plt.plot(listE,listy,'.',ms = ms,color = 'lightseagreen')
    # plt.yscale('log')
    plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
    if save_graphs==1:
        plt.tight_layout()
        os.chdir(path_save)
        plt.savefig( 'Re_many_Efield' + label1 + '.png', format='png')    

#%%
