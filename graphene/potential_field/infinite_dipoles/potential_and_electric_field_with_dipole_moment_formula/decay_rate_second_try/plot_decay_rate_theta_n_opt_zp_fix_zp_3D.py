
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
from scipy.interpolate import interp1d
#import seaborn as sns
#sns.set()
#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/potential_field/infinite_dipoles/potential_and_electric_field_with_dipole_moment_formula/decay_rate_second_try','')
path_save = path_basic + '/' + 'decay_rate_theta_n'
if not os.path.exists(path_save):
    print('Creating folder to save graphs')
    os.mkdir(path_save)

err = 'decay_rate_theta_n.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from decay_rate_theta_n import decay_rate_theta_inf_dipoles_ana_res,decay_rate_theta_inf_dipoles_ana_res_div_gamma0
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

epsi1, epsi2 = 1,1
hbmu, hbgama = 0.3,0.0001

print('Definir parametros del problema')


b = - 0.01

int_v = 20
int_v = 10

Nmax = 1

if Nmax!= 0:
    labelz =  r'$\Gamma_{\rm SP}/(\Gamma_{\rm 0}$N)'
else:
    labelz = r'$\Gamma_{\rm SP}/(\Gamma_{\rm 0}$N)'

tabla = np.loadtxt('zp_optimum_for_decay_rate_graphene_resonance_b-10nm.txt', delimiter='\t', skiprows=1)
tabla = np.transpose(tabla)
[listx,listy,listz] = tabla ## energy and zp and lambda_p 

zp_nano = listy[-20]

omegac0_1 = np.max(listx)*1e-3/(c*hb)
lambda_SP_1 = 2*np.pi/omegac0_1

omegac0_2 = np.min(listx)*1e-3/(c*hb)
lambda_SP_2 = 2*np.pi/omegac0_2


a_min = np.real(lambda_SP_1)*Nmax/(int_v - 1)
a_max = np.real(lambda_SP_2)*Nmax/(int_v + 1)


f1 = interp1d(listx, listy)
N = 500
lim1,lim2 = 10, -60
#lim1,lim2 = 16,-80
listx = np.linspace(listx[lim1], listx[lim2], N)


list_a_micros = np.array(np.linspace(200,15000,N))*1e-3
if Nmax == 2:
    list_EmeV = np.linspace(150,400,N)
    list_EmeV = listx
    list_EmeV = np.linspace(30,150,N)

    list_EmeV = np.linspace(60,100,N)
    list_a_micros = np.array(np.linspace(1000,9500,N))*1e-3

elif Nmax == 1:
    list_EmeV = np.linspace(40,400,N)
    list_EmeV = np.linspace(30,150,N)
    list_EmeV = np.linspace(50,350,N)
    list_a_micros = np.array(np.linspace(12000,18000,N))*1e-3

    list_EmeV = np.linspace(30,100,N)
    list_EmeV = np.linspace(48,76,N)

    list_EmeV = np.linspace(48,86,N)
    list_a_micros = np.array(np.linspace(3000,10500,N))*1e-3
    list_a_micros = np.array(np.linspace(2000,10500,N))*1e-3

#    list_EmeV = listx
elif Nmax == 0:
    list_EmeV = np.linspace(30,50,N)
    list_EmeV = listx
    list_EmeV = np.linspace(43.25,150,N)
#    list_EmeV = np.linspace(41,53,N)
#    list_a_micros = np.array(np.linspace(12000,14000,N))*1e-3

#lista de maximos [43.75561981083406, 50.88340900504751, 56.5714557053054, 61.45097706413012, 65.78438920110963] en Emev para n = 0,1,2,3,4
    
labelp = r'_zp%inm_N%i' %(zp_nano,Nmax)

#%%


def function_real_ana(energy0_meV,a_micros):
                
    omegac0 = energy0_meV*1e-3/(c*hb)  
    a = a_micros
    zp = zp_nano*1e-3
         
    rta = decay_rate_theta_inf_dipoles_ana_res_div_gamma0(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,zp,a,b,Nmax)
    if Nmax!=0:
        final_rta = np.log10(rta)
        print(final_rta)
    else:
        final_rta = np.log10(rta)
        print(final_rta)
    return rta 

#%%

X, Y = np.meshgrid(list_EmeV, list_a_micros)
f_num = np.vectorize(function_real_ana)
Z_num = f_num(X, Y)

list_omega_omegaD = np.array(list_EmeV)*1e-3/hbmu ## sin el hbar porque se cancela con el listx 
list_omegac = np.array(list_EmeV)*1e-3/(hb*c)
list_a_omega_v = np.array(list_a_micros)*np.array(list_omegac)*int_v

#%%

tamfig = [2.5, 2]
tamletra = 9
tamtitle  = 8
tamnum = 7
tamlegend = 6
labelpady = 1.5
labelpadx = 1.5
pad = 2.5
mk = 1
ms = 1
hp = 0.5
length_marker = 1.5
dpi = 500


#%%
    


def graph(labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad):
    plt.figure(figsize=tamfig)
    plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
    plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
    plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in", pad = pad)
#    plt.title(title,fontsize=int(tamtitle*0.9))

    return  
 
#%%

labelx = r'$\omega/\omega_{\rm D}$'  
labely = r'a$\omega/v$'

import matplotlib.colors as mcolors
limits = [np.min(list_omega_omegaD) , np.max(list_omega_omegaD), np.min(list_a_omega_v) , np.max(list_a_omega_v)]
graph(labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#for x in [x1,x2,x3,x4]:
#    plt.plot(ejey_aux1,x*np.ones(10),'--',color = 'grey' )
#norm = mcolors.DivergingNorm(vmin=np.min(Z_num), vmax = np.max(Z_num),vcenter = 0.5)

#plt.plot(maxis,0.18*np.ones(len(maxis)),'o',color = 'green' )




#cmap = plt.cm.Rdbu  # define the colormap
cmap = plt.cm.RdBu  # define the colormap
# extract all colors from the .jet map
cmaplist = [cmap(i) for i in range(cmap.N)]
import matplotlib as mpl
# create the new map
cmap = mpl.colors.LinearSegmentedColormap.from_list(
    'Custom cmap', cmaplist, cmap.N)
if Nmax == 1:
    ticks = [6,6.25,6.5,6.75,7,7.25,7.5,7.75]
    n_color = 21
elif Nmax == 2:
    ticks = [1.5,3.5,4.5,5.5,6.5,7.5,8.5]
    n_color = 14
else:
    Z_min = np.round(np.min(Z_num))
    Z_max = np.round(np.max(Z_num))
    ticks = np.linspace(Z_min,Z_max,8)
    n_color = 14
    
# define the bins and normalize
bounds = np.linspace(np.min(Z_num), np.max(Z_num), n_color)

norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

im = plt.imshow(Z_num, extent = limits, cmap=cmap, aspect='auto', interpolation = 'none',origin = 'lower', norm=norm) 
cbar = plt.colorbar(im, fraction=0.046, pad=0.04, orientation = 'vertical')
#plt.clim(np.min(ticks_z),np.max(ticks_z))

#ax2 = fig.add_axes([0.95, 0.1, 0.03, 0.8])
#cb = plt.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm, spacing='proportional', ticks=bounds, boundaries=bounds, format='%1i')
cbar.ax.set_yticklabels(ticks)
cbar.ax.tick_params(labelsize = tamnum, width=0.5, length = 0,pad = 1.5)

#plt.plot(list_omega_omegaD_lim,listy_2,'-',color = 'blue')
cbar.outline.set_visible(False)
cbar.ax.set_title(labelz,fontsize=tamlegend)
plt.tight_layout()
os.chdir(path_save)
plt.savefig('decay_rate_fix_zp_' + labelp + '_3D.png', format='png',bbox_inches='tight',pad_inches = 0.008,dpi = dpi)



#%%





import matplotlib.colors as mcolors
from matplotlib.colors import SymLogNorm
limits = [np.min(list_EmeV) , np.max(list_EmeV), np.min(list_a_micros) , np.max(list_a_micros)]
graph('Energy (meV)',r'Lattice period ($\mu$m)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#for x in [x1,x2,x3,x4]:
#    plt.plot(ejey_aux1,x*np.ones(10),'--',color = 'grey' )
#norm = mcolors.DivergingNorm(vmin=np.min(Z_num), vmax = np.max(Z_num),vcenter = 0.5)

#plt.plot(maxis,0.18*np.ones(len(maxis)),'o',color = 'green' )




#cmap = plt.cm.Rdbu  # define the colormap
cmap = plt.cm.RdBu  # define the colormap
# extract all colors from the .jet map
cmaplist = [cmap(i) for i in range(cmap.N)]
import matplotlib as mpl
# create the new map
cmap = mpl.colors.LinearSegmentedColormap.from_list(
    'Custom cmap', cmaplist, cmap.N)
if Nmax == 1:
    ticks = [6,6.25,6.5,6.75,7,7.25,7.5,7.75]
    ticks = [10**6,1.5*10**6,0.5*10**7,10**7,1.5*10**7,2*10**7]
    n_color = 21
elif Nmax == 2:
    ticks = [1.5,3.5,4.5,5.5,6.5,7.5,8.5]
    n_color = 14
else:
    Z_min = np.round(np.min(Z_num))
    Z_max = np.round(np.max(Z_num))
    ticks = np.linspace(Z_min,Z_max,8)
    n_color = 14


vmin,vmax = np.min(Z_num), np.max(Z_num)
maxlog=int(np.ceil( np.log10( np.abs(vmax) )))
minlog=int(np.ceil( np.log10( np.abs(vmin) )))
tick_locations = ( [(10.0**x) for x in np.linspace(minlog,maxlog-0.5,maxlog - minlog+2) ])       
tick_locations = [10**6,5*10**6,10**7]
 
# define the bins and normalize
bounds = np.linspace(np.min(Z_num), np.max(Z_num), n_color)

norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

im = plt.imshow(Z_num, extent = limits, cmap=cmap, aspect='auto', interpolation = 'none',origin = 'lower', norm = SymLogNorm(linthresh=0.03, 
                                        linscale=0.03, vmin=int(vmin), vmax=int(vmax))) 

cbar = plt.colorbar(im, fraction=0.046, pad=0.04, orientation = 'vertical')
#plt.clim(np.min(ticks_z),np.max(ticks_z))

#ax2 = fig.add_axes([0.95, 0.1, 0.03, 0.8])
#cb = plt.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm, spacing='proportional', ticks=bounds, boundaries=bounds, format='%1i')
cbar.ax.set_yticklabels(tick_locations)
cbar.set_ticks(tick_locations)
cbar.ax.tick_params(labelsize = tamnum, width=0.5, length = 0,pad = 1.5)

#plt.plot(list_omega_omegaD_lim,listy_2,'-',color = 'blue')
cbar.outline.set_visible(False)
cbar.ax.set_title(labelz,fontsize=tamlegend)
plt.tight_layout()
os.chdir(path_save)
plt.savefig('decay_rate_fix_zp_' + labelp + '_v2_3D.png', format='png',bbox_inches='tight',pad_inches = 0.008,dpi = dpi)





















#list_n = [0,1,2,3,4]   
#
#graph(title,labelx,labely ,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#for n in list_n:
#    
#    list_y_re = []
#
#
#    for ind in range(len(listx_2)):
##        zp = listy_2[ind]
#        x =  listx_2[ind]
##        x = 43 #meV
#        list_y_re.append(function_real_ana(x,n))
#        
#
#    listx_3 = np.array(listx_2)/np.array(listz_2)
#    plt.plot(listx_3,np.array(list_y_re)*1e-8,'-',ms = ms, label = 'n = %i'%(n))
#    
#plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=False,handletextpad=hp, handlelength=1)
##    plt.grid(1)
#plt.tight_layout()
#
##plt.yscale('log')
#os.chdir(path_save)
#plt.savefig('decay_rate_fix_zp_' + labelp + '.png', format='png',bbox_inches='tight',pad_inches = 0.008,dpi = dpi)


#%%

