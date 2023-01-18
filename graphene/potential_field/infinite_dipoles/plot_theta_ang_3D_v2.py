
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

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/potential_field/infinite_dipoles','')
path_save = path_basic + '/' + 'theta_ang'
if not os.path.exists(path_save):
    print('Creating folder to save graphs')
    os.mkdir(path_save)

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

epsi1,epsi2 = 1,1
hbmu,hbgama = 0.3,0.0001

#Ndipoles = 50
#list_xD = np.linspace(-0.01,0.01,Ndipoles)


int_v = 10
Nmax = 2
#def omega_n_THz(int_v,a_nm,Nmax):  ## omega puede valer esto o menos para que se generen plasmones
#    
#    a = a_nm*1e-3
#    
#    aux = alfac*int_v*hbmu/(hb)
#    rta = aux + 0.5*np.sqrt((2*aux)**2 + 16*alfac*c*hbmu*Nmax*np.pi/(a*hb))
#    
#    return rta*1.5

title = r'v = c/%i, n = %i' %(int_v,Nmax) 
labelx = r'$\omega/\omega_{\rm D}$'  
labely = r'a$\omega/v$'

def theta(E_meV,a_micros):
    

#    
#    omega_n = omega_n_THz(int_v,a_nm,Nmax)
#    omegac = omega_n/c
#    print('omega/c:', omegac)
    
 #   omegac = omegac - 0.5 ## ver si funciona con una freq menor 


    E = E_meV*1e-3  
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    alfa_p = 1j*(epsi1 + epsi2)/(cond)
    omegac = E/(hb*c)
    kp = alfa_p*omegac    
    
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
    a = a_micros
    
    theta0 = np.arccos((omegac*int_v - 2*np.pi*Nmax/a)/np.real(kp))
#    print(theta0)
    return theta0*180/np.pi


    
#%%
    
    
tamfig = [2.5, 2]
tamletra = 8
tamtitle  = 8
tamnum = 6
tamlegend = 6
labelpady = -0.5
labelpadx = 0.5
pad = 2.5
mk = 1
ms = 1
hp = 0.5
length_marker = 1.5
dpi = 500


def graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad):
    plt.figure(figsize=tamfig)
    plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
    plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
    plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in", pad = pad)
#    plt.title(title,fontsize=int(tamtitle*0.9))

    return   

#%% 
N2 = 800
if Nmax == 2:
    listx = np.linspace(150,400,N2)

elif Nmax == 1:
    listx = np.linspace(100,350,N2)
#    listx = np.linspace(200,350,N2)
#    listx = np.linspace(43.5,150,N2)
    listx = np.linspace(20,350,N2)
    listx = np.linspace(30,100,N2)
    listy = np.array(np.linspace(2600,10500,N2))*1e-3
elif Nmax == 0:
    listx = np.linspace(40,60,N2)
    listx = np.linspace(43.5,150,N2)
    listy = np.array(np.linspace(200,15000,N2))*1e-3







if Nmax == 2:
    list_EmeV = np.linspace(150,400,N2)
    list_EmeV = listx
    list_EmeV = np.linspace(30,150,N2)

    list_EmeV = np.linspace(55,100,N2)
    list_a_micros = np.array(np.linspace(2000,10500,N2))*1e-3

elif Nmax == 1:
    list_EmeV = np.linspace(40,400,N2)
    list_EmeV = np.linspace(30,150,N2)
    list_EmeV = np.linspace(50,350,N2)
    list_a_micros = np.array(np.linspace(12000,18000,N2))*1e-3

    list_EmeV = np.linspace(30,100,N2)
    list_EmeV = np.linspace(48,86,N2)
    list_a_micros = np.array(np.linspace(2000,10500,N2))*1e-3

#
#    list_EmeV = np.linspace(40,76,N2)
#    list_a_micros = np.array(np.linspace(2000,10500,N2))*1e-3

#    list_EmeV = listx
elif Nmax == 0:
    list_EmeV = np.linspace(30,50,N2)
    list_EmeV = listx
    list_EmeV = np.linspace(43.25,150,N2)

listx = list_EmeV
listy = list_a_micros


X, Y = np.meshgrid(listx, listy)
f_num = np.vectorize(theta)
Z_num = f_num(X, Y)




#%% 


#list_angle_lim = []
#list_omega_omegaD_lim = []
#list_a_omega_v = []
#
#for a_micros in listy :
#    list_angle_energy = []
#    for ind in range(len(listy)):
#        energy0_meV = listx[ind]
#
#        angle =  theta(energy0_meV,a_micros)
#
#        if np.isnan(angle) == True :
#            
#            list_angle_energy.append(energy0_meV)
#      
#    omega_omegaD_lim = np.max(list_angle_energy)*1e-3/hbmu
#    
#    list_omega_omegaD_lim.append(omega_omegaD_lim)
#    omegac =  np.max(list_angle_energy)*1e-3/(hb*c)
#    list_a_omega_v.append(a_micros*omegac*int_v)
#    
#    
#    list_angle_lim.append(np.max(list_angle_energy))


#%%
#ticks_z = [0,0.2,0.4,0.6,0.8,1]

omega_D = hbmu ## sin el hbar porque se cancela con el listx 

listx_2 = np.array(listx)*1e-3/hbmu ## sin el hbar porque se cancela con el listx 
#list_angle_lim_2_ejex = np.array(list_angle_lim)*1e-3/hbmu ## sin el hbar porque se cancela con el listx 

listy_2 = []
list_angle_lim_2_ejey = []

for k in range(len(listy)):
    a_micros = listy[k]
    omegac = listx[k]*1e-3/(hb*c)
    rta = a_micros*omegac*int_v
    listy_2.append(rta)
    
#    omegac2 = list_angle_lim[k]*1e-3/(hb*c)
#    rta2 = a_micros*omegac2*int_v
#    list_angle_lim_2_ejey.append(rta2)
    
    
    
#%% 
import matplotlib.colors as mcolors
    
list_omega_omegaD = np.array(list_EmeV)*1e-3/hbmu ## sin el hbar porque se cancela con el listx 
list_omegac = np.array(list_EmeV)*1e-3/(hb*c)
list_a_omega_v = np.array(list_a_micros)*np.array(list_omegac)*int_v    

labelz = r'$\theta$ ' +  '(' + u"\N{DEGREE SIGN}" + ')'


limits = [np.min(list_omega_omegaD) , np.max(list_omega_omegaD), np.min(list_a_omega_v) , np.max(list_a_omega_v)]
graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#for x in [x1,x2,x3,x4]:
#    plt.plot(ejey_aux1,x*np.ones(10),'--',color = 'grey' )
#norm = mcolors.DivergingNorm(vmin=np.min(Z_num), vmax = np.max(Z_num),vcenter = 0.5)



cmap = plt.cm.RdBu  # define the colormap
# extract all colors from the .jet map
cmaplist = [cmap(i) for i in range(cmap.N)]
import matplotlib as mpl
# create the new map
cmap = mpl.colors.LinearSegmentedColormap.from_list(
    'Custom cmap', cmaplist, cmap.N)

 
Z_min = np.round(np.min(Z_num))
Z_max = np.round(np.max(Z_num))
dist = round((Z_max -Z_min)/4)
if Nmax == 1: 
    ticks = np.array(np.linspace(Z_min, Z_max, 11  ))
else:
    ticks = np.array(np.linspace(Z_min, Z_max, int((Z_max - Z_min)/20)  + 6 ))
n_color = 2*(dist - 1)-1
ticks_2 = []
for tt in ticks: 
    ticks_2.append(np.round(tt,2))


# define the bins and normalize
bounds = np.linspace(np.min(Z_num), np.max(Z_num)  , dist + 1)

norm = mpl.colors.BoundaryNorm(bounds, cmap.N)


im = plt.imshow(Z_num, extent = limits, cmap=cmap, aspect='auto', interpolation = 'none',origin = 'lower', norm=norm) 
#plt.plot(maxis,0.18*np.ones(len(maxis)),'o',color = 'green' )
cbar = plt.colorbar(im, fraction=0.046, pad=0.04, orientation = 'vertical')
#plt.clim(np.min(ticks_z),np.max(ticks_z))
cbar.set_ticks(ticks_2)
cbar.ax.set_yticklabels(ticks_2)
cbar.ax.tick_params(labelsize = tamnum, width=0.5, length = 0,pad = 1.5)

#plt.plot(list_omega_omegaD_lim,listy_2,'-',color = 'blue')
cbar.outline.set_visible(False)
cbar.ax.set_title(labelz,fontsize=tamlegend)

plt.tight_layout()
os.chdir(path_save)
plt.savefig( 'theta_N%i' %(Nmax) + '.png', format='png',bbox_inches='tight',pad_inches = 0,dpi = dpi)   

#%%         
        



import matplotlib.colors as mcolors
limits = [np.min(listx) , np.max(listx), np.min(listy) , np.max(listy)]
graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#for x in [x1,x2,x3,x4]:
#    plt.plot(ejey_aux1,x*np.ones(10),'--',color = 'grey' )
#norm = mcolors.DivergingNorm(vmin=np.min(Z_num), vmax = np.max(Z_num),vcenter = 0.5)
im = plt.imshow(Z_num, extent = limits, cmap=plt.cm.hot, aspect='auto', interpolation = 'none',origin = 'lower') 
#plt.plot(maxis,0.18*np.ones(len(maxis)),'o',color = 'green' )
cbar = plt.colorbar(im, fraction=0.046, pad=0.04, orientation = 'vertical')
#plt.clim(np.min(ticks_z),np.max(ticks_z))
cbar.ax.tick_params(labelsize = tamnum, width=0.5, length = 2,pad = 1.5)

#plt.plot(list_angle_lim,listy,'-',color = 'blue')

cbar.set_label(r'$\theta$ ' +  '(' + u"\N{DEGREE SIGN}" + ')',fontsize=tamlegend,labelpad = 2)
plt.tight_layout()
os.chdir(path_save)





    