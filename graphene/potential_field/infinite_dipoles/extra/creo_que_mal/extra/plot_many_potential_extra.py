
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

plot_vs_xy = 1
plot_vs_yz = 0
plot_vs_Energy = 0

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/potential_field/many_potential','')
path_save = path_basic + '/' + 'extra'

if save_graphs==1:    
    if not os.path.exists(path_save):
        print('Creating folder to save graphs')
        os.mkdir(path_save)

err = 'many_potential.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from potential_many_dipoles_num_extra import phi_many_dipoles_num_extra
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
#omega = 0.7*1e12
cota = 1200 #nanometros
epsi1,epsi2 = 1,1
hbmu,hbgama = 0.3,0.0001
zp = 0.05 #micrometros 

px = 1
py = 0
pz = 0

omega0THz = 0.6
omega0 = omega0THz*1e12
lambdda0 = c/(2*np.pi*omega0)
E0 =omega0THz*1e12*hb  #eV
omegac0 = E0/aux   

N = 10
a = 0.5 

# list_xD = np.array([10])
    
x0 = 0
y0 = 0
z0 = 0

listky = np.linspace(0.1,800,1000)

title1 = r'v = c/%i $\mu$m/s, a = %inm, N = %i, E = %.2fmeV' %(int_v,a*1e3,N,E0*1e3)
title2 = r'$\hbar\mu$ = %.2feV, $\hbar\gamma$ = %.4feV, $z_p$=%inm' %(hbmu,hbgama,zp*1e3) 
title3 = r'x = %inm, y = %inm, z = %inm' %(x0*1e3,y0*1e3,z0*1e3) 
title4 = r'px = %i, py = %i, pz = %i' %(px,py,pz)

labelp = r'_px%i_py%i_pz%i' %(px,py,pz)
labelx = r'ky'
label1 = r'_yz_3D' + labelp  + '_E%.4f' %(E0)
labelx, labely1, labely2 = 'ky', 'Re(function inside integral)', 'Im(function inside integral)'

def function_real(ky):
    rta = phi_many_dipoles_num_extra(omegac0,epsi1,epsi2,hbmu,hbgama,x0,y0,z0,a,N,zp,int_v,ky,px,py,pz)
    return rta.real

def function_imag(ky):  
    rta = phi_many_dipoles_num_extra(omegac0,epsi1,epsi2,hbmu,hbgama,x0,y0,z0,a,N,zp,int_v,ky,px,py,pz)
    return rta.imag

title = title1 + '\n' + title2 + ', '  + title3 + '\n' + title4 

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

import seaborn as sns
sns.set()

listy_re = []
listy_im = []

for ky in listky:
    rta_re = function_real(ky)
    rta_im = function_imag(ky)
    listy_re.append(rta_re)
    listy_im.append(rta_im) 
#    print(j)
#    j = j + 1
    
graph(title,labelx,labely1,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(listky,listy_re,'.-',ms = ms,color = 'lightseagreen')
if save_graphs==1:
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'inside_integral_Re' + label1 + '.png', format='png')    


graph(title,labelx,labely2,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(listky,listy_im,'.-',ms = ms,color = 'lightseagreen')
if save_graphs==1:
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'inside_integral_Im' + label1 + '.png', format='png')    


# plot 2D 
    
def function_real2D(ky,N):
    rta = phi_many_dipoles_num_extra(omegac0,epsi1,epsi2,hbmu,hbgama,x0,y0,z0,a,N,zp,int_v,ky,px,py,pz)
    return rta.real

def function_imag2D(ky,N):  
    rta = phi_many_dipoles_num_extra(omegac0,epsi1,epsi2,hbmu,hbgama,x0,y0,z0,a,N,zp,int_v,ky,px,py,pz)
    return rta.imag    

Nlist = 101
list_x = np.linspace(1,200,Nlist)
list_y = np.linspace(1,101,Nlist)
X, Y = np.meshgrid(list_x, list_y)
limits = [min(list_x) , max(list_x), min(list_y) , max(list_y)]

aux_list = np.array(np.ones(Nlist))

f1_real = np.vectorize(function_real2D)
Z1_real = f1_real(X, Y)

f2_real = np.vectorize(function_imag2D)
Z2_real = f2_real(X, Y)
   
print('Graficar el External field' + label1)

graph(title,'ky','N',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
im = plt.imshow(Z1_real, extent = limits, cmap='RdBu_r', interpolation='bilinear')
cbar = plt.colorbar(im)
cbar.ax.tick_params(labelsize = tamnum)
cbar.set_label(r'Re(...)',fontsize=tamlegend,labelpad = 1)
if save_graphs == True:
    os.chdir(path_save)
    plt.tight_layout()
    plt.savefig('2Dinside_integral_Re' + label1 + '.png', format='png')

graph(title,'ky','N',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
im = plt.imshow(Z2_real, extent = limits, cmap='RdBu_r', interpolation='bilinear')
cbar = plt.colorbar(im)
cbar.ax.tick_params(labelsize = tamnum)
cbar.set_label(r'Im(...)',fontsize=tamlegend,labelpad = 1)
if save_graphs == True:
    os.chdir(path_save)
    plt.tight_layout()
    plt.savefig('2Dinside_integral_Im' + label1 + '.png', format='png')

    
#%%
