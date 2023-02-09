
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: leila

funciones dentro de
el campo externo reflejado 
"""

import numpy as np
import sys
import os 
import matplotlib.pyplot as plt

#%%

save_graphs = 1 #guardar los graficos 2D del campo 

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/External_Efield','')
path_save = path_basic + '/' + 'functions_fieldE_ref'

if save_graphs==1:    
    if not os.path.exists(path_save):
        print('Creating folder to save graphs')
        os.mkdir(path_save)

err = 'fieldE_ref_numerical.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from functions_fieldE_ref_numerical import Efield_NUM_QE, Efield_NUM
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
b = -0.01 #electron position in z < 0 (arriba del plano)
#omega = 0.7*1e12

epsi1,epsi2 = 1,1
hbmu,hbgama = 0.3,0.0001
zp = 0.05

E = 0.1
omegaTHz0 = E*1e-12/hb
omegac0 = omegaTHz0*aux2

z0 = 0
x0 = 0
y0 = 0

title1 = r'$\epsilon_1$ = %i, v = c/%i $\mu$m/s, E = %.4feV' %(epsi1,int_v,E) 
label1 = '_E%.4f' %(E)
labelx,labely = r'$\alpha_\parallel$', r'$k$ $x_D$'
title2 = r'$\epsilon_2$ = %i, $\hbar\mu$ = %.2feV, $\hbar\gamma$ = %.4feV, $z_p=$%inm' %(epsi2,hbmu,hbgama,zp*1e3) 
title3 = r'x = %inm, y = %inm, z = %inm, b = %inm' %(x0*1e3,y0*1e3,z0*1e3,b*1e3)
title = title1 + '\n'  + title2 + '\n'  + title3   

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

#%%

cota = 50
N = 100
list_x1 = np.linspace(0.05,0.95,int(N/10)) #alfa paralelo : evitar el 1 
list_x2 = np.linspace(1.95,40.95,int(N*9/10)) #alfa paralelo : evitar el 1
list_x = []
for x in list_x1:
    list_x.append(x)
for x in list_x2:
    list_x.append(x)

if E == 1: 
    list_x = np.linspace(0,20,N)
    list_y = np.linspace(-10,10,N) #x tilde 
elif E == 0.5:
    list_x = np.linspace(0,23,N)
    list_y = np.linspace(-10,10,N) #x tilde 
elif E == 0.1:
    list_x = np.linspace(22,24,N)
    list_y = np.linspace(-0.5,0.5,N) #x tilde 
elif E == 0.01:
    list_x = np.linspace(50,150,N)
    list_y = np.linspace(-50,50,N) #x tilde 
elif E == 0.001:
    list_x = np.linspace(5,45,N)
    list_y = np.linspace(-7.5,7.5,N) #x tilde 


X, Y = np.meshgrid(list_x, list_y)
limits = [min(list_x) , max(list_x), min(list_y) , max(list_y)]

def functionQE_real(alpha_parallel,xD_tilde):
    rta = Efield_NUM_QE(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,b,zp,x0,y0,z0,alpha_parallel,xD_tilde)
    return rta.real

def function_real(alpha_parallel,xD_tilde):
    rta = Efield_NUM(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,b,zp,x0,y0,z0,alpha_parallel,xD_tilde)
    return rta.real

f1_real = np.vectorize(functionQE_real)
Z1_real = f1_real(X, Y)

f2_real = np.vectorize(function_real)
Z2_real = f2_real(X, Y)

def functionQE_imag(alpha_parallel,xD_tilde):
    rta = Efield_NUM_QE(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,b,zp,x0,y0,z0,alpha_parallel,xD_tilde)
    return rta.imag

def function_imag(alpha_parallel,xD_tilde):
    rta = Efield_NUM(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,b,zp,x0,y0,z0,alpha_parallel,xD_tilde)
    return rta.imag

f1_imag = np.vectorize(functionQE_imag)
Z1_imag = f1_imag(X, Y)

f2_imag = np.vectorize(function_imag)
Z2_imag = f2_imag(X, Y)

os.chdir(path_save)   
print('Graficar el External field' + label1)

graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
im = plt.imshow(Z1_real, extent = limits, cmap=plt.cm.hot, interpolation='bilinear')
cbar = plt.colorbar(im)
cbar.ax.tick_params(labelsize = tamnum)
if save_graphs == True:
    plt.tight_layout()
    plt.savefig('Re_function_Eref_QE' + label1 + '.png', format='png')


graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
im = plt.imshow(Z2_real, extent = limits, cmap=plt.cm.hot, interpolation='bilinear')
cbar = plt.colorbar(im)
cbar.ax.tick_params(labelsize = tamnum)
if save_graphs == True:
    plt.tight_layout()
    plt.savefig('Re_function_Eref' + label1 + '.png', format='png')
    

graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
im = plt.imshow(Z1_imag, extent = limits, cmap=plt.cm.hot, interpolation='bilinear')
cbar = plt.colorbar(im)
cbar.ax.tick_params(labelsize = tamnum)
if save_graphs == True:
    plt.tight_layout()
    plt.savefig('Imag_function_Eref_QE' + label1 + '.png', format='png')


graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
im = plt.imshow(Z2_imag, extent = limits, cmap=plt.cm.hot, interpolation='bilinear')
cbar = plt.colorbar(im)
cbar.ax.tick_params(labelsize = tamnum)
if save_graphs == True:
    plt.tight_layout()
    plt.savefig('Imag_function_Eref'  + label1 + '.png', format='png')    

#%%
#%%