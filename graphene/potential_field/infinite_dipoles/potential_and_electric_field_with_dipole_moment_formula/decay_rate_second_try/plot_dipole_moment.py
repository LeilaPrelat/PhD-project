
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
import seaborn as sns
sns.set()

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/potential_field/potential/potential_and_electric_field_with_dipole_moment_formula','')
path_save = path_basic + '/' + 'dipole_moment_opt_values'
if not os.path.exists(path_save):
    print('Creating folder to save graphs')
    os.mkdir(path_save)

err = 'fieldE_direct_numerical.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from dipole_moment import dipole_moment_anav2_res
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

int_v = 10

epsi1,epsi2 = 1,1
hbmu,hbgama = 0.3,0.0001

b = -0.01

tabla = np.loadtxt('zp_optimum_for_decay_rate_graphene_resonance_b-10nm.txt', delimiter='\t', skiprows=1)
tabla = np.transpose(tabla)
[listx,listy,listz] = tabla ## energy and zp 

f1 = interp1d(listx, listy)
f2 = interp1d(listx, listz)

N = 5000

listx_2 = np.linspace(listx[0], listx[-1], N)


listy_2 = f1(listx_2)
listz_2 = f2(listx_2)  


#title2 = r'$\hbar\mu$ = %.2feV, $\hbar\gamma$ = %.4feV' %(hbmu,hbgama) 
#title3 = r'$z_p$=%inm, px=%i, py=%i, pz=%i' %(zp*1e3,px,py,pz)
title4 = r'v = c/%i b = %i nm' %(int_v,b*1e3)
labelp = r'_res_opt_values' 

N = 100

    # z0 = 0.06*1e3
labelx = r'Surface-dipole distance, $z_{\rm p}/\lambda_{\rm p}$' 
label1 = 'vs_zp_lambda_p'  + labelp


title =  title4 

def function_ana_res_x(energy0_meV,zp_nano):
    omegac0 = energy0_meV*1e-3/aux 
    zp = zp_nano*1e-3
    
    px_f,py_f,pz_f = dipole_moment_anav2_res(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,b,zp)
    # print(rta)
#    rta = np.abs(px_f)**2 + np.abs(py_f)**2 + np.abs(pz_f)**2 

    return px_f



def function_ana_res_y(energy0_meV,zp_nano):
    omegac0 = energy0_meV*1e-3/aux 
    zp = zp_nano*1e-3
    
    px_f,py_f,pz_f = dipole_moment_anav2_res(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,b,zp)
    # print(rta)
#    rta = np.abs(px_f)**2 + np.abs(py_f)**2 + np.abs(pz_f)**2 

    return py_f



def function_ana_res_z(energy0_meV,zp_nano):
    omegac0 = energy0_meV*1e-3/aux 
    zp = zp_nano*1e-3

    px_f,py_f,pz_f = dipole_moment_anav2_res(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,b,zp)
    # print(rta)
#    rta = np.abs(px_f)**2 + np.abs(py_f)**2 + np.abs(pz_f)**2 

    return pz_f

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

#%%

listy_re_ana_x = []
listy_im_ana_x = []
listy_abs_ana_x = []

listy_re_ana_y = []
listy_im_ana_y = []
listy_abs_ana_y = []

listy_re_ana_z = []
listy_im_ana_z = []
listy_abs_ana_z = []



for ind in range(len(listx_2)):
    
    zp_nano = listy_2[ind]
    energy_meV =  listx_2[ind]

   
    px = function_ana_res_x(energy_meV,zp_nano)  
    py = function_ana_res_y(energy_meV,zp_nano)  
    pz = function_ana_res_z(energy_meV,zp_nano)  
    
    
    listy_re_ana_x.append(np.real(px))
    listy_im_ana_x.append(np.imag(px))
    listy_abs_ana_x.append(np.abs(px))


    listy_re_ana_y.append(np.real(py))
    listy_im_ana_y.append(np.imag(py))
    listy_abs_ana_y.append(np.abs(py))


    listy_re_ana_z.append(np.real(pz))
    listy_im_ana_z.append(np.imag(pz))
    listy_abs_ana_z.append(np.abs(pz))



#%%
    
listx_3 = np.array(listx_2)/np.array(listz_2)

label_extra = r'/$\alpha_{eff}$'

#%%

graph(title,labelx,r'Re($p_x$)' ,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#    plt.plot(listx,listy_re_ana,'.',ms = ms,color = 'purple',label = 'analytical')
plt.plot(listx_3,listy_re_ana_x,'.',ms = ms,color = 'purple',label = 'PP analytical')
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
#plt.yscale('log')
os.chdir(path_save)
plt.savefig( 're_px' + label1 + '.png', format='png')   

#%%

graph(title,labelx,r'Re($p_y$)' ,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#    plt.plot(listx,listy_re_ana,'.',ms = ms,color = 'purple',label = 'analytical')
plt.plot(listx_3,listy_re_ana_y,'.',ms = ms,color = 'purple',label = 'PP analytical')
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
#plt.yscale('log')
os.chdir(path_save)
plt.savefig( 're_py' + label1 + '.png', format='png')   

#%%

graph(title,labelx,r'Re($p_z$)' ,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#    plt.plot(listx,listy_re_ana,'.',ms = ms,color = 'purple',label = 'analytical')
plt.plot(listx_3,listy_re_ana_z,'.',ms = ms,color = 'purple',label = 'PP analytical')
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
#plt.yscale('log')
os.chdir(path_save)
plt.savefig( 're_pz' + label1 + '.png', format='png')   

#%%

graph(title,labelx,r'Im($p_x$)' ,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#    plt.plot(listx,listy_re_ana,'.',ms = ms,color = 'purple',label = 'analytical')
plt.plot(listx_3,listy_im_ana_x,'.',ms = ms,color = 'purple',label = 'PP analytical')
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
#plt.yscale('log')
os.chdir(path_save)
plt.savefig( 'im_px' + label1 + '.png', format='png')   

#%%

graph(title,labelx,r'Im($p_y$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#    plt.plot(listx,listy_re_ana,'.',ms = ms,color = 'purple',label = 'analytical')
plt.plot(listx_3,listy_im_ana_y,'.',ms = ms,color = 'purple',label = 'PP analytical')
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
#plt.yscale('log')
os.chdir(path_save)
plt.savefig( 'im_py' + label1 + '.png', format='png')   

#%%

graph(title,labelx,r'Im($p_z$)' ,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#    plt.plot(listx,listy_re_ana,'.',ms = ms,color = 'purple',label = 'analytical')
plt.plot(listx_3,listy_im_ana_z,'.',ms = ms,color = 'purple',label = 'PP analytical')
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
#plt.yscale('log')
os.chdir(path_save)
plt.savefig( 'im_pz' + label1 + '.png', format='png')   

#%%


graph(title,labelx,r'|$p_x$|' ,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#    plt.plot(listx,listy_re_ana,'.',ms = ms,color = 'purple',label = 'analytical')
plt.plot(listx_3,listy_im_ana_x,'.',ms = ms,color = 'purple',label = 'PP analytical')
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
#plt.yscale('log')
os.chdir(path_save)
plt.savefig( 'abs_px' + label1 + '.png', format='png')   

#%%

graph(title,labelx,r'|$p_y$|',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#    plt.plot(listx,listy_re_ana,'.',ms = ms,color = 'purple',label = 'analytical')
plt.plot(listx_3,listy_im_ana_y,'.',ms = ms,color = 'purple',label = 'PP analytical')
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
#plt.yscale('log')
os.chdir(path_save)
plt.savefig( 'abs_py' + label1 + '.png', format='png')   

#%%

graph(title,labelx,r'|$p_z$|' ,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#    plt.plot(listx,listy_re_ana,'.',ms = ms,color = 'purple',label = 'analytical')
plt.plot(listx_3,listy_im_ana_z,'.',ms = ms,color = 'purple',label = 'PP analytical')
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
#plt.yscale('log')
os.chdir(path_save)
plt.savefig( 'abs_pz' + label1 + '.png', format='png')   






















