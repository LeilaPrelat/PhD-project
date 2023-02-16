
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
path_constants =  path_basic.replace('/potential_field/infinite_dipoles','')
path_save = path_basic + '/' + 'decay_rate_theta_n'
if not os.path.exists(path_save):
    print('Creating folder to save graphs')
    os.mkdir(path_save)

err = 'decay_rate_theta_n.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from decay_rate_theta_n import decay_rate_theta_inf_dipoles_ana_res,decay_rate_theta_inf_dipoles_ana_res_div_gamma0_v3
except ModuleNotFoundError:
    print(err)

try:
    sys.path.insert(1, path_constants)
    from Silica_epsilon import epsilon_Silica
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


list_n = [0,1,2,3,4]   
print('Definir parametros del problema')

int_v = 10
Nmax = list_n[-1]
b = - 0.01


d_nano_film = 1

D_disk_nano = 100
d_thickness_disk_nano = 1


labely = r'$\Gamma_{\rm SP, n}/\Gamma_{\rm EELS}$'
#labely = r'Emission probability (eV$^{-1}$)'
labelp = r'_dfilm%.1fnm_ddisk%.1fnm_v%i' %(d_nano_film,d_thickness_disk_nano,int_v)
tabla = np.loadtxt('zp_optimum_for_decay_rate_hBN_disks_resonance' + labelp + '.txt', delimiter='\t', skiprows=1)
tabla = np.transpose(tabla)
[listx,listy,listz] = tabla

zp_nano = listy[-20]
zp_nano = 3.5
#zp_nano = 50
#zp_nano = listy[-20]
omegac0_1 = np.max(listx)/(c*hb)
lambda_SP_1 = 2*np.pi/omegac0_1

omegac0_2 = np.min(listx)/(c*hb)
lambda_SP_2 = 2*np.pi/omegac0_2


a_min = np.real(lambda_SP_1)*Nmax/(int_v - 1)
a_max = np.real(lambda_SP_2)*Nmax/(int_v + 1)

a = np.mean([a_min,a_max])
a = 0.001*1e-3
a = 120*1e-3
#a = 5*1e-3

#a = 150*1e-3
#
#a_nm = a*1e3


labelx = r'Surface-dipole distance, $z_{\rm 0}/\lambda_{\rm p}$'  
    
title2 = r'v = c/%i, b = %i nm' %(int_v,b*1e3) 
title3 = r'a = %i nm' %(a*1e3)
title4 = r', $z_p$ = $z^{opt}_p$' 

labelp = r'_a%.2fnm_zp%.2fnm_d%.2fnm' %(a*1e3,zp_nano,d_nano_film)
#label1 = 'vs_zp_lambda_p'  + labelp


#elif theta_degree == 60:
#    lim1, lim2 = 2.6,2.9
    
#def find_nearest(array, value):
#    array = np.asarray(array)
#    idx = (np.abs(array - value)).argmin()
#    return array[idx],idx    
#   
#num1,ind1 = find_nearest(np.array(listx), value=lim1)
#num2,ind2 = find_nearest(np.array(listx), value=lim2)

f1 = interp1d(listx, listy)
f2 = interp1d(listx, listz)

N = 225
lim1,lim2 = 18,-60
lim1,lim2 = 0,-58
#lim1,lim2 = 14,-1
listx_2 = np.linspace(listx[lim1], listx[lim2], N)
#listx_2 = np.linspace(listx[lim1], 0.2, N)
#
listy_2 = f1(listx_2)
listz_2 = f2(listx_2)  



title2 = r'v = c/%i, b = %i nm' %(int_v,b*1e3) 
title3 = r'a = %i nm' %(a*1e3)
title4 = r', $z_p$ = $z^{opt}_p$'


#labelp = r'_a%inm_zp_opt' %(a*1e3)
#label1 = 'vs_zp'  + labelp

    

#title1 = r'$\kappa$ = %.2f$\omega_o$, $\kappa_r$ = %.2f$\kappa$, $\hbar\omega_o$ = %i meV' %(kappa_factor_omega0, kappa_r_factor, energy0_pol)   
 
title =  title2 + '\n' +  title3  + title4



#%%


def function_real_ana(energy0_meV,zp_nano,n):
    
#    a = a_nano*1e-3
    omegac0 = energy0_meV/(c*hb)  
    zp = zp_nano*1e-3
         
    rta = decay_rate_theta_inf_dipoles_ana_res_div_gamma0_v3(omegac0,epsilon_Silica,d_nano_film,d_thickness_disk_nano,D_disk_nano,int_v,zp,a,b,n)

    return rta


tamfig = [2.5, 2]
tamletra = 7
tamtitle  = 8
tamnum = 6
tamlegend = 6
labelpady = 2
labelpadx = 3
pad = 2.5
mk = 1
ms = 1
hp = 0.5
length_marker = 1.5
dpi = 500


#%%
    


def graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad):
    plt.figure(figsize=tamfig)
    plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
    plt.ylabel(labely,fontsize=tamletra,labelpad =labelpadx)
    plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in", pad = pad)
#    plt.title(title,fontsize=int(tamtitle*0.9))

    return  
 
#%%

maxis = []

#    if theta_degree != 0:     
#        listx_2 = listx
#        listy_2 = listy
#        listz_2 = listz
for n in list_n:
    
    list_y_re = []


    for ind in range(len(listx_2)):
#        zp_nano = listy_2[ind]
        x =  listx_2[ind]
#        x = 43 #meV
        list_y_re.append(function_real_ana(x,zp_nano,n))
        
    maxi = np.max(list_y_re)
    maxis.append(maxi)
    print(n,listx_2[int(np.argmax(list_y_re))],maxi)
#    list_y_re = np.array(list_y_re)/maxi
     
    #%%
    
maxis = []
list_y_re_tot = []

for n in list_n:
    
    list_y_re = []


    for ind in range(len(listx_2)):
#        zp_nano = listy_2[ind]
        x =  listx_2[ind]
#        x = 43 #meV
        list_y_re.append(function_real_ana(x,zp_nano,n))
        
    maxi = np.max(list_y_re)
    maxis.append(maxi)
#    print(n,maxi)
#    list_y_re = np.array(list_y_re)/np.max(maxis)

#    list_y_re = np.array(list_y_re)*1e14
    
    list_y_re_tot.append(list_y_re)

    #%%

listx_3 = []
for ind in range(len(listy_2)):
    listx_3.append(listy_2[ind]/listz_2[ind])
    
listx_4 = np.linspace(np.min(listx_3),np.max(listx_3),N)

 #%% 
listx_3 = np.array(listy_2)/np.array(listz_2)
graph(title,labelx,labely ,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
k = 0
for n in list_n:
    
    plt.plot(listx_4,np.array(list_y_re_tot[k]),'-',ms = ms, label = 'n = %i'%(n))
    k = k + 1
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=False,handletextpad=hp, handlelength=1)
#    plt.grid(1)
plt.tight_layout()

#plt.xlim([0,200])
os.chdir(path_save)
plt.savefig('decay_rate_fix_zp_' + labelp + '.png', format='png',bbox_inches='tight',pad_inches = 0.008,dpi = dpi)


#%%

