
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
    from decay_rate_theta_n import decay_rate_theta_inf_dipoles_ana_res,decay_rate_theta_inf_dipoles_ana_res_div_gamma0
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

epsi1, epsi3 = 1,1

print('Definir parametros del problema')


b = - 0.01

d_nano = 0.4

int_v = 10

Nmax = 4

labely = r'$\Gamma_{n,\rm SP}/\Gamma_{\rm 0}$ $\times$ $10^4$'
#labely = r'Emission probability (eV$^{-1}$)'

tabla = np.loadtxt('zp_optimum_for_decay_rate_hBN_disks_resonance_d%.2fnm_v%i.txt'%(d_nano,int_v), delimiter='\t', skiprows=1)
tabla = np.transpose(tabla)
[listx,listy,listz] = tabla

zp_nano = listy[0]
zp_nano = 3

#zp_nano = listy[-20]
omegac0_1 = np.max(listx)/(c*hb)
lambda_SP_1 = 2*np.pi/omegac0_1

omegac0_2 = np.min(listx)/(c*hb)
lambda_SP_2 = 2*np.pi/omegac0_2


a_min = np.real(lambda_SP_1)*Nmax/(int_v - 1)
a_max = np.real(lambda_SP_2)*Nmax/(int_v + 1)

a = np.mean([a_min,a_max])

#a = 8000*1e-3
a = 150*1e-3

a_nm = a*1e3


labelx = r'Surface-dipole distance, $z_{\rm 0}/\lambda_{\rm p}$'  
    
title2 = r'v = c/%i, b = %i nm' %(int_v,b*1e3) 
title3 = r'a = %i nm' %(a*1e3)
title4 = r', $z_p$ = $z^{opt}_p$' 

labelp = r'_a%inm_zp%inm_d%inm' %(a*1e3,zp_nano,d_nano)
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

N = 500
lim1,lim2 = 18,-60
lim1,lim2 = 0,-1
listx_2 = np.linspace(listx[lim1], listx[lim2], N)


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


def function_real_ana(energy0_meV,zp_nano,Nmax):
    
#    a = a_nano*1e-3
    omegac0 = energy0_meV/(c*hb)  
    zp = zp_nano*1e-3
         
    rta = decay_rate_theta_inf_dipoles_ana_res_div_gamma0(omegac0,epsilon_Silica,d_nano,int_v,zp,a,b,Nmax)

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
list_n = [0,1,2,3,4]   
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
    print(n,maxi)
#    list_y_re = np.array(list_y_re)/maxi
     
maxis = []
    
graph(title,labelx,labely ,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
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
    
    listx_3 = np.array(listx_2)/np.array(listz_2)
    plt.plot(listx_3,np.array(list_y_re)*1e-4,'-',ms = ms, label = 'n = %i'%(n))
    
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=False,handletextpad=hp, handlelength=1)
#    plt.grid(1)
plt.tight_layout()

#plt.yscale('log')
os.chdir(path_save)
plt.savefig('decay_rate_fix_zp_' + labelp + '.png', format='png',bbox_inches='tight',pad_inches = 0.008,dpi = dpi)


#%%

