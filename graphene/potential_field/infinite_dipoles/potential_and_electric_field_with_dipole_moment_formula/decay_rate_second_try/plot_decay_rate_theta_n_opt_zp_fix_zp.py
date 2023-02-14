
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
    from decay_rate_theta_n import decay_rate_theta_inf_dipoles_ana_res,decay_rate_theta_inf_dipoles_ana_res_div_gamma0_v3
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
list_n = [0,1,2,3,4] 

int_v = 20
int_v = 10

Nmax = list_n[-1]

labely = r'$\Gamma_{\rm SP, n}/\Gamma_{\rm EELS}$ $\times$ $10^8$'
labely = 'Decay rate of surface' + '\n' + r'plasmons $\Gamma_{n,\rm SP}/\Gamma_{\rm 0}$ $\times$ $10^8$'
labely = r'$\Gamma_{\rm SP, n}/\Gamma_{\rm EELS}$'
#labely = r'Emission probability (eV$^{-1}$)'

tabla = np.loadtxt('zp_optimum_for_decay_rate_graphene_resonance_b-10nm.txt', delimiter='\t', skiprows=1)
tabla = np.transpose(tabla)
[listx,listy,listz] = tabla ## energy and zp 

ind = 18
ind = -20
print(listx[ind])
zp_nano = listy[ind]

omegac0_1 = np.max(listx)*1e-3/(c*hb)
lambda_SP_1 = 2*np.pi/omegac0_1

omegac0_2 = np.min(listx)*1e-3/(c*hb)
lambda_SP_2 = 2*np.pi/omegac0_2


a_min = np.real(lambda_SP_1)*Nmax/(int_v - 1)
a_max = np.real(lambda_SP_2)*Nmax/(int_v + 1)

a = np.mean([a_min,a_max]) ##  micrones 

a = 1.5
a = 30
#a = 500*1e-3
#a = 5031*1e-3
print(a)


a_nm = a*1e3


labelx = r'Surface-dipole distance, $z_{\rm 0}/\lambda_{\rm p}$'  
    
title2 = r'v = c/%i, b = %i nm' %(int_v,b*1e3) 
title3 = r'a = %i nm' %(a*1e3)
title4 = r', $z_p$ = $z^{opt}_p$' 

labelp = r'_a%inm_zp%inm' %(a*1e3,zp_nano)
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

N = 5000
lim1,lim2 = 16,-80
if a == 15: 
    lim1,lim2 = 16,-80
#lim1,lim2 = 10,-63
    listx_2 = np.linspace(listx[lim1], listx[lim2], N)
    listx_2 = np.linspace(38,70,N)
elif a == 150:
    lim1,lim2 = 16,-75
    listx_2 = np.linspace(listx[lim1], listx[lim2], N)
else:
    lim1,lim2 = 16,-70
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


def function_real_ana(energy0_meV,Nmax):
                
    omegac0 = energy0_meV*1e-3/(c*hb)  
    zp = zp_nano*1e-3
         
    rta = decay_rate_theta_inf_dipoles_ana_res_div_gamma0_v3(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,zp,a,b,Nmax)

    return rta


tamfig = [2.5, 2]
tamletra = 9
tamtitle  = 9
tamnum = 7
tamlegend = 7
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
    plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
    plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in", pad = pad)
#    plt.tick_params(labelsize = tamnum, length = 4 , width=1, direction="in", pad = pad) ### para cleo europe
#    plt.title(title,fontsize=int(tamtitle*0.9))

    return  
 
#%%
from scipy.signal import find_peaks
maxis_zp_lambda_p = []
  
#    if theta_degree != 0:     
#        listx_2 = listx
#        listy_2 = listy
#        listz_2 = listz
for n in list_n:
    
    list_y_re = []


    for ind in range(len(listx_2)):
#        zp = listy_2[ind]
        x =  listx_2[ind]
#        x = 43 #meV
        list_y_re.append(function_real_ana(x,n))
        
    maxi, _ = find_peaks(list_y_re, height=0)
    if len(maxi)>1:
        list_maximos = []
        for maxis in maxi: 
            list_maximos.append(list_y_re[maxis])
        index_max = np.argmax(list_maximos)
        maxis_zp_lambda_p.append(listx_2[maxi[int(index_max)]])
    else:
        print(int(maxi))
        maxis_zp_lambda_p.append(listx_2[int(maxi)])
    print(n,maxis_zp_lambda_p)
#    list_y_re = np.array(list_y_re)/maxi
   
    
#%%
listx_3 = []
for ind in range(len(listy_2)):
    listx_3.append(listy_2[ind]/listz_2[ind])
    
listx_4 = np.linspace(np.min(listx_3),np.max(listx_3),N)

    
graph(title,'Energy [meV]',labely ,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
for n in list_n:
    
    list_y_re = []


    for ind in range(len(listx_2)):
#        zp = listy_2[ind]
        x =  listx_2[ind]
#        x = 43 #meV
        list_y_re.append(function_real_ana(x,n))
            
    plt.plot(listx_2,np.array(list_y_re)*1e2,'-',ms = ms, label = 'n = %i'%(n))
    
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=False,handletextpad=hp, handlelength=1)
#    plt.grid(1)
plt.tight_layout()
    
graph(title,labelx,labely ,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
for n in list_n:
    
    list_y_re = []


    for ind in range(len(listx_2)):
#        zp = listy_2[ind]
        x =  listx_2[ind]
#        x = 43 #meV
        list_y_re.append(function_real_ana(x,n))
            
    listx_3 = np.array(listy_2)/np.array(listz_2)
    plt.plot(listx_4,np.array(list_y_re),'-',ms = ms, label = 'n = %i'%(n))
    
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=False,handletextpad=hp, handlelength=1)
#    plt.grid(1)
plt.tight_layout()

#plt.yscale('log')
os.chdir(path_save)
#plt.savefig('decay_rate_fix_zp_' + labelp + '.png', format='png',bbox_inches='tight',pad_inches = 0.008,dpi = dpi)
plt.savefig('decay_rate_fix_zp_cleo_' + labelp + '.png', format='png',bbox_inches='tight',pad_inches = 0.008,dpi = dpi)

#%%

