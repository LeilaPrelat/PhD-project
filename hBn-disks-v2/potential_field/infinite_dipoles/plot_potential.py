
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
#import seaborn as sns
#sns.set()
from scipy.interpolate import interp1d

plot_vs_zp = 0
plot_vs_zp_lambda_p = 1

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/potential_field/infinite_dipoles','')
path_save = path_basic + '/' + 'potential'
if not os.path.exists(path_save):
    print('Creating folder to save graphs')
    os.mkdir(path_save)

err = 'decay_rate_theta_n.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from potential_n import potential_inf_dipoles_pole_aprox, potential_inf_dipoles_num, potential_inf_dipoles_ana
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


list_n = [0,1,2,3,4]

labely = r'$\phi_n$'

int_v = 10
Nmax = 4
b = - 0.01

d_nano_film = 1

D_disk_nano = 100
d_thickness_disk_nano = 1


#labely = r'Emission probability (eV$^{-1}$)'
labelp = r'_dfilm%.1fnm_ddisk%.1fnm_v%i' %(d_nano_film,d_thickness_disk_nano,int_v)
tabla = np.loadtxt('zp_optimum_for_decay_rate_hBN_disks_resonance' + labelp + '.txt', delimiter='\t', skiprows=1)
tabla = np.transpose(tabla)
[listx,listy,listz] = tabla

zp_nano = listy[0]
zp_nano = 0.05

#zp_nano = listy[-20]
omegac0_1 = np.max(listx)/(c*hb)
lambda_SP_1 = 2*np.pi/omegac0_1

omegac0_2 = np.min(listx)/(c*hb)
lambda_SP_2 = 2*np.pi/omegac0_2


a_min = np.real(lambda_SP_1)*Nmax/(int_v - 1)
a_max = np.real(lambda_SP_2)*Nmax/(int_v + 1)

a = np.mean([a_min,a_max])

#a = 8000*1e-3
a = 185*1e-3

a_nm = a*1e3


labelx = r'Surface-dipole distance, $z_{\rm 0}/\lambda_{\rm p}$'  
    
title2 = r'v = c/%i, b = %i nm' %(int_v,b*1e3) 
title3 = r'a = %i nm' %(a*1e3)
title4 = r', $z_p$ = $z^{opt}_p$' 

labelp = r'_a%.2fnm_zp%.2fnm_d%inm' %(a*1e3,zp_nano,d_nano_film)
#label1 = 'vs_zp_lambda_p'  + labelp


f1 = interp1d(listx, listy)
f2 = interp1d(listx, listz)

N = 200
lim1,lim2 = 18,-60
lim1,lim2 = 0,-1
listx_2 = np.linspace(listx[lim1], listx[lim2], N)


listy_2 = f1(listx_2)
listz_2 = f2(listx_2)  



title2 = r'v = c/%i, b = %i nm' %(int_v,b*1e3) 
title3 = r'a = %i nm' %(a*1e3)
title4 = r', $z_p$ = $z^{opt}_p$'

 
title =  title2 + '\n' +  title3  + title4


#%%

def function_ana(zp_nano,energy0_eV,Nmax):
                
    omegac0 = energy0_eV/(c*hb)  
    zp = zp_nano*1e-3
         
    rta = potential_inf_dipoles_ana(omegac0,epsilon_Silica,d_nano_film,d_thickness_disk_nano,D_disk_nano,int_v,zp,a,b,Nmax)

    return rta


def function_num(zp_nano,energy0_eV,Nmax):
                
    omegac0 = energy0_eV/(c*hb)  
    zp = zp_nano*1e-3
         
    rta = potential_inf_dipoles_num(omegac0,epsilon_Silica,d_nano_film,d_thickness_disk_nano,D_disk_nano,int_v,zp,a,b,Nmax)

    return rta


def function_pole_aprox(zp_nano,energy0_eV,Nmax):
                
    omegac0 = energy0_eV/(c*hb)  
    zp = zp_nano*1e-3
         
    rta = potential_inf_dipoles_pole_aprox(omegac0,epsilon_Silica,d_nano_film,d_thickness_disk_nano,D_disk_nano,int_v,zp,a,b,Nmax)

    return rta


#%%

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
    
list_y_re_tot = []
list_y_im_tot = []


listx_3 = []
for ind in range(len(listy_2)):
    listx_3.append(listy_2[ind]/listz_2[ind])
    
listx_4 = np.linspace(np.min(listx_3),np.max(listx_3),N)

if plot_vs_zp == 1:
    
    for n in list_n:
        print(n)
        
        list_y_re = []
        list_y_im = []
    
        for ind in range(len(listx)):
#            zp = listy[ind]
            x =  listx[ind]
            rta = function_ana(zp_nano,x,n)
            list_y_re.append(np.real(rta))
            list_y_im.append(np.imag(rta))
            

        list_y_re_tot.append(list_y_re)
        list_y_im_tot.append(list_y_im)
    
    
    graph(title,labelx,labely ,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    for n in list_n:
        plt.plot(listx_4,list_y_re_tot[n],'.-',ms = ms, label = 'n = %i'%(n))  
    plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=False,handletextpad=hp, handlelength=1)
#    plt.grid(1)
    plt.tight_layout()
    #plt.yscale('log')
    os.chdir(path_save)
    plt.savefig('re_phi_n_' + labelp + '.png', format='png',bbox_inches='tight',pad_inches = 0.008,dpi = dpi)



    graph(title,labelx,labely ,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    for n in list_n:
        plt.plot(listx_4,list_y_im_tot[n],'.-',ms = ms, label = 'n = %i'%(n))  
    plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=False,handletextpad=hp, handlelength=1)
#    plt.grid(1)
    plt.tight_layout()
    #plt.yscale('log')
    os.chdir(path_save)
    plt.savefig('im_phi_n_' + labelp + '.png', format='png',bbox_inches='tight',pad_inches = 0.008,dpi = dpi)

    for n in list_n:
        tabla = np.array([listx_4,list_y_re_tot[n],list_y_im_tot[n]])
        tabla = np.transpose(tabla)
        header1 = 'E [eV]     Re(phi,n)' + ', ' + title + ', ' + name_this_py
        np.savetxt( 'phi_inf_dip_hBN_disks_ana_' + labelp + 'n%i'%(n) + '.txt', tabla, fmt='%1.11e', delimiter='\t', header = header1)
    
    
elif plot_vs_zp_lambda_p == 1:
        
    for n in list_n:
        print(n)
        
        list_y_re = []
        list_y_im = []
    
        for ind in range(len(listx_2)):
#            zp = listy_2[ind]
            x =  listx_2[ind]  ## energy 
            rta = function_ana(zp_nano,x,n)
            list_y_re.append(np.real(rta))
            list_y_im.append(np.imag(rta))
            

        list_y_re_tot.append(list_y_re)
        list_y_im_tot.append(list_y_im)
    
    
    graph(title,labelx,r'Re{$\phi$}' ,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    for n in list_n:
        plt.plot(listx_4,list_y_re_tot[n],'.-',ms = ms, label = 'n = %i'%(n))  
    plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=False,handletextpad=hp, handlelength=1)
#    plt.grid(1)
    plt.tight_layout()
    #plt.yscale('log')
    os.chdir(path_save)
    plt.savefig('re_phi_n_' + labelp + '.png', format='png',bbox_inches='tight',pad_inches = 0.008,dpi = dpi)


    graph(title,labelx,r'Im{$\phi$}' ,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    for n in list_n:
        plt.plot(listx_4,list_y_im_tot[n],'.-',ms = ms, label = 'n = %i'%(n))  
    plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=False,handletextpad=hp, handlelength=1)
#    plt.grid(1)
    plt.tight_layout()
    #plt.yscale('log')
    os.chdir(path_save)
    plt.savefig('im_phi_n_' + labelp + '.png', format='png',bbox_inches='tight',pad_inches = 0.008,dpi = dpi)


    for n in list_n:
        tabla = np.array([listx_4,list_y_re_tot[n],list_y_im_tot[n]])
        tabla = np.transpose(tabla)
        header1 = 'E [eV]     Re(phi,n)' + ', ' + title + ', ' + name_this_py
        np.savetxt( 'phi_inf_dip_hBN_disks_ana_' + labelp + 'n%i'%(n) + '.txt', tabla, fmt='%1.11e', delimiter='\t', header = header1)

#%%

