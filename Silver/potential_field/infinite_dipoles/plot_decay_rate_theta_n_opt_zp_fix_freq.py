
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



plot_num = 0
plot_color_map = 0

plot_vs_theta = 0
plot_vs_zp = 0
plot_vs_zp_lambda_p = 1

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

d_nano = 1

int_v = 10

Nmax = 4

labely = r'$\Gamma_n$'
labely = r'Emission probability (arb. units)'


tabla = np.loadtxt('zp_optimum_for_decay_rate_resonance_d%inm_v%i.txt'%(d_nano,int_v), delimiter='\t', skiprows=1)
tabla = np.transpose(tabla)
[listx,listy,listz] = tabla

energy0_eV = listx[10]

omegac0_1 = np.max(listx)/(c*hb)
lambda_SP_1 = 2*np.pi/omegac0_1

omegac0_2 = np.min(listx)/(c*hb)
lambda_SP_2 = 2*np.pi/omegac0_2


a_min = np.real(lambda_SP_1)*Nmax/(int_v - 1)
a_max = np.real(lambda_SP_2)*Nmax/(int_v + 1)

a = np.mean([a_min,a_max])

a = 5031*1e-3
#a = 0.5

a_nm = a*1e3

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx],idx

#    listx_2 = zp/(2*np.pi/listx_3)
if plot_vs_zp == 1:    
    
    theta_degree = 30
    theta = theta_degree*np.pi/180

    labelx = r'$\hbar\omega$ [eV]'  
        
    title2 = r'v = c/%i, b = %i nm, $\theta$ = %.2fº' %(int_v,b*1e3,theta_degree) 
    title3 = r'a = %i nm' %(a*1e3)
    title4 = r', $z_p$ = $z^{opt}_p$, d = % i nm' %(d_nano)
    
    labelp = r'_a%inm_zp_opt_theta%i' %(a*1e3,theta_degree)
    label1 = 'vs_zp'  + labelp
    
elif plot_vs_zp_lambda_p == 1:        

    theta_degree = 45
    theta = theta_degree*np.pi/180

    labelx = r'Surface-dipole distance, $z_{\rm p}/\lambda_{\rm p}$'  
        
    title2 = r'v = c/%i, b = %i nm, $\theta$ = %.2fº' %(int_v,b*1e3,theta_degree) 
    title3 = r'a = %i nm' %(a*1e3)
    title4 = r', $z_p$ = $z^{opt}_p$, d = % i nm' %(d_nano)
    
    labelp = r'_a%inm_zp_opt_theta%i_d%inm' %(a*1e3,theta_degree,d_nano)
    label1 = 'vs_zp_lambda_p'  + labelp

    if theta_degree == 0:
        lim1, lim2 = 1.6,2.2
        
#    elif theta_degree == 60:
#        lim1, lim2 = 1.6,2.2
    elif theta_degree == 30:
        lim1, lim2 = 1.8,2.2


    elif theta_degree == 45:
        lim1, lim2 = 2.1,2.5

    elif theta_degree == 60:
        lim1, lim2 = 2.6,2.9
        
   
    num1,ind1 = find_nearest(np.array(listx), value=lim1)
    num2,ind2 = find_nearest(np.array(listx), value=lim2)
    
    f1 = interp1d(listx, listy)
    f2 = interp1d(listx, listz)
    
    listx_2 = np.linspace(listx[0], listx[-1], 500)
    
    
    listy_2 = f1(listx_2)
    listz_2 = f2(listx_2)
    
   
elif plot_vs_theta == 1:
    energy0 = 1
    labelx = r'$\theta$ [degree]'  
    
    title2 = r'v = c/%i, b = %i nm, $\hbar\omega$ = %i eV' %(int_v,b*1e3,energy0) 
    title3 = r'a = %i nm' %(a*1e3)
    title4 = r', $z_p$ = $z^{opt}_p$, d = % i nm' %(d_nano)
    
    labelp = r'_a%inm_zp_opt_theta%i' %(a*1e3,theta_degree)
    label1 = 'vs_zp'  + labelp 


    listx_2 = np.linspace(0,np.pi,len(listx))
    listx_3 = np.array(listx_2)*180/np.pi
    labelx = r'$\theta$ [degree]'  


title =  title2 + '\n' +  title3  + title4



#%%


def function_real_ana(zp_nano,Nmax):
                
    omegac0 = energy0_eV/(c*hb)  
    zp = zp_nano*1e-3
         
    rta = decay_rate_theta_inf_dipoles_ana_res(omegac0,epsi1,epsi3,d_nano,int_v,zp,a,b,Nmax)

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
    

if plot_vs_zp == 1:
    maxis = []
    list_n = [0,2,4]        
    list_y_re_tot = []
    for n in list_n:
        
        list_y_re = []
    
    
        for ind in range(len(listx)):
            zp = listy[ind]
            x =  listx[ind]
            list_y_re.append(function_real_ana(theta,zp,x,n))
            
        maxi = np.max(list_y_re)
        maxis.append(maxi)
        list_y_re_tot.append(list_y_re)
        print(n,maxi)
    #    list_y_re = np.array(list_y_re)/maxi
         
    maxis = []
        
    graph(title,labelx,labely ,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    for n in list_n:
        
        list_y_re = []
    
    
        for ind in range(len(listx)):
            zp = listy[ind]
            x =  listx[ind]
            list_y_re.append(function_real_ana(theta,zp,x,n))
            
        maxi = np.max(list_y_re)
        maxis.append(maxi)
        list_y_re_tot.append(list_y_re)
    #    print(n,maxi)
        list_y_re = np.array(list_y_re)/np.max(maxis)
        
        plt.plot(listx,list_y_re,'.-',ms = ms, label = 'n = %i'%(n))
        
    plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
    plt.grid(1)
    plt.tight_layout()
    #plt.yscale('log')
    os.chdir(path_save)
    plt.savefig('decay_rate_' + label1 + '.png', format='png')


elif plot_vs_zp_lambda_p == 1:
    maxis = []
    list_n = [3]   
#    if theta_degree != 0:     
#        listx_2 = listx
#        listy_2 = listy
#        listz_2 = listz
    for n in list_n:
        
        list_y_re = []
    
    
        for ind in range(len(listx_2)):
            zp = listy_2[ind]
            x =  listx_2[ind]
            list_y_re.append(function_real_ana(zp,n))
            
        maxi = np.max(list_y_re)
        maxis.append(maxi)
        print(n,maxi)
    #    list_y_re = np.array(list_y_re)/maxi
         
    maxis = []
        
    graph(title,labelx,labely ,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    for n in list_n:
        
        list_y_re = []
    
    
        for ind in range(len(listx_2)):
            zp = listy_2[ind]
            x =  listx_2[ind]
            list_y_re.append(function_real_ana(zp,n))
            
        maxi = np.max(list_y_re)
        maxis.append(maxi)
    #    print(n,maxi)
#        list_y_re = np.array(list_y_re)/np.max(maxis)
        
        listx_3 = np.array(listx_2)/np.array(listz_2)
        plt.plot(listx_3,list_y_re,'-',ms = ms, label = 'n = %i'%(n))
        
    plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=False,handletextpad=hp, handlelength=1)
#    plt.grid(1)
    plt.tight_layout()

    plt.yscale('log')
    os.chdir(path_save)
    plt.savefig('decay_rate_fix_freq' + label1 + '.png', format='png',dpi = dpi)


elif plot_vs_theta == 1:

    maxis = []
    list_n = [0,1,5]
    list_y_re_tot = []
    for n in list_n:
        
        list_y_re = []
    
    
        for ind in range(len(listx)):
            zp = listy[ind]
            x =  listx_3[ind]
            list_y_re.append(function_real_ana(x,zp,energy0,n))
            
        maxi = np.max(list_y_re)
        maxis.append(maxi)
        list_y_re_tot.append(list_y_re)
        print(n,maxi)
    #    list_y_re = np.array(list_y_re)/maxi
         
    maxis = []
        
    graph(title,labelx,labely ,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    for n in list_n:
        
        list_y_re = []
    
    
        for ind in range(len(listx)):
            zp = listy[ind]
            x =  listx_3[ind]
            list_y_re.append(function_real_ana(x,zp,energy0,n))
            
        maxi = np.max(list_y_re)
        maxis.append(maxi)
        list_y_re_tot.append(list_y_re)
    #    print(n,maxi)
        list_y_re = np.array(list_y_re)/np.max(maxis)
        
        plt.plot(listx_3,list_y_re,'.-',ms = ms, label = 'n = %i'%(n))
        
    plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
    plt.grid(1)
    plt.tight_layout()
    #plt.yscale('log')
    os.chdir(path_save)
    plt.savefig('decay_rate_' + label1 + '.png', format='png')    

#%%

