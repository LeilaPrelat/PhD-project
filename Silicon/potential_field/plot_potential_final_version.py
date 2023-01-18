
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
import seaborn as sns
sns.set()

#%%

plot_vs_R = 0
plot_vs_E = 1

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_ctes =  path_basic.replace('/' + 'potential_field','')
path_save = path_basic + '/' + 'potential_final_version'
if not os.path.exists(path_save):
    print('Creating folder to save graphs')
    os.mkdir(path_save)

err = 'fieldE_direct_numerical.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from potential_final_version import G1_num, G1_pole_aprox, G2_num, G2_pole_aprox
except ModuleNotFoundError:
    print(err)

try:
    sys.path.insert(1, path_ctes)
    from constants import constantes
except ModuleNotFoundError:
    print('constants.py no se encuentra en ' + path_ctes)

pi,hb,c,alfac,mu1,mu2,mu3 = constantes()
aux = c*hb
aux2 = 1e12/c

#%%

print('Definir parametros del problema')

#v = c/int_v
#omega = 0.7*1e12
cota = 75 #nanometros

epsilon2 = 12
d = 20*1e-3
zp = 0.05
z0 = zp

title = r'z = %inm, $z_p$=%inm, $\epsilon_2$ = %i, d = %inm' %(z0*1e3,zp*1e3,epsilon2,d*1e3)
labelp = r'_d%inm' %(d*1e3)

N = 10
if plot_vs_R == 1: 
    omega_omegaWG0 = 2 # meV 
    labelx = 'R [nm]'  
    title = title + ', ' + r'$\omega/\omega_{WG}$ = %i' %(omega_omegaWG0)
    label1 = 'vs_R' + labelp
    listx = np.linspace(100,4000,N)
else:
    R0 = 5000 #nanos
    R0 = 1000 #nanos
    # z0 = 0.06*1e3
    labelx = '$\omega/\omega_{WG}$'   
    title = title + ', ' + 'R = %inm' %(R0)
    label1 = 'vs_omega_omega_WG' + labelp
    listx = np.linspace(0.0001,2,N)
    listx = np.linspace(0.1,4,N)

omegaWG = (np.pi/d)/(np.sqrt(epsilon2-1)) # el c esta metido en aux

#%%

#def G1_real_ana(R_nano,energy0):
#    omegac0 = energy0*1e-3/aux 
#    R = R_nano*1e-3 
#
#    rta = G1_ana(omegac0,epsi1,epsi2,hbmu,hbgama,phi0,R,z0,zp,px,py,pz)
#    
#    return rta.real
#
#def G1_imag_ana(R_nano,energy0):
#    omegac0 = energy0*1e-3/aux 
#    R = R_nano*1e-3 
#    
#    rta = G1_ana(omegac0,epsi1,epsi2,hbmu,hbgama,phi0,R,z0,zp,px,py,pz)
#    return rta.imag
#
#def G2_real_ana(R_nano,energy0):
#    omegac0 = energy0*1e-3/aux 
#    R = R_nano*1e-3 
#
#    rta = G2_ana(omegac0,epsi1,epsi2,hbmu,hbgama,phi0,R,z0,zp,px,py,pz)
#    
#    return rta.real
#
#def G2_imag_ana(R_nano,energy0):
#    omegac0 = energy0*1e-3/aux 
#    R = R_nano*1e-3 
#    
#    rta = G2_ana(omegac0,epsi1,epsi2,hbmu,hbgama,phi0,R,z0,zp,px,py,pz)
#    return rta.imag

#%%

def G1_real_num(R_nano,omega_omegaWG):
#    omega0 = energy0*1e-3/aux 
#    omega_omegaWG = omega0/omegaWG
    R = R_nano*1e-3

    rta = G1_num(omega_omegaWG,d,R,z0,zp)
    
    return rta.real

def G1_imag_num(R_nano,omega_omegaWG):
#    omega0 = energy0*1e-3/aux    
#    omega_omegaWG = omega0/omegaWG
    R = R_nano*1e-3
    
    rta = G1_num(omega_omegaWG,d,R,z0,zp)
    return rta.imag

def G2_real_num(R_nano,omega_omegaWG):
#    omega0 = energy0*1e-3/aux    
#    omega_omegaWG = omega0/omegaWG
    R = R_nano*1e-3

    rta = G2_num(omega_omegaWG,d,R,z0,zp)
    
    return rta.real

def G2_imag_num(R_nano,omega_omegaWG):
#    omega0 = energy0*1e-3/aux    
#    omega_omegaWG = omega0/omegaWG
    R = R_nano*1e-3
    
    rta = G2_num(omega_omegaWG,d,R,z0,zp)
    return rta.imag

#%%


def G1_real_pole_aprox(R_nano,omega_omegaWG):
#    omega0 = energy0*1e-3/aux 
#    omega_omegaWG = omega0/omegaWG
    R = R_nano*1e-3

    rta = G1_pole_aprox(omega_omegaWG,d,R,z0,zp)
    
    return rta.real

def G1_imag_pole_aprox(R_nano,omega_omegaWG):
#    omega0 = energy0*1e-3/aux 
#    omega_omegaWG = omega0/omegaWG
    R = R_nano*1e-3
    
    rta = G1_pole_aprox(omega_omegaWG,d,R,z0,zp)
    return rta.imag


def G2_real_pole_aprox(R_nano,omega_omegaWG):
#    omega0 = energy0*1e-3/aux 
#    omega_omegaWG = omega0/omegaWG
    R = R_nano*1e-3

    rta = G2_pole_aprox(omega_omegaWG,d,R,z0,zp)
    
    return rta.real

def G2_imag_pole_aprox(R_nano,omega_omegaWG):
#    omega0 = energy0*1e-3/aux 
#    omega_omegaWG = omega0/omegaWG
    R = R_nano*1e-3
    
    rta = G2_pole_aprox(omega_omegaWG,d,R,z0,zp)
    return rta.imag

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

if plot_vs_R == 1: 
    
    listy_re_anaG1 = []
    listy_im_anaG1 = []

    listy_re_numG1 = []
    listy_im_numG1 = []
   
    listy_re_pole_aproxG1 = []
    listy_im_pole_aproxG1 = []
   
    listy_re_anaG2 = []
    listy_im_anaG2 = []

    listy_re_numG2 = []
    listy_im_numG2 = []
   
    listy_re_pole_aproxG2 = []
    listy_im_pole_aproxG2 = []
       
    
    for value in listx: 

#        y_re_anaG1 = G1_real_ana(value,E0)
#        y_im_anaG1 = G1_imag_ana(value,E0)        


#        y_re_numG1 = G1_real_num(value,omega_omegaWG0)
#        y_im_numG1 = G1_imag_num(value,omega_omegaWG0)  
        
        
        y_re_pole_aproxG1 = G1_real_pole_aprox(value,omega_omegaWG0)
        y_im_pole_aproxG1 = G1_imag_pole_aprox(value,omega_omegaWG0)  
        
        
#        listy_re_anaG1.append(y_re_anaG1)
#        listy_im_anaG1.append(y_im_anaG1)
        
#        listy_re_numG1.append(y_re_numG1)
#        listy_im_numG1.append(y_im_numG1)   
        
        listy_re_pole_aproxG1.append(y_re_pole_aproxG1)
        listy_im_pole_aproxG1.append(y_im_pole_aproxG1)   


#        y_re_anaG2 = G2_real_ana(value,E0)
#        y_im_anaG2 = G2_imag_ana(value,E0)        


#        y_re_numG2 = G2_real_num(value,omega_omegaWG0)
#        y_im_numG2 = G2_imag_num(value,omega_omegaWG0)  
        
        
        y_re_pole_aproxG2 = G2_real_pole_aprox(value,omega_omegaWG0)
        y_im_pole_aproxG2 = G2_imag_pole_aprox(value,omega_omegaWG0)  
        
        
#        listy_re_anaG2.append(y_re_anaG2)
#        listy_im_anaG2.append(y_im_anaG2)
        
#        listy_re_numG2.append(y_re_numG2)
#        listy_im_numG2.append(y_im_numG2)   
        
        listy_re_pole_aproxG2.append(y_re_pole_aproxG2)
        listy_im_pole_aproxG2.append(y_im_pole_aproxG2)   
        
    
    graph(title,labelx,'Re($G_1$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#    plt.plot(listx,listy_re_anaG1,'.',ms = ms,color = 'purple',label = 'PP analytical')
    plt.plot(listx,listy_re_pole_aproxG1,'.',ms = ms,color = 'darkred',label = 'PP numerical')
#    plt.plot(listx,listy_re_numG1,'.-',ms = ms,color = 'lightseagreen',label = 'full numerical')
    plt.legend(loc = 'lower left',markerscale=mk,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'Re_G1' + label1 + '.png', format='png')   

    graph(title,labelx,'Im($G_1$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#    plt.plot(listx,listy_im_anaG1,'.',ms = ms,color = 'purple',label = 'PP analytical')
    plt.plot(listx,listy_im_pole_aproxG1,'.',ms = ms,color = 'darkred',label = 'PP numerical')
#    plt.plot(listx,listy_im_numG1,'.-',ms = ms,color = 'lightseagreen',label = 'full numerical')
    plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'Im_G1' + label1 + '.png', format='png')   


    
    graph(title,labelx,'Re($G_2$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#    plt.plot(listx,listy_re_anaG2,'.',ms = ms,color = 'purple',label = 'PP analytical')
    plt.plot(listx,listy_re_pole_aproxG2,'.',ms = ms,color = 'darkred',label = 'PP numerical')
#    plt.plot(listx,listy_re_numG2,'.-',ms = ms,color = 'lightseagreen',label = 'full numerical')
    plt.legend(loc = 'lower left',markerscale=mk,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'Re_G2' + label1 + '.png', format='png')   

    graph(title,labelx,'Im($G_2$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#    plt.plot(listx,listy_im_anaG2,'.',ms = ms,color = 'purple',label = 'PP analytical')
    plt.plot(listx,listy_im_pole_aproxG2,'.',ms = ms,color = 'darkred',label = 'PP numerical')
#    plt.plot(listx,listy_im_numG2,'.-',ms = ms,color = 'lightseagreen',label = 'full numerical')
    plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'Im_G2' + label1 + '.png', format='png')   
    

if plot_vs_E == 1: 
    listy_re_anaG1 = []
    listy_im_anaG1 = []

    listy_re_numG1 = []
    listy_im_numG1 = []
   
    listy_re_pole_aproxG1 = []
    listy_im_pole_aproxG1 = []
   
    listy_re_anaG2 = []
    listy_im_anaG2 = []

    listy_re_numG2 = []
    listy_im_numG2 = []
   
    listy_re_pole_aproxG2 = []
    listy_im_pole_aproxG2 = []
   
    
    j = 0
    for value in listx: 

#        y_re_anaG1 = G1_real_ana(R0,value)
#        y_im_anaG1 = G1_imag_ana(R0,value)        


        y_re_numG1 = G1_real_num(R0,value)
        y_im_numG1 = G1_imag_num(R0,value)  
        
        
        y_re_pole_aproxG1 = G1_real_pole_aprox(R0,value)
        y_im_pole_aproxG1 = G1_imag_pole_aprox(R0,value)  
        
        
#        listy_re_anaG1.append(y_re_anaG1)
#        listy_im_anaG1.append(y_im_anaG1)
        
        listy_re_numG1.append(y_re_numG1)
        listy_im_numG1.append(y_im_numG1)   
        
        listy_re_pole_aproxG1.append(y_re_pole_aproxG1)
        listy_im_pole_aproxG1.append(y_im_pole_aproxG1)   
        

#        y_re_anaG2 = G2_real_ana(R0,value)
#        y_im_anaG2 = G2_imag_ana(R0,value)        


        y_re_numG2 = G2_real_num(R0,value)
        y_im_numG2 = G2_imag_num(R0,value)  
        
        
        y_re_pole_aproxG2 = G2_real_pole_aprox(R0,value)
        y_im_pole_aproxG2 = G2_imag_pole_aprox(R0,value)  
        
        
#        listy_re_anaG2.append(y_re_anaG2)
#        listy_im_anaG2.append(y_im_anaG2)
        
        listy_re_numG2.append(y_re_numG2)
        listy_im_numG2.append(y_im_numG2)   
        
        listy_re_pole_aproxG2.append(y_re_pole_aproxG2)
        listy_im_pole_aproxG2.append(y_im_pole_aproxG2) 

        err1_num = np.abs(y_re_numG1 - y_re_pole_aproxG1)
        err2_num = np.abs(y_re_numG2 - y_re_pole_aproxG2)
        
        print('error 1 real', err1_num)
        print('error 2 real', err2_num)
        j = j + 1
   
    graph(title,labelx,'Re($G_1$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#    plt.plot(listx,listy_re_anaG1,'.',ms = ms,color = 'purple',label = 'PP analytical')
    plt.plot(listx,listy_re_pole_aproxG1,'.',ms = ms,color = 'darkred',label = 'PP numerical')
    plt.plot(listx,listy_re_numG1,'.-',ms = ms,color = 'lightseagreen',label = 'full numerical')
    plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'Re_G1' + label1 + '.png', format='png')   

    graph(title,labelx,'Im($G_1$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#    plt.plot(listx,listy_im_anaG1,'.',ms = ms,color = 'purple',label = 'PP analytical')
    plt.plot(listx,listy_im_pole_aproxG1,'.',ms = ms,color = 'darkred',label = 'PP numerical')
    plt.plot(listx,listy_im_numG1,'.-',ms = ms,color = 'lightseagreen',label = 'full numerical')
    plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'Im_G1' + label1 + '.png', format='png')   



    graph(title,labelx,'Re($G_2$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#    plt.plot(listx,listy_re_anaG2,'.',ms = ms,color = 'purple',label = 'PP analytical')
    plt.plot(listx,listy_re_pole_aproxG2,'.',ms = ms,color = 'darkred',label = 'PP numerical')
    plt.plot(listx,listy_re_numG2,'.-',ms = ms,color = 'lightseagreen',label = 'full numerical')
    plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'Re_G2' + label1 + '.png', format='png')   

    graph(title,labelx,'Im($G_2$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#    plt.plot(listx,listy_im_anaG2,'.',ms = ms,color = 'purple',label = 'PP analytical')
    plt.plot(listx,listy_im_pole_aproxG2,'.',ms = ms,color = 'darkred',label = 'PP numerical')
    plt.plot(listx,listy_im_numG2,'.-',ms = ms,color = 'lightseagreen',label = 'full numerical')
    plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'Im_G2' + label1 + '.png', format='png') 
    
#%%
