
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

plot_vs_R = 1
plot_vs_E = 0

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_ctes =  path_basic.replace('/' + 'E_field/closing_contourns/separating_integrals','')
path_save = path_basic + '/' + 'E_pols'
#if not os.path.exists(path_save):
#    print('Creating folder to save graphs')
#    os.mkdir(path_save)

err = 'E_pols.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from E_pols import Ex_pols_num_G1, Ex_pols_num_G2, Ex_pols_num_G4, Ex_pols_num_G5, Ex_pols_ana
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

px,py,pz = 1,0,0

cte_phi0 = 8
if cte_phi0!=0:
    phi0 = np.pi/cte_phi0
    title0 = r'px = %i, py = %i, pz = %i, $\varphi$=$\pi$/%i' %(px,py,pz,cte_phi0)
else:
    phi0 = 0 
    title0 = r'px = %i, py = %i, pz = %i, $\varphi$=0' %(px,py,pz)

title1 = r'z = %inm, $z_p$=%inm, $\epsilon_2$ = %i, d = %inm' %(z0*1e3,zp*1e3,epsilon2,d*1e3)
labelp = r'_d%inm_px%i_py%i_pz%i_phi%i' %(d*1e3,px,py,pz,phi0)
title = title1 + '\n' + title0

N = 100
if plot_vs_R == 1: 
    omega_omegaWG0 = 0.8 # meV 
    labelx = 'R [nm]'  
    title = title + ', ' + r'$\omega/\omega_{WG}$ = %.2f' %(omega_omegaWG0)
    label1 = 'vs_R' + labelp + '_omega_omegaWG%i' %(omega_omegaWG0)
    listx = np.linspace(50,3000,N)
else:
    R0 = 500 #nanos
    # z0 = 0.06*1e3
    labelx = '$\omega/\omega_{WG}$'   
    title = title + ', ' + 'R = %inm' %(R0)
    label1 = 'vs_omega_omega_WG' + labelp
    listx = np.linspace(0.01,4,N)

#omegaWG = (np.pi/d)/(np.sqrt(epsilon2-1)) # el c esta metido en aux


#%%

def Ex_real_numG1(R_nano,omega_omegaWG):
#    omega0 = energy0*1e-3/aux 
#    omega_omegaWG = omega0/omegaWG
    R = R_nano*1e-3

    rta = Ex_pols_num_G1(omega_omegaWG,d,R,phi0,z0,zp,px,py,pz)
    
    return rta.real

def Ex_imag_numG1(R_nano,omega_omegaWG):
#    omega0 = energy0*1e-3/aux    
#    omega_omegaWG = omega0/omegaWG
    R = R_nano*1e-3
    
    rta = Ex_pols_num_G1(omega_omegaWG,d,R,phi0,z0,zp,px,py,pz)
    return rta.imag

def Ex_real_numG2(R_nano,omega_omegaWG):
#    omega0 = energy0*1e-3/aux    
#    omega_omegaWG = omega0/omegaWG
    R = R_nano*1e-3

    rta = Ex_pols_num_G2(omega_omegaWG,d,R,phi0,z0,zp,px,py,pz)
    
    return rta.real

def Ex_imag_numG2(R_nano,omega_omegaWG):
#    omega0 = energy0*1e-3/aux    
#    omega_omegaWG = omega0/omegaWG
    R = R_nano*1e-3
    
    rta = Ex_pols_num_G2(omega_omegaWG,d,R,phi0,z0,zp,px,py,pz)
    return rta.imag


def Ex_real_numG4(R_nano,omega_omegaWG):
#    omega0 = energy0*1e-3/aux    
#    omega_omegaWG = omega0/omegaWG
    R = R_nano*1e-3

    rta = Ex_pols_num_G4(omega_omegaWG,d,R,phi0,z0,zp,px,py,pz)
    
    return rta.real

def Ex_imag_numG4(R_nano,omega_omegaWG):
#    omega0 = energy0*1e-3/aux    
#    omega_omegaWG = omega0/omegaWG
    R = R_nano*1e-3
    
    rta = Ex_pols_num_G4(omega_omegaWG,d,R,phi0,z0,zp,px,py,pz)
    return rta.imag


def Ex_real_numG5(R_nano,omega_omegaWG):
#    omega0 = energy0*1e-3/aux    
#    omega_omegaWG = omega0/omegaWG
    R = R_nano*1e-3

    rta = Ex_pols_num_G5(omega_omegaWG,d,R,phi0,z0,zp,px,py,pz)
    
    return rta.real

def Ex_imag_numG5(R_nano,omega_omegaWG):
#    omega0 = energy0*1e-3/aux    
#    omega_omegaWG = omega0/omegaWG
    R = R_nano*1e-3
    
    rta = Ex_pols_num_G5(omega_omegaWG,d,R,phi0,z0,zp,px,py,pz)
    return rta.imag

#%%

def Ex_real_anaG2(R_nano,omega_omegaWG):
#    omega0 = energy0*1e-3/aux 
#    omega_omegaWG = omega0/omegaWG
    R = R_nano*1e-3

    rta = Ex_pols_ana(omega_omegaWG,d,R,phi0,z0,zp,px,py,pz)
    
    return rta.real

def Ex_imag_anaG2(R_nano,omega_omegaWG):
#    omega0 = energy0*1e-3/aux    
#    omega_omegaWG = omega0/omegaWG
    R = R_nano*1e-3
    
    rta = Ex_pols_ana(omega_omegaWG,d,R,phi0,z0,zp,px,py,pz)
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
    


    listy_re_numG1 = []
    listy_im_numG1 = []

    listy_re_numG2 = []
    listy_im_numG2 = []

    listy_re_numG4 = []
    listy_im_numG4 = []

    listy_re_numG5 = []
    listy_im_numG5 = []
    
    listy_re_Gtot = []
    listy_im_Gtot = []

    listy_re_anaG2 = []
    listy_im_anaG2 = []
    
    j = 0
    for value in listx: 



        y_re_numG1 = Ex_real_numG1(value,omega_omegaWG0)
        y_im_numG1 = Ex_imag_numG1(value,omega_omegaWG0)  
        

        y_re_numG2 = Ex_real_numG2(value,omega_omegaWG0)
        y_im_numG2 = Ex_imag_numG2(value,omega_omegaWG0)  


        y_re_numG4 = Ex_real_numG4(value,omega_omegaWG0)
        y_im_numG4 = Ex_imag_numG4(value,omega_omegaWG0)  


        y_re_numG5 = Ex_real_numG5(value,omega_omegaWG0)
        y_im_numG5 = Ex_imag_numG5(value,omega_omegaWG0)  



        y_re_anaG2 = Ex_real_anaG2(value,omega_omegaWG0)
        y_im_anaG2 = Ex_imag_anaG2(value,omega_omegaWG0)  
        

        listy_re_numG1.append(y_re_numG1)
        listy_im_numG1.append(y_im_numG1)   
                
        
        listy_re_numG2.append(y_re_numG2)
        listy_im_numG2.append(y_im_numG2)   
        

        listy_re_anaG2.append(y_re_anaG2)
        listy_im_anaG2.append(y_im_anaG2)
        
        listy_re_numG4.append(y_re_numG4)
        listy_im_numG4.append(y_im_numG4)   

        
        listy_re_numG5.append(y_re_numG5)
        listy_im_numG5.append(y_im_numG5)   
        
        
        listy_re_Gtot.append( y_re_numG1  + y_re_numG4 + y_re_numG5)
        listy_im_Gtot.append( y_im_numG1  + y_re_numG4 + y_re_numG5)

#       
        print(j)
#        
        j = j + 1
        

        

if plot_vs_E == 1: 
    
    listy_re_numG1 = []
    listy_im_numG1 = []

    listy_re_numG2 = []
    listy_im_numG2 = []

    listy_re_numG4 = []
    listy_im_numG4 = []

    listy_re_numG5 = []
    listy_im_numG5 = []
   
    listy_re_Gtot = []
    listy_im_Gtot = []

    listy_re_anaG2 = []
    listy_im_anaG2 = []
    
    j = 0
    for value in listx: 



        y_re_numG1 = Ex_real_numG1(R0,value)
        y_im_numG1 = Ex_imag_numG1(R0,value)  
        

        y_re_numG2 = Ex_real_numG2(R0,value)
        y_im_numG2 = Ex_imag_numG2(R0,value)  


        y_re_numG4 = Ex_real_numG4(R0,value)
        y_im_numG4 = Ex_imag_numG4(R0,value)  


        y_re_numG5 = Ex_real_numG5(R0,value)
        y_im_numG5 = Ex_imag_numG5(R0,value)  



        y_re_anaG2 = Ex_real_anaG2(R0,value)
        y_im_anaG2 = Ex_imag_anaG2(R0,value)  
        

        listy_re_numG1.append(y_re_numG1)
        listy_im_numG1.append(y_im_numG1)   
                
        
        listy_re_numG2.append(y_re_numG2)
        listy_im_numG2.append(y_im_numG2)   
        

        listy_re_anaG2.append(y_re_anaG2)
        listy_im_anaG2.append(y_im_anaG2)
        
        listy_re_numG4.append(y_re_numG4)
        listy_im_numG4.append(y_im_numG4)   

        
        listy_re_numG5.append(y_re_numG5)
        listy_im_numG5.append(y_im_numG5)   
       
        listy_re_Gtot.append(  y_re_numG4 + y_re_numG5)
        listy_im_Gtot.append( y_re_numG4 + y_re_numG5)
        
   #%%
   
graph(title,labelx,'Re($G_1$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(listx,listy_re_numG1,'.-',ms = ms,color = 'lightseagreen',label = 'full numerical')
#plt.plot(listx,listy_re_pole_aproxG1,'.-',ms = ms,color = 'darkorange',label = 'PP numerical')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
os.chdir(path_save)
plt.savefig( 'Re_G1' + label1 + '.png', format='png')   

graph(title,labelx,'Im($G_1$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(listx,listy_im_numG1,'.-',ms = ms,color = 'lightseagreen',label = 'full numerical')
#plt.plot(listx,listy_im_pole_aproxG1,'.-',ms = ms,color = 'darkorange',label = 'PP numerical')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
plt.savefig( 'Im_G1' + label1 + '.png', format='png')   
    

#%%

graph(title,labelx,'Re($G_2$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(listx,listy_re_numG2,'.-',ms = ms,color = 'lightseagreen',label = 'full numerical')
plt.plot(listx,listy_re_anaG2,'.-',ms = ms,color = 'purple',label = 'PP analytical')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
plt.savefig( 'Re_G2' + label1 + '.png', format='png')   

graph(title,labelx,'Im($G_2$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(listx,listy_im_numG2,'.-',ms = ms,color = 'lightseagreen',label = 'full numerical')
plt.plot(listx,listy_im_anaG2,'.-',ms = ms,color = 'purple',label = 'PP analytical')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
plt.savefig( 'Im_G2' + label1 + '.png', format='png')   


#%%


graph(title,labelx,'Re($G_4$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(listx,listy_re_numG4,'.-',ms = ms,color = 'lightseagreen',label = 'full numerical')
#plt.plot(listx,listy_re_pole_aproxG1,'.-',ms = ms,color = 'darkorange',label = 'PP numerical')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
os.chdir(path_save)
plt.savefig( 'Re_G4' + label1 + '.png', format='png')   

graph(title,labelx,'Im($G_4$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(listx,listy_im_numG4,'.-',ms = ms,color = 'lightseagreen',label = 'full numerical')
#plt.plot(listx,listy_im_pole_aproxG1,'.-',ms = ms,color = 'darkorange',label = 'PP numerical')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
plt.savefig( 'Im_G4' + label1 + '.png', format='png')   


#%%


graph(title,labelx,'Re($G_5$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(listx,listy_re_numG5,'.-',ms = ms,color = 'lightseagreen',label = 'full numerical')
#plt.plot(listx,listy_re_pole_aproxG1,'.-',ms = ms,color = 'darkorange',label = 'PP numerical')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
os.chdir(path_save)
plt.savefig( 'Re_G5' + label1 + '.png', format='png')   

graph(title,labelx,'Im($G_5$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(listx,listy_im_numG5,'.-',ms = ms,color = 'lightseagreen',label = 'full numerical')
#plt.plot(listx,listy_im_pole_aproxG1,'.-',ms = ms,color = 'darkorange',label = 'PP numerical')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
plt.savefig( 'Im_G5' + label1 + '.png', format='png')   

#%%



graph(title,labelx,'Re($G_{tot}$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(listx,listy_re_Gtot,'.-',ms = ms,color = 'lightseagreen',label = 'full numerical')
plt.plot(listx,listy_re_pole_aproxG1,'.-',ms = ms,color = 'darkorange',label = 'PP numerical')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
os.chdir(path_save)
plt.savefig( 'Re_Gtot' + label1 + '.png', format='png')   

graph(title,labelx,'Im($G_{tot}$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(listx,listy_im_Gtot,'.-',ms = ms,color = 'lightseagreen',label = 'full numerical')
plt.plot(listx,listy_im_pole_aproxG1,'.-',ms = ms,color = 'darkorange',label = 'PP numerical')
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
plt.savefig( 'Im_Gtot' + label1 + '.png', format='png')   