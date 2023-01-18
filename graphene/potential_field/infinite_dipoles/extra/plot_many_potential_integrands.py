
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

plot_vs_y = 0
plot_vs_E = 1

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/potential_field/infinite_dipoles','')
path_save = path_basic + '/' + 'integrands'
if not os.path.exists(path_save):
    print('Creating folder to save graphs')
    os.mkdir(path_save)

err = 'many_potential_integrals.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from many_potential_integrands import G1_num, G1_pole_aprox, G2_num, G2_pole_aprox,G3_num_odd, G3_pole_aprox_odd
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

#v = c/int_v
#omega = 0.7*1e12
cota = 75 #nanometros
epsi1,epsi2 = 1,1
hbmu,hbgama = 0.3,0.0001
zp = 0.05

y0 = 5
int_v = 10

title2 = r'$\hbar\mu$ = %.2feV, $\hbar\gamma$ = %.4feV, v = c/%i' %(hbmu,hbgama,int_v) 
title4 = r'$z_p$=%inm' %(zp*1e3)  + ', ' + 'y = %inm' %(y0*1e3)

N = 2000
    # z0 = 0.06*1e3
labelx = r'$\alpha_y$'   
label1 = 'vs_alpha_y_y0%inm' %(y0*1e3) 
cota_alpha_y = 600
list_alpha_y = np.linspace(-cota_alpha_y,cota_alpha_y,N)
list_E = [30,60]


title = title2 + '\n' + title4 

#%%

def function_G1_real_num(alpha_y,energy0):
    omegac0 = energy0*1e-3/aux 
    rta = G1_num(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,zp,y0,alpha_y)
    # print(rta)
    return rta.real

def function_G1_imag_num(alpha_y,energy0):
    omegac0 = energy0*1e-3/aux 
    rta = G1_num(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,zp,y0,alpha_y)
    return rta.imag

#%%
    

def function_G2_real_num(alpha_y,energy0):
    omegac0 = energy0*1e-3/aux 
    rta = G2_num(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,zp,y0,alpha_y)
    # print(rta)
    return rta.real

def function_G2_imag_num(alpha_y,energy0):
    omegac0 = energy0*1e-3/aux 
    rta = G2_num(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,zp,y0,alpha_y)
    return rta.imag

#%%

def function_G3_real_num(alpha_y,energy0):
    omegac0 = energy0*1e-3/aux 
    rta = G3_num_odd(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,zp,y0,alpha_y)
    # print(rta)
    return rta.real

def function_G3_imag_num(alpha_y,energy0):
    omegac0 = energy0*1e-3/aux 
    rta = G3_num_odd(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,zp,y0,alpha_y)
    return rta.imag

#%%
    
def function_G1_real_pole_aprox(alpha_y,energy0):
    omegac0 = energy0*1e-3/aux 


    rta = G1_pole_aprox(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,zp,y0,alpha_y)
    # print(rta)
    return rta.real

def function_G1_imag_pole_aprox(alpha_y,energy0):
    omegac0 = energy0*1e-3/aux 
    
    
    rta = G1_pole_aprox(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,zp,y0,alpha_y)
    return rta.imag


def function_G2_real_pole_aprox(alpha_y,energy0):
    omegac0 = energy0*1e-3/aux 
    rta = G2_pole_aprox(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,zp,y0,alpha_y)
    # print(rta)
    return rta.real

def function_G2_imag_pole_aprox(alpha_y,energy0):
    omegac0 = energy0*1e-3/aux     
    rta = G2_pole_aprox(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,zp,y0,alpha_y)
    return rta.imag


def function_G3_real_pole_aprox(alpha_y,energy0):
    omegac0 = energy0*1e-3/aux 
    rta = G3_pole_aprox_odd(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,zp,y0,alpha_y)
    # print(rta)
    return rta.real

def function_G3_imag_pole_aprox(alpha_y,energy0):
    omegac0 = energy0*1e-3/aux     
    rta = G3_pole_aprox_odd(omegac0,epsi1,epsi2,hbmu,hbgama,int_v,zp,y0,alpha_y)
    return rta.imag

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
length_marker = 1

def graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad):
    plt.figure(figsize=tamfig)
    plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
    plt.ylabel(labely,fontsize=tamletra,labelpad =labelpadx)
    plt.tick_params(labelsize = tamnum, pad = pad)
    plt.title(title,fontsize=int(tamtitle*0.9))

    return   

#%%

#listy_re_anaG1_tot = []
#listy_im_anaG1_tot = [] 

listy_re_numG1_tot = []
listy_im_numG1_tot = []

listy_re_pole_approx_G1_tot = []   
listy_im_pole_approx_G1_tot = []   



listy_re_numG2_tot = []
listy_im_numG2_tot = []

listy_re_pole_approx_G2_tot = []   
listy_im_pole_approx_G2_tot = []   

zeros_tot_G3 = []

listy_re_numG3_tot = []
listy_im_numG3_tot = []

listy_re_pole_approx_G3_tot = []   
listy_im_pole_approx_G3_tot = []   

for energy0 in list_E:
    
    listy_re_anaG1 = []
    listy_im_anaG1 = []
    
    listy_re_numG1 = []
    listy_im_numG1 = []
 
    listy_re_pole_approx_G1 = []   
    listy_im_pole_approx_G1 = []   


    listy_re_anaG2 = []
    listy_im_anaG2 = []
    
    listy_re_numG2 = []
    listy_im_numG2 = []
 
    listy_re_pole_approx_G2 = []   
    listy_im_pole_approx_G2 = []   
    

    listy_re_anaG3 = []
    listy_im_anaG3 = []
    
    listy_re_numG3 = []
    listy_im_numG3 = []
 
    listy_re_pole_approx_G3 = []   
    listy_im_pole_approx_G3 = []   
    
    zeros_G3 = []
    omegac0 = energy0*1e-3/aux 
    
    nmax = int(cota_alpha_y*y0*omegac0/np.pi - 0.5)
    
    n = 0

    for alpha_y in list_alpha_y:     
#        
        y_re_numG1 = function_G1_real_num(alpha_y,energy0)
        y_im_numG1 = function_G1_imag_num(alpha_y,energy0)
        
        y_re_pole_approxG1 = function_G1_real_pole_aprox(alpha_y,energy0)
        y_im_pole_approxG1 = function_G1_imag_pole_aprox(alpha_y,energy0)
        
        listy_re_numG1.append(y_re_numG1)
        listy_im_numG1.append(y_im_numG1)

        listy_re_pole_approx_G1.append(y_re_pole_approxG1)
        listy_im_pole_approx_G1.append(y_im_pole_approxG1)
        
        
        y_re_numG2 = function_G2_real_num(alpha_y,energy0)
        y_im_numG2 = function_G2_imag_num(alpha_y,energy0)
        
        y_re_pole_approxG2 = function_G2_real_pole_aprox(alpha_y,energy0)
        y_im_pole_approxG2 = function_G2_imag_pole_aprox(alpha_y,energy0)
        
        listy_re_numG2.append(y_re_numG2)
        listy_im_numG2.append(y_im_numG2)

        listy_re_pole_approx_G2.append(y_re_pole_approxG2)
        listy_im_pole_approx_G2.append(y_im_pole_approxG2)
                

        y_re_numG3 = function_G3_real_num(alpha_y,energy0)
        y_im_numG3 = function_G3_imag_num(alpha_y,energy0)
        
        y_re_pole_approxG3 = function_G3_real_pole_aprox(alpha_y,energy0)
        y_im_pole_approxG3 = function_G3_imag_pole_aprox(alpha_y,energy0)
        
        listy_re_numG3.append(y_re_numG3)
        listy_im_numG3.append(y_im_numG3)

        listy_re_pole_approx_G3.append(y_re_pole_approxG3)
        listy_im_pole_approx_G3.append(y_im_pole_approxG3)

        if n <= nmax:
            calcular_zeros_G3 = (2*n+1)*np.pi*0.5/(omegac0*y0)
            zeros_G3.append(calcular_zeros_G3)
        
        n = n + 1
        
    zeros_tot_G3.append(zeros_G3)


        
    listy_re_pole_approx_G1_tot.append(listy_re_pole_approx_G1)
    listy_im_pole_approx_G1_tot.append(listy_im_pole_approx_G1)  
    
    listy_re_numG1_tot.append(listy_re_numG1)
    listy_im_numG1_tot.append(listy_im_numG1)



    listy_re_pole_approx_G2_tot.append(listy_re_pole_approx_G2)
    listy_im_pole_approx_G2_tot.append(listy_im_pole_approx_G2)  
    
    listy_re_numG2_tot.append(listy_re_numG2)
    listy_im_numG2_tot.append(listy_im_numG2)


    listy_re_pole_approx_G3_tot.append(listy_re_pole_approx_G3)
    listy_im_pole_approx_G3_tot.append(listy_im_pole_approx_G3)  
    
    listy_re_numG3_tot.append(listy_re_numG3)
    listy_im_numG3_tot.append(listy_im_numG3)


#%%        
        
graph(title,labelx,'Re($G_1$) E = 30meV',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)

plt.plot(list_alpha_y,listy_re_pole_approx_G1_tot[0],'.',ms = 4,color = 'lightseagreen',label = 'PP numerical')
plt.plot(list_alpha_y,listy_re_numG1_tot[0],'.-',ms = 2,color = 'darkred',label = 'full numerical')
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
os.chdir(path_save)
plt.savefig( 'Re_G1_E30' + label1 + '.png', format='png')   

graph(title,labelx,'Re($G_1$) E = 60meV',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(list_alpha_y,listy_re_pole_approx_G1_tot[1],'.',ms = 4,color = 'lightseagreen',label = 'PP numerical')
plt.plot(list_alpha_y,listy_re_numG1_tot[1],'.-',ms = 2,color = 'darkred',label = 'full numerical')
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
os.chdir(path_save)
plt.savefig( 'Re_G1_E60' + label1 + '.png', format='png')   
   
graph(title,labelx,'Im($G_1$) E = 30meV',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(list_alpha_y,listy_im_pole_approx_G1_tot[0],'.',ms = 4,color = 'lightseagreen',label = 'PP numerical')
plt.plot(list_alpha_y,listy_im_numG1_tot[0],'.-',ms = 2,color = 'darkred',label = 'full numerical')
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
os.chdir(path_save)
plt.savefig( 'Im_G1_E30' + label1 + '.png', format='png')   



graph(title,labelx,'Im($G_1$) E = 60meV',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(list_alpha_y,listy_im_pole_approx_G1_tot[1],'.',ms = 4,color = 'lightseagreen',label = 'PP numerical')
plt.plot(list_alpha_y,listy_im_numG1_tot[1],'.-',ms = 2,color = 'darkred',label = 'full numerical')
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
os.chdir(path_save)
plt.savefig( 'Im_G1_E60' + label1 + '.png', format='png') 

#%%
graph(title,labelx,'Re($G_2$) E = 30meV',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(list_alpha_y,listy_re_pole_approx_G2_tot[0],'.',ms = 4,color = 'lightseagreen',label = 'PP numerical')
plt.plot(list_alpha_y,listy_re_numG2_tot[0],'.-',ms = 2,color = 'darkred',label = 'full numerical')
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
os.chdir(path_save)
plt.savefig( 'Re_G2_E30' + label1 + '.png', format='png')   

graph(title,labelx,'Re($G_2$) E = 60meV',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(list_alpha_y,listy_re_pole_approx_G2_tot[1],'.',ms = 4,color = 'lightseagreen',label = 'PP numerical')
plt.plot(list_alpha_y,listy_re_numG2_tot[1],'.-',ms = 2,color = 'darkred',label = 'full numerical')
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
os.chdir(path_save)
plt.savefig( 'Re_G2_E60' + label1 + '.png', format='png')   
   
graph(title,labelx,'Im($G_2$) E = 30meV',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(list_alpha_y,listy_im_pole_approx_G2_tot[0],'.',ms = 4,color = 'lightseagreen',label = 'PP numerical')
plt.plot(list_alpha_y,listy_im_numG2_tot[0],'.-',ms = 2,color = 'darkred',label = 'full numerical')
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
os.chdir(path_save)
plt.savefig( 'Im_G2_E30' + label1 + '.png', format='png')   


graph(title,labelx,'Im($G_2$) E = 60meV',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(list_alpha_y,listy_im_pole_approx_G2_tot[1],'.',ms = 4,color = 'lightseagreen',label = 'PP numerical')
plt.plot(list_alpha_y,listy_im_numG2_tot[1],'.-',ms = 2,color = 'darkred',label = 'full numerical')
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
os.chdir(path_save)
plt.savefig( 'Im_G2_E60' + label1 + '.png', format='png') 

#%%

graph(title,labelx,'Re($G_3$) E = 30meV',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(list_alpha_y,listy_re_pole_approx_G3_tot[0],'.',ms = 4,color = 'lightseagreen',label = 'PP numerical')
plt.plot(list_alpha_y,listy_re_numG3_tot[0],'.-',ms = 2,color = 'darkred',label = 'full numerical')
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
os.chdir(path_save)
plt.savefig( 'Re_G3_E30' + label1 + '.png', format='png')   

graph(title,labelx,'Re($G_3$) E = 60meV',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(list_alpha_y,listy_re_pole_approx_G3_tot[1],'.',ms = 4,color = 'lightseagreen',label = 'PP numerical')
plt.plot(list_alpha_y,listy_re_numG3_tot[1],'.-',ms = 2,color = 'darkred',label = 'full numerical')
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
os.chdir(path_save)
plt.savefig( 'Re_G3_E60' + label1 + '.png', format='png')   
   
graph(title,labelx,'Im($G_3$) E = 30meV',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(list_alpha_y,listy_im_pole_approx_G3_tot[0],'.',ms = 4,color = 'lightseagreen',label = 'PP numerical')
plt.plot(list_alpha_y,listy_im_numG3_tot[0],'.-',ms = 2,color = 'darkred',label = 'full numerical')
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
os.chdir(path_save)
plt.savefig( 'Im_G3_E30' + label1 + '.png', format='png')   


graph(title,labelx,'Im($G_3$) E = 60meV',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(list_alpha_y,listy_im_pole_approx_G3_tot[1],'.',ms = 4,color = 'lightseagreen',label = 'PP numerical')
plt.plot(list_alpha_y,listy_im_numG3_tot[1],'.-',ms = 2,color = 'darkred',label = 'full num')
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
os.chdir(path_save)
plt.savefig( 'Im_G3_E60' + label1 + '.png', format='png') 


