
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: leila

"""

import numpy as np
import sys
import os 
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()

#%%

plot_vs_y = 1
plot_vs_x = 0

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/potential_field/infinite_dipoles/num_vs_ana','')
path_2 =  path_basic.replace('/num_vs_ana','')

path_save = path_basic + '/' + 'many_potential2'
if not os.path.exists(path_save):
    print('Creating folder to save graphs')
    os.mkdir(path_save)

err = 'potential_many_dipoles.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from potential_many_dipoles2 import phi_many_dipoles_num2, phi_many_dipoles_ana2,phi_many_dipoles_pole_approx,integral_neglecting
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
#v = c/int_v
#omega = 0.7*1e12
cota = 75 #nanometros
epsi1,epsi2 = 1,1
hbmu,hbgama = 0.3,0.0001
zp = 0.05
z0 = 0

#Ndipoles = 50
#list_xD = np.linspace(-0.01,0.01,Ndipoles)
a = 0.05

px = 1
py = 0
pz = 0

#omega0THz = 600
#omega0 = omega0THz*1e12 
 
title1 = r'v = c/%i, a = %inm' %(int_v,a*1e3)     
title2 = r'$\hbar\mu$ = %.2feV, $\hbar\gamma$ = %.1fmeV' %(hbmu,hbgama*1e3) 
title3 = r'$z_p$=%inm, px=%i, py=%i, pz=%i, z = %inm' %(zp*1e3,px,py,pz,z0*1e3)


N = 100
if plot_vs_x == 1: 
    energy0 = 75*1e-3 # meV 
    y0 = 0.01
    omegac0 = energy0/(hb*c)
#    lambbda = 2*np.pi/omegac
    labelp = r'_px%i_py%i_pz%i_a%inm_E%imeV' %(px,py,pz,a*1e3,energy0*1e3)
    labelx = 'x [nm]'  
    title4 = 'E = %.2f meV, y=%inm' %(energy0*1e3,y0*1e3)
    label1 = 'vs_x' + labelp
else:
    energy0 = 43*1e-3 # meV 
    x0 = 0
    
    omegac0 = energy0/(hb*c)
#    lambbda = 2*np.pi/omegac
    labelp = r'_px%i_py%i_pz%i_a%inm_E%imeV' %(px,py,pz,a*1e3,energy0*1e3)
    labelx = 'y [nm]'  
    title4 = 'E = %.2f meV, x=%inm' %(energy0*1e3,x0*1e3)
    label1 = '_vs_y' + labelp


listx = np.linspace(10,3000,N)
title = title1 + ', '  + title2 + '\n' + title3 + '\n'  + title4 

def function_real_num(x0,y0):
    y0 = y0*1e-3 
    x0 = x0*1e-3 

    rta = phi_many_dipoles_num2(omegac0,epsi1,epsi2,hbmu,hbgama,x0,y0,z0,a,zp,int_v,px,py,pz)
    # print(rta)
    return rta.real

def function_imag_num(x0,y0):
    y0 = y0*1e-3 
    x0 = x0*1e-3 
    
    rta = phi_many_dipoles_num2(omegac0,epsi1,epsi2,hbmu,hbgama,x0,y0,z0,a,zp,int_v,px,py,pz)
    return rta.imag

#%%

def function_real_ana(x0,y0):
    y0 = y0*1e-3 
    x0 = x0*1e-3 

    rta = phi_many_dipoles_ana2(omegac0,epsi1,epsi2,hbmu,hbgama,x0,y0,z0,a,zp,int_v,px,py,pz)
    
    return rta.real

def function_imag_ana(x0,y0):
    y0 = y0*1e-3 
    x0 = x0*1e-3 
    
    rta = phi_many_dipoles_ana2(omegac0,epsi1,epsi2,hbmu,hbgama,x0,y0,z0,a,zp,int_v,px,py,pz)
    return rta.imag

#%%
    
def function_real_pole_approx(x0,y0):
    y0 = y0*1e-3 
    x0 = x0*1e-3 

    rta = phi_many_dipoles_pole_approx(omegac0,epsi1,epsi2,hbmu,hbgama,x0,y0,z0,a,zp,int_v,px,py,pz)
    
    return rta.real

def function_imag_pole_approx(x0,y0):
    y0 = y0*1e-3 
    x0 = x0*1e-3 
    
    rta = phi_many_dipoles_pole_approx(omegac0,epsi1,epsi2,hbmu,hbgama,x0,y0,z0,a,zp,int_v,px,py,pz)
    return rta.imag


#%%
    
def function_real_integral_neg(x0,y0):
    y0 = y0*1e-3 
    x0 = x0*1e-3 

    rta = integral_neglecting(omegac0,epsi1,epsi2,hbmu,hbgama,x0,y0,z0,a,zp,int_v,px,py,pz)
    
    return rta.real

def function_imag_integral_neg(x0,y0):
    y0 = y0*1e-3 
    x0 = x0*1e-3 
    
    rta = integral_neglecting(omegac0,epsi1,epsi2,hbmu,hbgama,x0,y0,z0,a,zp,int_v,px,py,pz)
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
ms = 5
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

if plot_vs_x == 1: 
    listy_re_ana = []
    listy_im_ana = []
    
    listy_re_num = []
    listy_im_num = []
    
    listy_re_pp_aprox = []
    listy_im_pp_aprox = []
    
    for value in listx: 

        y_re_ana = function_real_ana(value,y0)
        y_im_ana = function_imag_ana(value,y0)        

        y_re_num = function_real_num(value,y0)
        y_im_num = function_imag_num(value,y0)

        y_re_pp = function_real_pole_approx(value,y0)
        y_im_pp = function_imag_pole_approx(value,y0)

       
        listy_re_ana.append(y_re_ana)
        listy_im_ana.append(y_im_ana)
        
        listy_re_num.append(y_re_num)
        listy_im_num.append(y_im_num)   

        listy_re_pp_aprox.append(y_re_pp)
        listy_im_pp_aprox.append(y_im_pp)   

 
    graph(title,labelx,'Re($\phi$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    plt.plot(listx,listy_re_ana,'.-',ms = ms,color = 'purple',label = 'analytical')
    plt.plot(listx,listy_re_num,'.-',ms = ms,color = 'lightseagreen',label = 'sin aprox')
    plt.plot(listx,listy_re_pp_aprox,'.',ms = 3,color = 'darkred',label = 'pole aprox')
    plt.legend(loc = 'best',markerscale=1.5,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
#    plt.yscale('log')
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'Re_phi' + label1 + '.png', format='png')   
   
    graph(title,labelx,'Im($\phi$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    plt.plot(listx,listy_im_ana,'.',ms = ms,color = 'purple',label = 'analytical')
#    plt.plot(listx,listy_im_num,'.-',ms = 3,color = 'lightseagreen',label = 'sin aprox')
    plt.legend(loc = 'best',markerscale=1.5,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
#    plt.yscale('log')
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'Im_phi_ana' + label1 + '.png', format='png')   

    graph(title,labelx,'Im($\phi$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    plt.plot(listx,listy_im_ana,'.',ms = ms,color = 'purple',label = 'analytical')
    plt.plot(listx,listy_im_num,'.-',ms = ms,color = 'lightseagreen',label = 'sin aprox')
    plt.plot(listx,listy_im_pp_aprox,'.',ms = 3,color = 'darkred',label = 'pole aprox')
    plt.legend(loc = 'best',markerscale=1.5,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
#    plt.yscale('log')
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'Im_phi_num' + label1 + '.png', format='png')       

    
if plot_vs_y == 1: 
    
    listy_re_ana = []
    listy_im_ana = []
    listy_abs_ana = []
    
    listy_re_num = []
    listy_im_num = []
    listy_abs_num = []

    listy_re_pp_aprox = []
    listy_im_pp_aprox = []
    listy_abs_pp_aprox = []
    
    listy_int_neg = []
    
    

    for value in listx: 

        y_re_ana = function_real_ana(x0,value)
        y_im_ana = function_imag_ana(x0,value)        

        y_re_num = function_real_num(x0,value)
        y_im_num = function_imag_num(x0,value)

        y_re_pp = function_real_pole_approx(x0,value)
        y_im_pp = function_imag_pole_approx(x0,value)
        
        listy_re_ana.append(y_re_ana)
        listy_im_ana.append(y_im_ana)
        listy_abs_ana.append(np.abs(y_im_ana + 1j*y_im_ana))
        
        
        listy_re_num.append(y_re_num)
        listy_im_num.append(y_im_num)   
        listy_abs_num.append(np.abs(y_re_num + 1j*y_im_num))


        listy_re_pp_aprox.append(y_re_pp)
        listy_im_pp_aprox.append(y_im_pp)   
        listy_abs_pp_aprox.append(np.abs(y_re_pp + 1j*y_im_pp))
 
    
        y_re_int = function_real_integral_neg(x0,value)
        y_im_int = function_imag_integral_neg(x0,value)    
        
        listy_int_neg.append(np.abs(y_re_int + 1j*y_im_int))
    
    
    graph(title,labelx,'Re($\phi$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    plt.plot(listx,listy_re_ana,'.-',ms = ms,color = 'purple',label = 'analytical')
    plt.plot(listx,listy_re_num,'.-',ms = ms,color = 'lightseagreen',label = 'numerical')
    plt.plot(listx,listy_re_pp_aprox,'.',ms = 3,color = 'darkred',label = 'pole aprox')
    plt.legend(loc = 'best',markerscale=1.5,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
#    plt.yscale('log')
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'Re_phi' + label1 + '.png', format='png')   

    graph(title,labelx,'Im($\phi$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    plt.plot(listx,listy_im_ana,'.',ms = ms,color = 'purple',label = 'analytical')
    plt.plot(listx,listy_im_num,'.-',ms = ms,color = 'lightseagreen',label = 'numerical')
    plt.plot(listx,listy_im_pp_aprox,'.',ms = 3,color = 'darkred',label = 'pole aprox')
    plt.legend(loc = 'best',markerscale=1.5,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
#    plt.yscale('log')
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'Im_phi_num' + label1 + '.png', format='png')       


    graph(title,labelx,'$|\phi|$',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    plt.plot(listx,listy_abs_ana,'.',ms = ms,color = 'purple',label = 'analytical')
    plt.plot(listx,listy_abs_num,'.-',ms = ms,color = 'lightseagreen',label = 'numerical')
    plt.plot(listx,listy_abs_pp_aprox,'.',ms = 3,color = 'darkred',label = 'pole aprox')
    plt.legend(loc = 'best',markerscale=1.5,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
#    plt.yscale('log')
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'Abs_phi_num' + label1 + '.png', format='png')       


    graph(title,labelx,'|integral neglected|',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    plt.plot(listx,listy_int_neg,'.',ms = ms,color = 'purple',label = 'analytical')
    plt.legend(loc = 'best',markerscale=1.5,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
#    plt.yscale('log')
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'Abs_integral' + label1 + '.png', format='png')       


#%%
