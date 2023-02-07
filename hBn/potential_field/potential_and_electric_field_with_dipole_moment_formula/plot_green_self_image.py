
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

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/potential_field/potential_and_electric_field_with_dipole_moment_formula','')
path_save = path_basic + '/' + 'Green_self_image'
if not os.path.exists(path_save):
    print('Creating folder to save graphs')
    os.mkdir(path_save)

err = 'fieldE_direct_numerical.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from green_self_image import green_self_ana,green_self_pole_aprox_parallel_num,green_self_pole_aprox, green_self_ana2, green_self_num, green_self_ana3,green_self_ana_parallel_num
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

epsi1,epsi3 = 1,1
zp = 0.05


d_nano = 10

title = r'$z_p$=%inm, d = %.2f nm' %(zp*1e3,d_nano)

N = 100

    # z0 = 0.06*1e3
labelx = r'$\hbar\omega$ (eV)'   
label1 = 'vs_E_d%.2fnm' %(d_nano) 

x1 = 0.09260651629072682 
x2 = 0.10112781954887218


x3 = 0.17030075187969923
x4 = 0.19937343358395992

listx = np.linspace(0.09,0.2,N)


def function_num_xx_re(energy0):
    omegac0 = energy0/aux 

    rtaself_x, rtaself_y, rtaself_z = green_self_num(omegac0,epsi1,epsi3,d_nano,zp)
    # print(rta    
    return np.real(rtaself_x)

def function_num_xx_im(energy0):
    omegac0 = energy0/aux 

    rtaself_x, rtaself_y, rtaself_z = green_self_num(omegac0,epsi1,epsi3,d_nano,zp)
    # print(rta    
    return np.imag(rtaself_x)


def function_ana_xx_re(energy0):
    omegac0 = energy0/aux 

    rtaself_x, rtaself_y, rtaself_z  = green_self_ana(omegac0,epsi1,epsi3,d_nano,zp)
    
    return np.real(rtaself_x)


def function_ana_xx_re_parallel_num(energy0):
    omegac0 = energy0/aux 

    rtaself_x, rtaself_y, rtaself_z  = green_self_ana_parallel_num(omegac0,epsi1,epsi3,d_nano,zp)
    
    return np.real(rtaself_x)


def function_ana2_xx_re(energy0):
    omegac0 = energy0/aux 

    rtaself_x, rtaself_y, rtaself_z  = green_self_ana2(omegac0,epsi1,epsi3,d_nano,zp)
    
    return np.real(rtaself_x)



def function_ana3_xx_re(energy0):
    omegac0 = energy0/aux 

    rtaself_x, rtaself_y, rtaself_z  = green_self_ana3(omegac0,epsi1,epsi3,d_nano,zp)
    
    return np.real(rtaself_x)

def function_ana_xx_im(energy0):
    omegac0 = energy0/aux 

    rtaself_x, rtaself_y, rtaself_z  = green_self_ana(omegac0,epsi1,epsi3,d_nano,zp)
    
    return np.imag(rtaself_x)

def function_ana_xx_im_parallel_num(energy0):
    omegac0 = energy0/aux 

    rtaself_x, rtaself_y, rtaself_z  = green_self_ana_parallel_num(omegac0,epsi1,epsi3,d_nano,zp)
    
    return np.imag(rtaself_x)


def function_ana2_xx_im(energy0):
    omegac0 = energy0/aux 

    rtaself_x, rtaself_y, rtaself_z  = green_self_ana2(omegac0,epsi1,epsi3,d_nano,zp)
    
    return np.imag(rtaself_x)



def function_ana3_xx_im(energy0):
    omegac0 = energy0/aux 

    rtaself_x, rtaself_y, rtaself_z  = green_self_ana3(omegac0,epsi1,epsi3,d_nano,zp)
    
    return np.imag(rtaself_x)

def function_pole_aprox_xx_re(energy0):
    omegac0 = energy0/aux 

    rtaself_x, rtaself_y, rtaself_z  = green_self_pole_aprox(omegac0,epsi1,epsi3,d_nano,zp)
    
    return np.real(rtaself_x)



def function_pole_aprox_xx_im(energy0):
    omegac0 = energy0/aux 

    rtaself_x, rtaself_y, rtaself_z  = green_self_pole_aprox(omegac0,epsi1,epsi3,d_nano,zp)
    
    return np.imag(rtaself_x)




def function_pole_aprox_xx_re_parallel_num(energy0):
    omegac0 = energy0/aux 

    rtaself_x, rtaself_y, rtaself_z  = green_self_pole_aprox_parallel_num(omegac0,epsi1,epsi3,d_nano,zp)
    
    return np.real(rtaself_x)



def function_pole_aprox_xx_im_parallel_num(energy0):
    omegac0 = energy0/aux 

    rtaself_x, rtaself_y, rtaself_z  = green_self_pole_aprox_parallel_num(omegac0,epsi1,epsi3,d_nano,zp)
    
    return np.imag(rtaself_x)


#%%
    
tamfig = [3.5,3]
tamletra = 9
tamtitle  = 9
tamnum = 7
tamlegend = 7
labelpady = 2
labelpadx = 3
pad = 2.5
mk = 1
ms = 3
hp = 0.5
length_marker = 1.5
dpi = 500



def graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad):
    plt.figure(figsize=tamfig)
    plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
    plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
    plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in", pad = pad)
    plt.title(title,fontsize=int(tamtitle*0.9))

    return   

#%%

listy_re_ana = []
listy_re_ana2 = []
listy_re_ana3 = []
listy_re_ana_parallel_num = []
listy_re_num = []
listy_re_pole_aprox = []
listy_re_pole_aprox_parallel_num = []

listy_im_ana = []
listy_im_ana2 = []
listy_im_ana3 = []
listy_im_ana_parallel_num = []
listy_im_num = []
listy_im_pole_aprox = []
listy_im_pole_aprox_parallel_num = []

for value in listx: 

    y_re_ana = function_ana_xx_re(value)       
    y_re_ana2 = function_ana2_xx_re(value)       
    y_re_ana3 = function_ana3_xx_re(value)  
    y_re_ana_parallel_num = function_ana_xx_re_parallel_num(value)       
    y_re_num = function_num_xx_re(value)
    y_re_pole_aprox = function_pole_aprox_xx_re(value)
    y_re_pole_aprox_parallel_num = function_pole_aprox_xx_re_parallel_num(value)
    
    listy_re_ana.append(y_re_ana)
    listy_re_ana2.append(y_re_ana2)
    listy_re_ana3.append(y_re_ana3)
    listy_re_ana_parallel_num.append(y_re_ana_parallel_num)
    listy_re_num.append(y_re_num)
    listy_re_pole_aprox.append(y_re_pole_aprox)
    listy_re_pole_aprox_parallel_num.append(y_re_pole_aprox_parallel_num)

    y_im_ana = function_ana_xx_im(value)       
    y_im_ana2 = function_ana2_xx_im(value)   
    y_im_ana3 = function_ana3_xx_im(value)   
    y_im_ana_parallel_num = function_ana_xx_im_parallel_num(value)   
    y_im_num = function_num_xx_im(value)
    y_im_pole_aprox = function_pole_aprox_xx_im(value)
    y_im_pole_aprox_parallel_num = function_pole_aprox_xx_im_parallel_num(value)

    listy_im_ana.append(y_im_ana)
    listy_im_ana2.append(y_im_ana2)
    listy_im_ana3.append(y_im_ana3)
    listy_im_ana_parallel_num.append(y_im_ana_parallel_num)
    listy_im_num.append(y_im_num)
    listy_im_pole_aprox.append(y_im_pole_aprox)
    listy_im_pole_aprox_parallel_num.append(y_im_pole_aprox_parallel_num)
    
    
#%%
    
label_png = r'$r_{\rm p} = R_p k_\parallel/(k_\parallel - k_{\rm p})$'    
#label2 = r'$r_{\rm p} = k_{\rm p}/(k_\parallel - k_{\rm p})$'      

graph(title,labelx,r'Re{G$_{self}$} ($\mu$m)$^{-3}$',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#plt.plot(listx,listy_re_ana3,'.',ms = ms,color = 'blue',label = 'PP analytical 3')
#plt.plot(listx[0:-2],listy_re_ana[0:-2],'.',ms = ms,color = 'purple',label = 'PP analytical')
plt.plot(listx,listy_re_ana2,'.',ms = ms,color = 'red',label = 'PP ana Exp func')
#plt.plot(listx,listy_re_ana_parallel_num,'.',ms = ms+1,color = 'darkred',label = 'PP ana' + label_png)
plt.plot(listx,listy_re_num,'.',ms = ms,color = 'orange',label = 'full numerical')
plt.plot(listx,listy_re_pole_aprox,'.-',ms = ms,color = 'blue',label = 'PP numerical')
#plt.plot(listx,listy_re_pole_aprox_parallel_num,'.-',ms = 3,color = 'lightseagreen',label = 'PP num ' +label_png)
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0.1,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
#    plt.yscale('log')
os.chdir(path_save)
plt.savefig( 'Re_Gself' + label1 + '.png', format='png')   


graph(title,labelx,r'Im{G$_{self}$} ($\mu$m$^{-3}$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#plt.plot(listx,listy_im_ana,'.',ms = ms,color = 'blue',label = 'PP analytical 3')
#plt.plot(listx[0:-2],listy_im_ana[0:-2],'.',ms = ms,color = 'purple',label = 'PP analytical')
plt.plot(listx,listy_im_ana2,'.',ms = ms,color = 'red',label = 'PP ana Exp func')
#plt.plot(listx,listy_im_ana_parallel_num,'.',ms = ms+1,color = 'darkred',label = 'PP ana' + label_png)
plt.plot(listx,listy_im_num,'.',ms = ms,color = 'orange',label = 'full numerical')
plt.plot(listx,listy_im_pole_aprox,'.-',ms = ms,color = 'blue',label = 'PP numerical')
#plt.plot(listx,listy_im_pole_aprox_parallel_num,'.-',ms = 3,color = 'lightseagreen',label = 'PP num ' +label_png)
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0.1,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
#    plt.yscale('log')
os.chdir(path_save)
plt.savefig( 'Im_Gself' + label1 + '.png', format='png')   



graph(title,labelx,r'Re{G$_{self}$} ($\mu$m)$^{-3}$',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#plt.plot(listx,listy_re_ana3,'.',ms = ms,color = 'blue',label = 'PP analytical 3')
#plt.plot(listx[0:-2],listy_re_ana[0:-2],'.',ms = ms,color = 'purple',label = 'PP analytical')
plt.plot(listx,listy_re_ana2,'.',ms = ms,color = 'red',label = 'PP ana Exp func')
plt.plot(listx,listy_re_ana_parallel_num,'.',ms = ms,color = 'darkred',label = 'PP ana' + label_png)
plt.plot(listx,listy_re_num,'.',ms = ms,color = 'orange',label = 'full numerical')
plt.plot(listx,listy_re_pole_aprox,'.-',ms = ms,color = 'blue',label = 'PP numerical')
plt.plot(listx,listy_re_pole_aprox_parallel_num,'.-',ms = ms,color = 'lightseagreen',label = 'PP num ' +label_png)
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0.1,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
#    plt.yscale('log')
os.chdir(path_save)
plt.savefig( 'Re_Gself_all' + label1 + '.png', format='png',bbox_inches='tight',pad_inches = 0.008,dpi = dpi) 


graph(title,labelx,r'Im{G$_{self}$} ($\mu$m)$^{-3}$',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#plt.plot(listx,listy_im_ana,'.',ms = ms,color = 'blue',label = 'PP analytical 3')
#plt.plot(listx[0:-2],listy_im_ana[0:-2],'.',ms = ms,color = 'purple',label = 'PP analytical')
plt.plot(listx,listy_im_ana2,'.',ms = ms+1,color = 'red',label = 'PP ana Exp func')
plt.plot(listx,listy_im_ana_parallel_num,'.',ms = ms+1,color = 'darkred',label = 'PP ana' + label_png)
plt.plot(listx,listy_im_num,'.',ms = ms,color = 'orange',label = 'full numerical')
plt.plot(listx,listy_im_pole_aprox,'.-',ms = 3,color = 'blue',label = 'PP numerical')
plt.plot(listx,listy_im_pole_aprox_parallel_num,'.-',ms = 3,color = 'lightseagreen',label = 'PP num ' +label_png)
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0.1,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
#    plt.yscale('log')
os.chdir(path_save)
plt.savefig( 'Im_Gself_all' + label1 + '.png', format='png',bbox_inches='tight',pad_inches = 0.008,dpi = dpi)   


tabla = np.array([listx,listy_re_num,listy_im_num])
tabla = np.transpose(tabla)
header1 = 'E [eV]     Re{Gself}     Im{Gself}' + ', ' + title + ', ' + name_this_py
np.savetxt( 'Gself_hBN_num_' + label1 + '.txt', tabla, fmt='%1.11e', delimiter='\t', header = header1)


tabla = np.array([listx,listy_re_ana2,listy_im_ana2])
tabla = np.transpose(tabla)
header1 = 'E [eV]     Re{Gself}     Im{Gself}' + ', ' + title + ', ' + name_this_py
np.savetxt( 'Gself_hBN_ana_' + label1 + '.txt', tabla, fmt='%1.11e', delimiter='\t', header = header1)


tabla = np.array([listx,listy_re_pole_aprox,listy_im_pole_aprox])
tabla = np.transpose(tabla)
header1 = 'E [eV]     Re{Gself}     Im{Gself}' + ', ' + title + ', ' + name_this_py
np.savetxt( 'Gself_hBN_pole_aprox_' + label1 + '.txt', tabla, fmt='%1.11e', delimiter='\t', header = header1)


#%%


# for gm meeting 




graph(title,labelx,r'Re{G$_{self}$} ($\mu$m$^{-3})$',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#plt.plot(listx,listy_re_ana3,'.',ms = ms,color = 'blue',label = 'PP analytical 3')
#plt.plot(listx[0:-2],listy_re_ana[0:-2],'.',ms = ms,color = 'purple',label = 'PP analytical')
#plt.plot(listx,listy_re_ana2,'.',ms = ms,color = 'red',label = 'PP ana Exp func')
plt.plot(listx,listy_re_num,'.',ms = ms,color = 'orange',label = 'full num')
plt.plot(listx,listy_re_pole_aprox_parallel_num,'.-',ms = ms,color = 'lightseagreen',label = 'PP num ' +label_png)
plt.plot(listx,listy_re_ana_parallel_num,'.',ms = ms,color = 'darkred',label = 'PP ana ' + label_png)

#plt.plot(listx,listy_re_pole_aprox,'.-',ms = ms,color = 'blue',label = 'PP numerical')
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0.1,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
#    plt.yscale('log')
os.chdir(path_save)
plt.savefig( 'Re_Gself_GM' + label1 + '.png', format='png',bbox_inches='tight',pad_inches = 0.008,dpi = dpi) 


graph(title,labelx,r'Im{G$_{self}$} ($\mu$m$^{-3}$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#plt.plot(listx,listy_im_ana,'.',ms = ms,color = 'blue',label = 'PP analytical 3')
#plt.plot(listx[0:-2],listy_im_ana[0:-2],'.',ms = ms,color = 'purple',label = 'PP analytical')
#plt.plot(listx,listy_im_ana2,'.',ms = ms+1,color = 'red',label = 'PP ana Exp func')
plt.plot(listx,listy_im_num,'.',ms = ms,color = 'orange',label = 'full num')
plt.plot(listx,listy_im_pole_aprox_parallel_num,'.-',ms = ms,color = 'lightseagreen',label = 'PP num ' +label_png)
plt.plot(listx,listy_im_ana_parallel_num,'.',ms = ms,color = 'darkred',label = 'PP ana' + label_png)
#plt.plot(listx,listy_im_pole_aprox,'.-',ms = 3,color = 'blue',label = 'PP numerical')

plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0.1,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
#    plt.yscale('log')
os.chdir(path_save)
plt.savefig( 'Im_Gself_GM' + label1 + '.png', format='png',bbox_inches='tight',pad_inches = 0.008,dpi = dpi)  









 
label_png2 = r'$r_{\rm p} = R_p k_{\rm p}/(k_\parallel - k_{\rm p})$'    


graph(title,labelx,r'Re{G$_{self}$} ($\mu$m$^{-3}$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#plt.plot(listx,listy_re_ana3,'.',ms = ms,color = 'blue',label = 'PP analytical 3')
#plt.plot(listx[0:-2],listy_re_ana[0:-2],'.',ms = ms,color = 'purple',label = 'PP analytical')
plt.plot(listx,listy_re_num,'.',ms = ms,color = 'orange',label = 'full num')
plt.plot(listx,listy_re_pole_aprox,'.-',ms = ms,color = 'blue',label = 'PP num '  + label_png2)
plt.plot(listx,listy_re_ana2,'.',ms = ms,color = 'red',label = 'PP ana ' + label_png2)
#plt.plot(listx,listy_re_pole_aprox_parallel_num,'.-',ms = ms,color = 'lightseagreen',label = 'PP num ' +label_png)
#plt.plot(listx,listy_re_ana_parallel_num,'.',ms = ms,color = 'darkred',label = 'PP ana ' + label_png)


plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0.1,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
#    plt.yscale('log')
os.chdir(path_save)
plt.savefig( 'Re_Gself_GM_v2' + label1 + '.png', format='png',bbox_inches='tight',pad_inches = 0.008,dpi = dpi) 


graph(title,labelx,r'Im{G$_{self}$} ($\mu$m$^{-3}$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#plt.plot(listx,listy_im_ana,'.',ms = ms,color = 'blue',label = 'PP analytical 3')
#plt.plot(listx[0:-2],listy_im_ana[0:-2],'.',ms = ms,color = 'purple',label = 'PP analytical')
plt.plot(listx,listy_im_num,'.',ms = ms,color = 'orange',label = 'full num')
plt.plot(listx,listy_im_pole_aprox,'.-',ms = ms,color = 'blue',label = 'PP num '+ label_png2)
plt.plot(listx,listy_im_ana2,'.',ms = ms,color = 'red',label = 'PP ana' + label_png2)

#plt.plot(listx,listy_im_pole_aprox_parallel_num,'.-',ms = ms,color = 'lightseagreen',label = 'PP num ' +label_png)
#plt.plot(listx,listy_im_ana_parallel_num,'.',ms = ms,color = 'darkred',label = 'PP ana' + label_png)


plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0.1,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
#    plt.yscale('log')
os.chdir(path_save)
plt.savefig( 'Im_Gself_GM_v2' + label1 + '.png', format='png',bbox_inches='tight',pad_inches = 0.008,dpi = dpi)  