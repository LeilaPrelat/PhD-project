
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila

Revisar la integral
Ref \cite{Gradshteyn} page 718 6.671 8.
"""
#from scipy import integrate
from mpmath import *
#from scipy import special
import numpy as np
import sys
import os 

mp.dps = 15 # precision del metodo de integracion
mp.pretty = True

graficar = 0

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/External_Efield/fieldE_ref/extra','')
try:
    sys.path.insert(1, path_constants)
    from graphene_sigma import sigma_DL
except ModuleNotFoundError:
    print('graphene_sigma.py no se encuentra en ' + path_basic)

try:
    sys.path.insert(1, path_constants)
    from constants import constantes
except ModuleNotFoundError:
    print('constants.py no se encuentra en ' + path_constants)

pi,hb,c,alfac,mu1,mu2 = constantes()
aux = c*hb

#%%

def revisarG_num(omegac,hbmu,hbgama,int_v):
    E = omegac*aux
    epsi1 = 1
    epsi2 = 1
    mu1 = 1
    n1 = epsi1*mu1
    cte1 = np.sqrt(n1)
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    alfa_p = 1j*(epsi1 + epsi2)/(cond*cte1)
    a = alfa_p
    b = int_v/n1

#    a = 1e5
#    b = 1e-3
    
#    cota_sup2 = 2000
    j0_change = lambda x: j0(a*x)
    Int = lambda x: j0_change(x)*cos(b*x)
    
    periodo1 = 2*pi/b
    
#    periodo2 = besseljzero(0,1)/a : x es real entonces el periodo tambien es real     
#    periodo = np.min([periodo1,periodo2])
    
    
 #   int_final = quadosc(Int, [0, inf], period = periodo1) 
 #   bessel function no es periodica (lo es asintoticamente) pero se pueden usar los ceros
 
    j0zero = lambda n: findroot(j0_change, pi*(n-0.25)/a) 
    int_final = quadosc(Int, [0, inf], zeros = j0zero) 

    return int_final

#%%


def revisarG_ana(omegac,hbmu,hbgama,int_v,Gradshtein_correcto):
    E = omegac*aux
    epsi1 = 1
    epsi2 = 1
    mu1 = 1
    n1 = epsi1*mu1
    cte1 = np.sqrt(n1)
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    alfa_p = 1j*(epsi1 + epsi2)/(cond*cte1)
    a = alfa_p
    b = int_v/n1

#    a = 1e5
 #   b = 1e-3

    if Gradshtein_correcto == 1:
        final = 1/np.sqrt(a**2 - b**2)
    
        if 0 < np.abs(b) < np.abs(a): 
            return final
        else:
            print('caso cero')
            return 0

    else: 
        final = 1/np.sqrt(b**2 - a**2)
    
        if 0 < np.abs(a) <  np.abs(b): 
            return final
        else:
            print('caso cero')
            return 0        
        
#%%

if graficar == 1:
    import matplotlib.pyplot as plt
    import seaborn as sns
    sns.set()   

    tamfig = (4.5,3.5)
    tamlegend = 10
    tamletra = 11
    tamtitle = 10
    tamnum = 9
    labelpady = -1.5
    labelpadx = 0.8
    pad = 0
    mk = 2
    ms = 4
    hp = 0.3

    def graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad):
        plt.figure(figsize=tamfig)
        plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
        plt.ylabel(labely,fontsize=tamletra,labelpad =labelpadx)
        plt.tick_params(labelsize = tamnum, pad = pad)
        plt.title(title,fontsize=int(tamtitle*0.9))

        return   


    path_constants =  path_basic.replace('/External_Efield/fieldE_ref/extra','')       

    int_v = 400
    v = c/int_v

#omega = 0.7*1e12
    cota = 2
    epsi1,epsi2 = 1,1
    hbmu,hbgama = 0.3,0.0001

    n1 = epsi1*mu1
    cte1 = np.sqrt(n1)    
    omegac = 0.1
    E = omegac*aux
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    alfa_p = 1j*(epsi1 + epsi2)/(cond*cte1)
    a = alfa_p
    b = int_v/n1
    
    label1 = '_vs_x'
    title1 = r'$\epsilon_1$ = %i, v = c/%i $\mu$m/s, E = %.2feV' %(epsi1,int_v,E) 
    labelx = r'x'
    title2A = r'$\epsilon_2$ = %i, $\hbar\mu$ = %.2feV, $\hbar\gamma$ = %.4feV' %(epsi2,hbmu,hbgama) 
    title = title1 + '\n'  + title2A   


    def function(x):
        alfa_p = 1j*(epsi1 + epsi2)/(cond*cte1)
        a = alfa_p
        b = int_v/n1

        return j0(a*x)*cos(b*x)


    cota = 4
    list_x = np.linspace(0,cota,41)
    listf_re = []
    listf_im = []
     
    for x in list_x:
        
        rta = function(x)
        listf_re.append(rta.real)
        listf_im.append(rta.imag)

    zeros_re = []
    zeros_im = []
    zeros_de_verdad = []
    periodo1 = 2*pi/b
    
    cant_periodos = int(cota/periodo1)
    for m in range(50):
#        periodo2 = besseljzero(0,m)/a 
#        periodo_re = np.min([periodo1,np.real(periodo2)])
#        periodo_im = np.min([periodo1,np.imag(periodo2)])
#        zeros_re.append(periodo_re)
#        zeros_im.append(periodo_im)
        zeros_de_verdad.append(periodo1*m)                
    
    graph(title,labelx,'function inside integral',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    plt.plot(list_x,listf_re,'.',ms = ms,color = 'blue',label = 'analytical')
    plt.plot(zeros_de_verdad,np.array(np.ones(50))*0,'-',ms = ms,color = 'magenta',label = 'numeric')
    plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=True,handletextpad=hp, handlelength=1)
    plt.tight_layout()
    os.chdir(path_save)
    plt.savefig( 'Re_function_Gradshtein_' + label1 + '.png', format='png')
    

    