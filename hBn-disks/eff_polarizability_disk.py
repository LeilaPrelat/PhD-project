#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
conductividad del grafeno  
ref:
@BOOK{kubo2,
   author       = {R. A. Depine}, 
   year         = {2017},
   title        = {Graphene Optics: Electromagnetic solution of canonical problems}, 
   publisher    = {IOP Concise Physics.\ San Rafael, CA, USA: Morgan and Claypool Publishers}
}

"""
import numpy as np
import sys
import os 


name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')

graficar = 1

try:
    sys.path.insert(1, path_basic)
    from green_self_image import green_self_num, green_self_pole_aprox, green_self_ana2
except ModuleNotFoundError:
    print('green_self_image.py no se encuentra en ' + path_basic)


try:
    sys.path.insert(1, path_basic)
    from constants import constantes
except ModuleNotFoundError:
    print('constants.py no se encuentra en ' + path_basic)

pi,hb,c,alfac,mu1,mu2 = constantes()


#%%

def epsilon_x(hbw):
    
    epsi_inf = 4.87
    hbgamma = 0.87*1e-3
    f_x = 1.83
    
    num = (170.1*1e-3)**2
    
    den = hbw*(hbw + 1j*hbgamma ) - num
    
    return epsi_inf - f_x*(num/den)


def epsilon_z(hbw):
    
    epsi_inf = 2.95
    hbgamma = 0.25*1e-3
    f_x = 0.61
    
    num = (92.5*1e-3)**2
    
    den = hbw*(hbw + 1j*hbgamma ) - num
    
    return epsi_inf - f_x*(num/den)


def polarizability_parallel(hbw,D_nano,epsilon):
    
    a_zeta, b_zeta, c_zeta = -0.01267, -45.34, 0.8635
    a_eta, b_eta, c_eta = 0.03801, -8.569, -0.1108

    D = D_nano*1e-3 ## tiene que estar en 1/micrones³
    
    t = D
    x = t/D
    
    zeta1 = a_zeta*np.exp(b_zeta*x)  + c_zeta    
    eta1 = a_eta*np.exp(b_eta*x)  + c_eta
    
#    omegac = hbw/(c*hb)
#    
    eta_parallel = 1j*np.imag(epsilon_x(hbw))/(D*epsilon)
    
    num = zeta1**2
    den = 1/eta_parallel - 1/eta1
    
    D_3 = D**3
    
    return D_3*epsilon*num/den
    

def polarizability_perp(hbw,D_nano,epsilon):
    
    a_zeta, b_zeta, c_zeta = -0.01267, -45.34, 0.8635
    a_eta, b_eta, c_eta = 0.03801, -8.569, -0.1108

    D = D_nano*1e-3 ## tiene que estar en 1/micrones³
   
    t = D
    x = t/D
    
    zeta1 = a_zeta*np.exp(b_zeta*x)  + c_zeta
    eta1 = a_eta*np.exp(b_eta*x)  + c_eta
    
#    omegac = hbw/(c*hb)
#    
    eta_perp = 1j*np.imag(epsilon_x(hbw))/(D*epsilon)
    
    num = zeta1**2
    den = 1/eta_perp - 1/eta1
    
    D_3 = D**3
    
    return D_3*epsilon*num/den




#%%

def eff_polarizability_ana(hbw,epsi1,epsi3,d_nano,zp,D,epsilon):
    
    alffa_x = polarizability_parallel(hbw,D,epsilon)
    alffa_y = polarizability_parallel(hbw,D,epsilon)
    alffa_z = polarizability_perp(hbw,D,epsilon)
  
    omegac = hbw/(hb*c)
    
    rtaself_x, rtaself_y, rtaself_z  =  green_self_ana2(omegac,epsi1,epsi3,d_nano,zp)
    alffa_eff = (1/alffa_x -  rtaself_x + 1/alffa_y -  rtaself_y +  1/alffa_z -  rtaself_z)**(-1)
    
    return alffa_eff

def eff_polarizability_num(hbw,epsi1,epsi3,d_nano,zp,D,epsilon):
    
    alffa_x = polarizability_parallel(hbw,D,epsilon)
    alffa_y = polarizability_parallel(hbw,D,epsilon)
    alffa_z = polarizability_perp(hbw,D,epsilon)
  
    omegac = hbw/(hb*c)
    
    rtaself_x, rtaself_y, rtaself_z  =  green_self_num(omegac,epsi1,epsi3,d_nano,zp)
    alffa_eff = (1/alffa_x -  rtaself_x + 1/alffa_y -  rtaself_y +  1/alffa_z -  rtaself_z)**(-1)
    
    return alffa_eff


def eff_polarizability_pole_aprox(hbw,epsi1,epsi3,d_nano,zp,D,epsilon):
    
    alffa_x = polarizability_parallel(hbw,D,epsilon)
    alffa_y = polarizability_parallel(hbw,D,epsilon)
    alffa_z = polarizability_perp(hbw,D,epsilon)
  
    omegac = hbw/(hb*c)
    
    rtaself_x, rtaself_y, rtaself_z  =  green_self_pole_aprox(omegac,epsi1,epsi3,d_nano,zp)
    alffa_eff = (1/alffa_x -  rtaself_x + 1/alffa_y -  rtaself_y +  1/alffa_z -  rtaself_z)**(-1)

    return alffa_eff


#%%

if graficar == 1:

    import matplotlib.pyplot as plt
    import seaborn as sns
    
    sns.set()
    
    name_this_py = os.path.basename(__file__)
    path = os.path.abspath(__file__) #path absoluto del .py actual
    path_basic = path.replace('/' + name_this_py,'')
    path_save = path_basic + '/' + 'eff_polarizalibity_disks'
    
    try:
        sys.path.insert(1, path_basic)
        from constants import constantes
    except ModuleNotFoundError:
        print('constants.py no se encuentra en ' + path_basic)
    
    pi,hb,c,alfac,mu1,mu2 = constantes()
    
    aux = hb*c
    tamfig = (4.5,3.5)
    tamlegend = 10
    tamletra = 10
    tamtitle = 10
    tamnum = 9
    loc2 = [0,1]
    pad = -2
    labelpady = -1.5
    labelpadx = 2
    lw = 1.5
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

    
    colors = ['darkred','steelblue','coral','yellowgreen']
    list_D = [8.8,15,10,20,120]
    
    epsi1,epsi3 = 1,1
    D = list_D[-1]
    
    d_nano = 10 
    zp = 0.05
    epsilon = 1
    
    x1 = 0.09260651629072682 
    x2 = 0.10112781954887218
    x3 = 0.17030075187969923
    x4 = 0.19937343358395992
    
    title = r'$\epsilon_1$ = %i, $\epsilon_3$ = %i, D = %i nm, d = %i nm, $z_p$ = %i nm' %(epsi1,epsi3,D,d_nano,zp*1e3)
    
    xmin = 0.172
    xmax = 0.197
    list_E = np.linspace(xmin,xmax,100)

    list_y_num_re = []
    list_y_ana_re = []
    list_y_pole_aprox_re = []

    list_y_num_im = []
    list_y_ana_im = []
    list_y_pole_aprox_im = []

    
    for energy in list_E:
        y_num = eff_polarizability_num(energy,epsi1,epsi3,d_nano,zp,D,epsilon)
        y_ana = eff_polarizability_ana(energy,epsi1,epsi3,d_nano,zp,D,epsilon)
        y_pole_aprox = eff_polarizability_pole_aprox(energy,epsi1,epsi3,d_nano,zp,D,epsilon)

        list_y_num_re.append(np.real(y_num))
        list_y_ana_re.append(np.real(y_ana))
        list_y_pole_aprox_re.append(np.real(y_pole_aprox))
    
    
        list_y_num_im.append(np.imag(y_num))
        list_y_ana_im.append(np.imag(y_ana))
        list_y_pole_aprox_im.append(np.imag(y_pole_aprox))
    
    
    
    graph(title,'$\hbar\omega$ [eV]',r'Re{$\alpha_{eff}$}',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#    plt.plot(listx,listy_re_ana,'.',ms = ms,color = 'purple',label = 'analytical')
    plt.plot(list_E,list_y_num_re,'.-',ms = ms,color = 'lightseagreen',label = 'full numerical')
    plt.plot(list_E,list_y_ana_re,'.',ms = ms+ 1,color = 'purple',label = 'analytical')
    plt.plot(list_E,list_y_pole_aprox_re,'.-',ms = ms,color = 'darkred',label = 'PP numerical')

    plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0.05,handletextpad=0.2, handlelength=length_marker)
    plt.tight_layout()
    #plt.yscale('log')
    os.chdir(path_save)
    plt.savefig( 'Re_alfa_eff_d%inm' %(d_nano)  + '.png', format='png')   
        
    
    graph(title,'$\hbar\omega$ [eV]',r'Im{$\alpha_{eff}$}',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#    plt.plot(listx,listy_re_ana,'.',ms = ms,color = 'purple',label = 'analytical')
    plt.plot(list_E,list_y_num_im,'.-',ms = ms,color = 'lightseagreen',label = 'full numerical')
    plt.plot(list_E,list_y_ana_im,'.',ms = ms + 1,color = 'purple',label = 'analytical')
    plt.plot(list_E,list_y_pole_aprox_im,'.-',ms = ms,color = 'darkred',label = 'PP numerical')

    plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0.05,handletextpad=0.2, handlelength=length_marker)
    plt.tight_layout()
    #plt.yscale('log')
    os.chdir(path_save)
    plt.savefig( 'Im_alfa_eff_d%inm' %(d_nano) + '.png', format='png')       
    
#%%
