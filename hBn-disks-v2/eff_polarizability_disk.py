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

graficar = 0



try:
    sys.path.insert(1, path_basic)
    from constants import constantes
except ModuleNotFoundError:
    print('constants.py no se encuentra en ' + path_basic)

pi,hb,c,alfac,mu1,mu2 = constantes()


############# bands ######################

x1 = 0.09260651629072682 
x2 = 0.10112781954887218


x3 = 0.17030075187969923
x4 = 0.19937343358395992

############# bands ######################

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


def polarizability_parallel(hbw,d_thickness_disk_nano,D_disk_nano,epsilon_silica):
    ""
    a_zeta, b_zeta, c_zeta = -0.01267, -45.34, 0.8635
    a_eta, b_eta, c_eta = 0.03801, -8.569, -0.1108

    
    t = 0.4  ### nm 
    x = t/D_disk_nano
    
    zeta1 = a_zeta*np.exp(b_zeta*x)  + c_zeta    
    eta1 = a_eta*np.exp(b_eta*x)  + c_eta
    

    cte_eta_aux = (epsilon_x(hbw) - 1)/(epsilon_silica(hbw) + 1) 
    eta_parallel = - d_thickness_disk_nano*cte_eta_aux/(2*np.pi*D_disk_nano)
    
    num = zeta1**2
    den = 1/eta_parallel - 1/eta1
    
    D_micro = D_disk_nano*1e-3 ## tiene que estar en 1/micrones³
    D_3 = D_micro**3
    
    
    cte_epsilon = (epsilon_silica(hbw) + 1)/2
    
    return D_3*cte_epsilon*num/den
    

def polarizability_perp(hbw,d_thickness_disk_nano,D_disk_nano,epsilon_silica):
    
    a_zeta, b_zeta, c_zeta = -0.01267, -45.34, 0.8635
    a_eta, b_eta, c_eta = 0.03801, -8.569, -0.1108


    t = 0.4  ### nm
    x = t/D_disk_nano
    
    
    zeta1 = a_zeta*np.exp(b_zeta*x)  + c_zeta
    eta1 = a_eta*np.exp(b_eta*x)  + c_eta
    
    cte_eta_aux = (epsilon_z(hbw) - 1)/(epsilon_silica(hbw) + 1) 
    eta_perp = - d_thickness_disk_nano*cte_eta_aux/(2*np.pi*D_disk_nano)
    
    num = zeta1**2
    den = 1/eta_perp - 1/eta1
    
    D_micro = D_disk_nano*1e-3 ## tiene que estar en 1/micrones³
    D_3 = D_micro**3
    
    cte_epsilon = (epsilon_silica(hbw) + 1)/2
    
    return D_3*cte_epsilon*num/den


#%%

if graficar == 1:

    import matplotlib.pyplot as plt
    import seaborn as sns
    
    sns.set()
    
    name_this_py = os.path.basename(__file__)
    path = os.path.abspath(__file__) #path absoluto del .py actual
    path_basic = path.replace('/' + name_this_py,'')
    path_save = path_basic + '/' + 'polarizalibity_disks'
    
    try:
        sys.path.insert(1, path_basic)
        from constants import constantes
    except ModuleNotFoundError:
        print('constants.py no se encuentra en ' + path_basic)
    
    
    try:
        sys.path.insert(1, path_basic)
        from Silica_epsilon import epsilon_Silica
    except ModuleNotFoundError:
        print('Silica_epsilon.py no se encuentra en ' + path_basic)
    
    pi,hb,c,alfac,mu1,mu2 = constantes()
    
    aux = hb*c
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

    
    colors = ['darkred','steelblue','coral','yellowgreen']
    list_D = [8.8,15,10,20,120]
    
    
    D_nano = 50
    d_nano = 0.1
    
    labelpng = 'd%.2fnm_D%.2fnm' %(d_nano,D_nano)

    
    x1 = 0.09260651629072682 
    x2 = 0.10112781954887218
    x3 = 0.17030075187969923
    x4 = 0.19937343358395992
    
    title = r'$d$ = %.2f nm, D = %.2f nm' %(d_nano,D_nano)
    
    xmin = 0.172
    xmax = 0.197
    list_E = np.linspace(x1,x4,100)

    list_y_perp_re = []
    list_y_parallel_re = []
    
    list_y_perp_im = []
    list_y_parallel_im = []


    for energy in list_E:
        y_perp = polarizability_perp(energy,d_nano,D_nano,epsilon_Silica)
        y_parallel = polarizability_parallel(energy,d_nano,D_nano,epsilon_Silica)


        list_y_perp_re.append(np.real(y_perp))
        list_y_parallel_re.append(np.real(y_parallel))
        list_y_perp_im.append(np.imag(y_perp))
        list_y_parallel_im.append(np.imag(y_parallel))
        
    
    
    graph(title,r'$\hbar\omega$ (eV)',r'Re{$\alpha$} ($\mu$m$^3$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#    plt.plot(listx,listy_re_ana,'.',ms = ms,color = 'purple',label = 'analytical')
    plt.plot(list_E,list_y_perp_re,'.-',ms = ms,color = 'lightseagreen',label = r'Re{$\alpha_{\perp}$}')
    plt.plot(list_E,list_y_parallel_re,'.-',ms = ms,color = 'purple',label = r'Re{$\alpha_{\parallel}$}')
    ejey_aux = np.linspace(np.min([np.min(list_y_perp_re),np.min(list_y_parallel_re)]), np.max([np.max(list_y_perp_re),np.max(list_y_parallel_re)]) , 10)
    for x in [x1,x2,x3]:
        plt.plot(x*np.ones(10), ejey_aux,'--',color = 'grey' )
    
    plt.plot(list_E[-1]*np.ones(10), ejey_aux,'--',color = 'grey' )  
        
    plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0.05,handletextpad=0.2, handlelength=length_marker)
    plt.tight_layout()
    #plt.yscale('log')
    os.chdir(path_save)
    plt.savefig( 'Re_alfa'  + labelpng  + '.png', format='png',bbox_inches='tight',pad_inches = 0.008,dpi = dpi)     
        
    
    graph(title,r'$\hbar\omega$ (eV)',r'Im{$\alpha$} ($\mu$m$^3$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#    plt.plot(listx,listy_re_ana,'.',ms = ms,color = 'purple',label = 'analytical')
    plt.plot(list_E,list_y_perp_im,'.-',ms = ms,color = 'lightseagreen',label = r'Im{$\alpha_{\perp}$}')
    plt.plot(list_E,list_y_parallel_im,'.-',ms = ms ,color = 'purple',label = r'Im{$\alpha_{\parallel}$}')
    ejey_aux = np.linspace(np.min([np.min(list_y_perp_im),np.min(list_y_parallel_im)]), np.max([np.max(list_y_perp_im),np.max(list_y_parallel_im)]) , 10)
    for x in [x1,x2,x3]:
        plt.plot(x*np.ones(10), ejey_aux,'--',color = 'grey' )
    
    plt.plot(list_E[-1]*np.ones(10), ejey_aux,'--',color = 'grey' )  
        

    plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0.05,handletextpad=0.2, handlelength=length_marker)
    plt.tight_layout()
    #plt.yscale('log')
    os.chdir(path_save)
    plt.savefig( 'Im_alfa' + labelpng + '.png', format='png',bbox_inches='tight',pad_inches = 0.008,dpi = dpi)       
    
#%%
