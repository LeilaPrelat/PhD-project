#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
"""

import numpy as np

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


#%%    

def hBn_lambda_p(hbw,epsi1,epsi3):
    """
    Parameters
    ----------
    hbw : energia = hbar*omega en eV
    carrier_density : carrier density
    d_nano : thickness in nm
    decay_rate : 
    masa_eff : 

    Returns
    -------
    devuelve la conductividad de un metal
    """
    epsi_x = epsilon_x(hbw)
    epsi_z = epsilon_z(hbw)

    epsi_m = np.sqrt(epsi_x*epsi_z)
    
    if np.imag(epsi_m) > 0:
        epsi_m = epsi_m
    else:
        epsi_m = -epsi_m
    
    num = (epsi1 - epsi_m)*(epsi3 - epsi_m)
    den = (epsi1 + epsi_m)*(epsi3 + epsi_m)
    
    tot = np.log(num/den)

    cte_epsi = np.sqrt(epsi_x/epsi_z)

########################################################################
### el self image green tensor no funciona bien con este cambio de signo 
### pero sin este cambio, el rp de la pole aprox en el primer polo es negativo 
    
    if np.imag(cte_epsi) > 0: 
        cte_epsi = cte_epsi
    else:
        cte_epsi = -cte_epsi
########################################################################    


    return 4*np.pi*cte_epsi/tot




########################################################################  
########################################################################  
########################################################################  
########################################################################  


def hBn_Rp(hbw,epsi1,epsi3):
    """
    Parameters
    ----------
    hbw : energia = hbar*omega en eV
    carrier_density : carrier density
    d_nano : thickness in nm
    decay_rate : 
    masa_eff : 

    Returns
    -------
    devuelve la conductividad de un metal
    """
    epsi_x = epsilon_x(hbw)
    epsi_z = epsilon_z(hbw)
    
    lambda_p_v = hBn_lambda_p(hbw,epsi1,epsi3) ## lambda/d

    
    
#    c = 3*10**(14)                ### light velocity in micron/seg
#    alfac = 1/137.0359            ### fine structure
    
    num = -epsi_z*epsi1
    den = epsi_x*epsi_z - epsi1**2
    
    tot = num/den

    return tot*lambda_p_v/np.pi




#%%






def hBn_lambda_p_Gself_image(hbw,epsi1,epsi3):
    """
    Parameters
    ----------
    hbw : energia = hbar*omega en eV
    carrier_density : carrier density
    d_nano : thickness in nm
    decay_rate : 
    masa_eff : 

    Returns
    -------
    devuelve la conductividad de un metal
    """
    epsi_x = epsilon_x(hbw)
    epsi_z = epsilon_z(hbw)

    epsi_m = np.sqrt(epsi_x*epsi_z)
    
#    if np.imag(epsi_m) > 0:
#        epsi_m = epsi_m
#    else:
#        epsi_m = -epsi_m
    
    num = (epsi1 - epsi_m)*(epsi3 - epsi_m)
    den = (epsi1 + epsi_m)*(epsi3 + epsi_m)
    
    tot = np.log(num/den)

    cte_epsi = np.sqrt(epsi_x/epsi_z)

########################################################################
### el self image green tensor no funciona bien con este cambio de signo 
### pero sin este cambio, el rp de la pole aprox en el primer polo es negativo 
    
#    if np.imag(cte_epsi) > 0: 
#        cte_epsi = cte_epsi
#    else:
#        cte_epsi = -cte_epsi
########################################################################    


    return 4*np.pi*cte_epsi/tot







def hBn_Rp_Gself_image(hbw,epsi1,epsi3):
    """
    Parameters
    ----------
    hbw : energia = hbar*omega en eV
    carrier_density : carrier density
    d_nano : thickness in nm
    decay_rate : 
    masa_eff : 

    Returns
    -------
    devuelve la conductividad de un metal
    """
    epsi_x = epsilon_x(hbw)
    epsi_z = epsilon_z(hbw)
    
    lambda_p_v = hBn_lambda_p_Gself_image(hbw,epsi1,epsi3) ## lambda/d

    
    
#    c = 3*10**(14)                ### light velocity in micron/seg
#    alfac = 1/137.0359            ### fine structure
    
    num = -epsi_z*epsi1
    den = epsi_x*epsi_z - epsi1**2
    
    tot = num/den

    return tot*lambda_p_v/np.pi










