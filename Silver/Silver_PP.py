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

#%%

def epsilon_m(hbw):
    
    epsi_b = 4
    hbgamma = 21*1e-3
    
    
    num = (9.17)**2
    
    den = hbw*(hbw + 1j*hbgamma ) 
    
    return epsi_b - num/den

#%%    

def Silver_lambda_p(hbw,epsi1,epsi3):
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
    
    epsi_m = epsilon_m(hbw)
    

    
    num = (epsi1 - epsi_m)*(epsi3 - epsi_m)
    den = (epsi1 + epsi_m)*(epsi3 + epsi_m)
    
    tot = np.log(num/den)


    return 4*np.pi/tot ##lambda dividido d


def Silver_Rp(hbw,epsi1,epsi3):
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

    
    lambda_p_v = Silver_lambda_p(hbw,epsi1,epsi3)

    epsi_m = epsilon_m(hbw)
    
#    c = 3*10**(14)                ### light velocity in micron/seg
#    alfac = 1/137.0359            ### fine structure
    
    num = -epsi_m*epsi1
    den = epsi_m**2 - epsi1**2
    
    tot = num/den

    return tot*lambda_p_v/np.pi



#%%

