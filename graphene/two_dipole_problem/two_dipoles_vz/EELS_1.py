#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila
"""
import numpy as np
import sys
import os 

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/two_dipoles_vx','')
#print('Importar modulos necesarios para este codigo')

try:
    sys.path.insert(1, path_constants)
    from dipole_moment_z import dipole_moment_z_num,dipole_moment_z_ana
except ModuleNotFoundError:
    print('dipole_moment_z.py no se encuentra en ' + path_basic)


try:
    sys.path.insert(1, path_constants)
    from fieldE_ref_z import Efield_ref_ana,Efield_ref_fresnel
except ModuleNotFoundError:
    print('fieldE_ref.py no se encuentra en ' + path_basic)
    


try:
    sys.path.insert(1, path_constants)
    from fieldE_direct_z import Efield_ANA
except ModuleNotFoundError:
    print('fieldE_direct_z.py no se encuentra en ' + path_basic)
      


try:
    sys.path.insert(1, path_constants)
    from constants import constantes
except ModuleNotFoundError:
    print('constants.py no se encuentra en ' + path_constants)

pi,hb,c,alfac,mu1,mu2 = constantes()
aux = c*hb


#%%

def EELS_1_ana(omegac,epsi1,epsi2,hbmu,hbgama,int_v,zp,x1,z1,x2,z2,omega01,kappa_factor_omega01,kappa_r_factor1,omega02,kappa_factor_omega02,kappa_r_factor2):     
    """    
    Parameters
    ----------
    omegac : omega/c = k0 en 1/micrometros    
    epsi1 : epsilon del medio de arriba del plano
    epsi2 : epsilon del medio de abajo del plano
    hbmu : chemical potential in eV  
    hbgama : collision frequency in eV
    z : coordenada z
    xD : coordenada x del dipolo 
    yD : coordenada y del dipolo
    zD : coordenada z del dipolo 
    zp : posicion del plano (>0)
    px : coordenada x del dipolo 
    py : coordenada y del dipolo
    pz : coordenada z del dipolo
    Returns
    -------
    formula del potencial electric con QE approximation, rp con 
    aproximacion del polo y con aprox de principal value para las integrales
    con rp
    """
    E = omegac*aux
    E_meV = E*1e3
    
    p1 = dipole_moment_z_ana(omegac,epsi1,epsi2,hbmu,hbgama,int_v,zp,x1,z1,x2,z2,omega01,kappa_factor_omega01,kappa_r_factor1,omega02,kappa_factor_omega02,kappa_r_factor2)
    
    Eref1 = Efield_ref_ana(omegac,epsi1,epsi2,hbmu,hbgama,x2,z2,x1,zp,int_v)

    
    R1 = np.abs(x1)
    Edir1 = Efield_ANA(omegac,epsi1,R1,z1,int_v)

    rta = 1j*p1*(Eref1 + Edir1)

    print(p1)
    
    print('Edir:', -1j*p1*(Edir1))
    print('Eref:', -1j*p1*(Eref1))
    
    cte_final = alfac*c ## aparece e^2/hbar  (aparece e dentro de los campos)
    
    
    cte_que_no_puse_en_los_campos = (omegac*c)**2
#    print('1:', 1j*p1)
#    print('2:', (Eref1 + Edir1))
#    
    expression_on_overleaf = rta*cte_final/(np.pi*cte_que_no_puse_en_los_campos) #queda en unidad de segundos
    
    #queremos en unidad de energia 
    
    final = expression_on_overleaf*E_meV*omegac*c
    
    return rta*E_meV*alfac/(np.pi*omegac) 


#%%

def EELS_1_num(omegac,epsi1,epsi2,hbmu,hbgama,int_v,zp,x1,z1,x2,z2,omega01,kappa_factor_omega01,kappa_r_factor1,omega02,kappa_factor_omega02,kappa_r_factor2):     
    """    
    Parameters
    ----------
    omegac : omega/c = k0 en 1/micrometros    
    epsi1 : epsilon del medio de arriba del plano
    epsi2 : epsilon del medio de abajo del plano
    hbmu : chemical potential in eV  
    hbgama : collision frequency in eV
    z : coordenada z
    xD : coordenada x del dipolo 
    yD : coordenada y del dipolo
    zD : coordenada z del dipolo 
    zp : posicion del plano (>0)
    px : coordenada x del dipolo 
    py : coordenada y del dipolo
    pz : coordenada z del dipolo
    Returns
    -------
    formula del potencial electric con QE approximation, rp con 
    aproximacion del polo y con aprox de principal value para las integrales
    con rp
    """
    E = omegac*aux
    E_meV = E*1e3
    
    p1 = dipole_moment_z_num(omegac,epsi1,epsi2,hbmu,hbgama,int_v,zp,x1,z1,x2,z2,omega01,kappa_factor_omega01,kappa_r_factor1,omega02,kappa_factor_omega02,kappa_r_factor2)
    
    Eref1 = Efield_ref_fresnel(omegac,epsi1,epsi2,hbmu,hbgama,x2,z2,x1,zp,int_v)
    R1 = np.abs(x1)
    Edir1 = Efield_ANA(omegac,epsi1,R1,z1,int_v)

    rta = 1j*p1*(Eref1 + Edir1)
    
    

    
    cte_final = alfac*c ## aparece e^2/hbar  (aparece e dentro de los campos)

 
    
    cte_que_no_puse_en_los_campos = (omegac*c)**2
#    print('1:', 1j*p1)
#    print('2:', (Eref1 + Edir1))
#    
    expression_on_overleaf = rta*cte_final/(np.pi*cte_que_no_puse_en_los_campos) #queda en unidad de segundos
    
    #queremos en unidad de energia 
    
    final = expression_on_overleaf*E_meV*omegac*c
    
    return rta*E_meV*alfac/(np.pi*omegac) 



#%%
