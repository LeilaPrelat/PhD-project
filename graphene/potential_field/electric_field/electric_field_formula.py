#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila

derivadas que necesito para el campo electrico

derivadas del potential electrico sin usar hankel 

"""
from scipy import special
import sys
import os 
import numpy as np

#init_printing(use_unicode=True)

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/potential_field/electric_field','')
#print('Importar modulos necesarios para este codigo')

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

def electric_field_function(omegac,epsi1,epsi2,hbmu,hbgama,x,y,z,xD,yD,zD,zp,px,py,pz):

    E = omegac*aux
#    k0 = omegac #=omega/c
    # x_y = ky/k0
    n1 = epsi1*mu1
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1
    
    cond = 4*pi*alfac*sigma_DL(E,hbmu,hbgama)
    Rp = 2*epsi1/(epsi1 + epsi2)
    alfa_p = 1j*(epsi1 + epsi2)/(cond*cte1)
    kp = alfa_p*k1
    # kp_2 = kp**2 

    R2 = np.sqrt((x + xD)**2 + (y + yD)**2)
    R2_2 = (x + xD)**2 + (y + yD)**2
    R2_4 = ((x + xD)**2 + (y + yD)**2)**2
    
#    R1 = np.sqrt((x - xD)**2 + (y - yD)**2)
    R1_2 = (x - xD)**2 + (y - yD)**2
    R1_4 = ((x - xD)**2 + (y - yD)**2)**2
       
    J0 = special.jn(0,kp*R2)
    J1 = special.jn(1,kp*R2)
    J2 = special.jn(2,kp*R2)
    
    z_dip_barra = k1*(np.abs(z) + 2*zp + np.abs(zD)) 
    exp0 = np.exp(-alfa_p*z_dip_barra)
    
    term0_ref = px*(x + xD)**2 + py*(y + yD)*(x + xD)
    
    term1A = k1*Rp*(kp**3)*0.5*term0_ref/R2_2
    term1B = (J0 - J2)*exp0
    term1 = term1A*term1B
    
    termaux_x = px*(y + yD) - py/(x+xD) + py*(y+yD)**2/(R2_2*(x+xD))
    term2A = k1*Rp*(kp**2)*(y + yD)*termaux_x
    term2B = J1*exp0/R2
    term2 = term2A*term2B
    
    term3A = 3*pz*(x-xD)*np.abs(z-zD)*np.sign(z)
    term3B = (R1_2 + np.abs(z-zD)**2)**(5/2)
    term3 = term3A*term3B
    
    termaux4 = (np.abs(z))**2
    term4A = term0_ref/R2_4
    term4B = -termaux4/((R2_2 + termaux4)**(3/2)) + (R2_2 + termaux4)**(-1/2)
    term4 = term4A*term4B

    term0_dir =  px*(x - xD)**2 + py*(y - yD)*(x - xD)
    
    termaux5 = np.abs(z-zD)**2
    term5A = termaux5/R1_4
    term5B = -termaux5/((R1_2 + termaux5)**(3/2)) + (R1_2 + termaux5)**(-1/2)
    term5 = term5A*term5B
    
    
    term6A = -term0_dir/R1_2
    term6B = 3*(x-xD)*termaux5/((R1_2 + termaux5)**(5/2)) - (x-xD)/((R1_2 + termaux5)**(3/2))
    term6 = term6A*term6B
    
    term7A = -term0_ref/R2_2
    term7B = 3*(x+xD)*termaux4/((R2_2 + termaux4)**(5/2)) - (x+xD)/((R2_2 + termaux4)**(3/2))
    term7 = term7A*term7B
    
    term8A = -termaux4*(y + yD)/((R2_2 + termaux4)**(3/2)) + (R2_2 + termaux4)**(-1/2)
    term8B = px*(y + yD)/R2_2 - py/(x+xD) + py*(y+yD)**2/(R2_2*(x+xD))
    term8 = -term8A*term8B/R2_2
    
    term9A =  -termaux5*(y - yD)/((R1_2 + termaux4)**(3/2)) + (R1_2 + termaux4)**(-1/2)
    term9B = px*(y - yD)/R1_2 - py/(x-xD) + py*(y-yD)**2/(R1_2*(x-xD))
    term9 = -term9A*term9B/R1_2
    
    
    return  term1  + term2 + term3 + term4 + term5 + term6 + term7 + term8 + term9
     
#print_latex(I2_6_aux2)


# para el potential 
#print_latex(term_tot_diff_x)
#print('')
#print_latex(term_tot_diff_z)

#%%


def electric_field_function_y(omegac,epsi1,epsi2,hbmu,hbgama,x,y,z,xD,yD,zD,zp,px,py,pz):

    E = omegac*aux
#    k0 = omegac #=omega/c
    # x_y = ky/k0
    n1 = epsi1*mu1
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1
    
    cond = 4*pi*alfac*sigma_DL(E,hbmu,hbgama)
    Rp = 2*epsi1/(epsi1 + epsi2)
    alfa_p = 1j*(epsi1 + epsi2)/(cond*cte1)
    kp = alfa_p*k1
    # kp_2 = kp**2 

    R2 = np.sqrt((x + xD)**2 + (y + yD)**2)
    R2_2 = (x + xD)**2 + (y + yD)**2
    R2_4 = ((x + xD)**2 + (y + yD)**2)**2
    
    R1 = np.sqrt((x - xD)**2 + (y - yD)**2)
    R1_2 = (x - xD)**2 + (y - yD)**2
    R1_4 = ((x - xD)**2 + (y - yD)**2)**2
       
    J0 = special.jn(0,kp*R2)
    J1 = special.jn(1,kp*R2)
    J2 = special.jn(2,kp*R2)
    
    z_dip_barra = k1*(np.abs(z) + 2*zp + np.abs(zD)) 
    exp0 = np.exp(-alfa_p*z_dip_barra)
    
    term0_ref = py*(y + yD)**2 + px*(y + yD)*(x + xD)
    
    term1A = k1*Rp*(kp**3)*0.5*term0_ref/R2_2
    term1B = (J0 - J2)*exp0
    term1 = term1A*term1B
    
    termaux_x = py*(x + xD) - px/(y + yD) + px*(x + xD)**2/(R2_2*(y + yD))
    term2A = k1*Rp*(kp**2)*(x + xD)*termaux_x
    term2B = J1*exp0/R2
    term2 = term2A*term2B
    
    term3A = 3*pz*(y-yD)*np.abs(z-zD)*np.sign(z)
    term3B = (R1_2 + np.abs(z-zD)**2)**(5/2)
    term3 = term3A*term3B
    
    termaux4 = (np.abs(z))**2
    term4A = term0_ref/R2_4
    term4B = -termaux4/((R2_2 + termaux4)**(3/2)) + (R2_2 + termaux4)**(-1/2)
    term4 = term4A*term4B

    term0_dir =  py*(y - yD)**2 + px*(y - yD)*(x - xD)
    
    termaux5 = np.abs(z-zD)**2
    term5A = termaux5/R1_4
    term5B = -termaux5/((R1_2 + termaux5)**(3/2)) + (R1_2 + termaux5)**(-1/2)
    term5 = term5A*term5B
    
    
    term6A = -term0_dir/R1_2
    term6B = 3*(y-yD)*termaux5/((R1_2 + termaux5)**(5/2)) - (y-yD)/((R1_2 + termaux5)**(3/2))
    term6 = term6A*term6B
    
    term7A = -term0_ref/R2_2
    term7B = 3*(y+yD)*termaux4/((R2_2 + termaux4)**(5/2)) - (y+yD)/((R2_2 + termaux4)**(3/2))
    term7 = term7A*term7B
    
    term8A = -termaux4*(y + yD)/((R2_2 + termaux4)**(3/2)) + (R2_2 + termaux4)**(-1/2)
    term8B = py*(x + xD)/R2_2 - px/(y+yD) + px*(x+xD)**2/(R2_2*(y+yD))
    term8 = -term8A*term8B/R2_2
    
    term9A =  -termaux5*(x - xD)/((R1_2 + termaux4)**(3/2)) + (R1_2 + termaux4)**(-1/2)
    term9B = py*(x - xD)/R1_2 - px/(y-yD) + px*(x-xD)**2/(R1_2*(y-yD))
    term9 = -term9A*term9B/R1_2
    
    
    return  term1  + term2 + term3 + term4 + term5 + term6 + term7 + term8 + term9

#%%

def electric_field_function_z(omegac,epsi1,epsi2,hbmu,hbgama,x,y,z,xD,yD,zD,zp,px,py,pz):

    E = omegac*aux
#    k0 = omegac #=omega/c
    # x_y = ky/k0
    n1 = epsi1*mu1
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1
    
    cond = 4*pi*alfac*sigma_DL(E,hbmu,hbgama)
    Rp = 2*epsi1/(epsi1 + epsi2)
    alfa_p = 1j*(epsi1 + epsi2)/(cond*cte1)
    kp = alfa_p*k1
    # kp_2 = kp**2 

    R2 = np.sqrt((x + xD)**2 + (y + yD)**2)
    R2_2 = (x + xD)**2 + (y + yD)**2
    # R2_4 = ((x + xD)**2 + (y + yD)**2)**2
    
    # R1 = np.sqrt((x - xD)**2 + (y - yD)**2)
    R1_2 = (x - xD)**2 + (y - yD)**2
    # R1_4 = ((x - xD)**2 + (y - yD)**2)**2
       
    # J0 = special.jn(0,kp*R2)
    J1 = special.jn(1,kp*R2)
    # J2 = special.jn(2,kp*R2)
    
    z_dip_barra = k1*(np.abs(z) + 2*zp + np.abs(zD)) 
    exp0 = np.exp(-alfa_p*z_dip_barra)
    
    term0_ref = px*(x + xD)**2 + py*(y + yD)*(x + xD)   
    
    term1A = k1*Rp*(kp**3)*0.5*term0_ref/R2_2
    term1B = J1*exp0
    term1 = -term1A*term1B*np.sign(z)
    
    term2A =  3*pz*(np.abs(z-zD)**2)*np.sign(z)
    term2B = (R1_2 + np.abs(z-zD)**2)**(5/2)
    term2 = term2A/term2B
    
    term3A =  3*pz*np.sign(z-zD)*np.sign(z)
    term3B = (R1_2 + np.abs(z-zD)**2)**(3/2)
    term3 = term3A/term3B    
    

    term0_dir =  py*(y - yD) + px*(x - xD)
    term1_ref = py*(y + yD) + px*(x + xD)
    
    term4A = 3*(np.abs(z-zD)**3)*np.sign(z-zD)
    term4B = (R1_2 + np.abs(z-zD)**2)**(5/2)
    term4C = 3*np.abs(z-zD)*np.sign(z-zD)
    term4D = (R1_2 + np.abs(z-zD)**2)**(3/2)
    term4 = -term0_dir*(term4A/term4B -  term4C/term4D)/R1_2
    
    
    term5A = 3*(np.abs(z) + np.abs(zD))**3*np.sign(z)
    term5B = (R2_2 + (np.abs(z) + np.abs(zD))**2)**(5/2)
    term5C = 3*(np.abs(z) + np.abs(zD))*np.sign(z)
    term5D = (R2_2 + (np.abs(z) + np.abs(zD))**2)**(3/2)
    term5 = -term1_ref*(term5A/term5B -  term5C/term5D)/R2_2
    
    return term1 + term2 + term3 + term4 + term5


#%%


def electric_field_function_tot(omegac,epsi1,epsi2,hbmu,hbgama,x,y,z,xD,yD,zD,zp,px,py,pz):
    
    Ez = electric_field_function_z(omegac,epsi1,epsi2,hbmu,hbgama,x,y,z,xD,yD,zD,zp,px,py,pz)
        
    Ex = electric_field_function(omegac,epsi1,epsi2,hbmu,hbgama,x,y,z,xD,yD,zD,zp,px,py,pz)
    
    Ey = electric_field_function_y(omegac,epsi1,epsi2,hbmu,hbgama,x,y,z,xD,yD,zD,zp,px,py,pz)
        
    return np.abs(Ez)**2 + np.abs(Ex)**2 + np.abs(Ey)**2

#%%
