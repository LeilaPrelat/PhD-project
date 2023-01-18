#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila

derivadas que necesito para el campo electrico

derivadas del potential electrico sin usar hankel 

"""
from sympy import *
import sys
import os 

#init_printing(use_unicode=True)

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/potential_field','')
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

def electric_field_function(omegac,epsi1,epsi2,hbmu,hbgama,x,y,z,xD,yD,zD,b,zp,px,py,pz):

    E = omegac*aux
#    k0 = omegac #=omega/c
    # x_y = ky/k0
    n1 = epsi1*mu1
    cte1 = sqrt(n1)
    k1 = omegac*cte1
    
    cond = 4*pi*alfac*sigma_DL(E,hbmu,hbgama)
    Rp = 2*epsi1/(epsi1 + epsi2)
    alfa_p = 1j*(epsi1 + epsi2)/(cond*cte1)
    kp = alfa_p*k1
    kp_2 = kp**2 

    x_var, y_var, z_var = symbols('x_var y_var z_var')
    
    R = sqrt(x_var**2 + y_var**2)
    phi = atan(y_var/x_var)
    J1 = besselj(1,kp*R)
    J0 = besselj(0,kp*R)
    expo_f = exp(-kp*(2*zp + Abs(z_var) + Abs(b)))
    
    term0_ref = (Abs(z_var) + Abs(b))**2
    term1_ref = term0_ref*((term0_ref + R**2)**(-3/2)) - (term0_ref + R**2)**(-1/2)
    
    term0_dir = Abs(z_var - zD)**2
    term1_dir = term0_dir*((term0_dir + R**2)**(-3/2)) - (term0_dir + R**2)**(-1/2)
    term_px_py = px*cos(phi) + py*sin(phi)
    
    term1_final = -(term1_dir + term1_ref)*term_px_py/R 
    
    term2_final = -Rp*kp_2*term_px_py*J1*expo_f
    
    term3_final = sign(z)*pz*(term0_ref*((term0_ref + R**2)**(-3/2)) + term0_dir*((term0_dir + R**2)**(-3/2)))
    term4_final = sign(z)*pz*(Rp*kp_2*J0*expo_f)
    
    term_tot_final = term1_final + term2_final + term3_final + term4_final 
    
    Efield_x = -diff( term_tot_final, x_var)
    Efield_y = -diff( term_tot_final, y_var)
    Efield_z = -diff( term_tot_final, z_var)
    
    rta1_Efield_x = Efield_x.subs(x_var,x)
    rta2_Efield_x = rta1_Efield_x.subs(y_var,y)
    rta3_Efield_x = rta2_Efield_x.subs(z_var,z)
    
    rta1_Efield_y = Efield_y.subs(x_var,x)
    rta2_Efield_y = rta1_Efield_y.subs(y_var,y)
    rta3_Efield_y = rta2_Efield_y.subs(z_var,z)    
    
    rta1_Efield_z = Efield_z.subs(x_var,x)
    rta2_Efield_z = rta1_Efield_z.subs(y_var,y)
    rta3_Efield_z = rta2_Efield_z.subs(z_var,z)    
    
    rta_final_final = Abs(rta3_Efield_x)**2 +  Abs(rta3_Efield_y)**2 +  Abs(rta3_Efield_z)**2  
    
    return rta_final_final 
    
#print_latex(I2_6_aux2)


# para el potential 
#print_latex(term_tot_diff_x)
#print('')
#print_latex(term_tot_diff_z)

